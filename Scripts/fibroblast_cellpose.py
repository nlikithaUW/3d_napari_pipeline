"""Count endothelial cells from the VE-cadherin channel with Cellpose, and
derive per-image fibroblast counts as (total nuclei - endothelial cells).

Why Cellpose: VE-cadherin is a junctional/membrane marker and is speckly here,
which breaks threshold+watershed segmentation. Cellpose's cyto model is shape-
aware and robust to that, and is the standard tool for counting cells from a
membrane/junctional stain.

For each OME-TIFF in --dir it:
  1. reads the VE-Cadherin channel plane-by-plane (skips any corrupt Z-slice),
  2. segments endothelial cells with Cellpose (2D per-Z stitched into 3D by
     default, or true 3D with --do-3d),
  3. counts cells, pulls total nuclei from <stem>_saved/fibroblast_count.npz,
  4. writes a row to endothelial_counts.csv and a QC overlay PNG per image.

Setup (in the image_analysis env):
    pip install cellpose          # torch+CUDA already present for your GPU

Run (from the 3d_napari_pipeline folder):
    python Scripts/fibroblast_cellpose.py --dir "C:\\...\\images\\WT"

Smoke-test one image first and eyeball its QC overlay before trusting the batch:
    python Scripts/fibroblast_cellpose.py --dir "...\\WT" --only 006

Notes:
  * --diameter defaults to auto (Cellpose estimates the endothelial cell size).
    If auto looks wrong on the QC overlay, set it explicitly in microns-in-pixels
    (it's in *pixels*; ~ endothelial diameter / xy_um_per_px).
  * --ds 2 downsamples XY 2x for speed/memory (recommended; 3D at full res is
    heavy). Counts are scale-invariant to modest downsampling.
  * --do-3d runs true volumetric Cellpose (slower, more memory) instead of
    2D+stitch. Stitch is the robust default for a thin monolayer.
"""

from __future__ import annotations

import os

# Windows PyTorch/MKL ship their own OpenMP runtimes; when both load you get
# "OMP: Error #15 ... libiomp5md.dll already initialized" and the process aborts.
# Allow the duplicate (the documented workaround) and cap threads for stability.
# Must be set BEFORE numpy/torch/cellpose import.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import argparse
import csv
import glob
import json
import re
from pathlib import Path

import numpy as np
import tifffile

DEFAULT_DIR = (Path(__file__).resolve().parent.parent / "3d_images_folder" /
               "images" / "WT")


def _channel_index(xml: str, want: str) -> int | None:
    names = re.findall(r'<Channel[^/>]*\sName="([^"]+)"', xml)
    if not names:
        return None
    for i, n in enumerate(names):
        if n == want:
            return i
    for i, n in enumerate(names):
        if "cad" in n.lower():
            return i
    return None


def _read_channel(tif: Path, channel: str, ds: int):
    """Read one channel as a (Z, Y, X) float32 volume, plane-by-plane."""
    with tifffile.TiffFile(str(tif)) as tf:
        s = tf.series[0]
        xml = tf.ome_metadata or ""
        cidx = _channel_index(xml, channel)
        if cidx is None:
            return None
        axes = s.axes
        if axes == "ZCYX":
            Z, C = s.shape[0], s.shape[1]
            page_of = lambda z: z * C + cidx
        elif axes == "CZYX":
            C, Z = s.shape[0], s.shape[1]
            page_of = lambda z: cidx * Z + z
        else:
            raise ValueError(f"{tif.name}: unsupported axes {axes}")
        H, W = s.shape[-2], s.shape[-1]
        planes = []
        for z in range(Z):
            try:
                planes.append(tf.pages[page_of(z)].asarray()[::ds, ::ds])
            except Exception:
                planes.append(np.zeros((H // ds + (H % ds > 0),
                                        W // ds + (W % ds > 0)), np.uint16))
    return np.stack(planes).astype(np.float32)


def _total_nuclei(saved: Path) -> int | None:
    fc = saved / "fibroblast_count.npz"
    if fc.exists():
        try:
            return int(np.load(fc)["total"])
        except Exception:
            return None
    return None


def _qc_overlay(vol, masks, out_png, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries
    z = vol.shape[0] // 2
    img = vol[z]
    m = masks[z] if masks.ndim == 3 else masks
    lo, hi = np.percentile(img, [1, 99.5])
    disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.4))
    ax[0].imshow(disp, cmap="gray"); ax[0].set_title(f"VE-cad Z-slice {z}")
    ax[0].axis("off")
    ax[1].imshow(disp, cmap="gray")
    b = find_boundaries(m, mode="outer")
    ax[1].imshow(np.ma.masked_where(~b, b), cmap="autumn", alpha=0.8)
    ax[1].set_title(f"Cellpose cells on this slice: {int(m.max())}")
    ax[1].axis("off")
    fig.suptitle(title)
    fig.tight_layout(); fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _tissue_mask(vol):
    """3D mask of real tissue (endothelial walls) vs background noise, so noise
    can't be mistaken for a wall or amplified by CLAHE. Otsu on a smoothed copy,
    then keep the largest connected component (drops scattered speckle)."""
    from scipy.ndimage import gaussian_filter, label, binary_closing
    from skimage.filters import threshold_otsu
    v = gaussian_filter(vol, (1.0, 2.0, 2.0))
    m = v > threshold_otsu(v)
    m = binary_closing(m, iterations=1)
    lbl, n = label(m)
    if n > 1:
        sizes = np.bincount(lbl.ravel()); sizes[0] = 0
        m = lbl == int(sizes.argmax())
    return m


def _split_z(vol, mask=None):
    """Depth that separates the tube's two walls: the valley in the per-Z tissue
    profile, taken *within the tissue depth range* so a tube sitting high (or
    low) in the stack still splits between its two walls rather than in the empty
    part of the volume."""
    from scipy.ndimage import gaussian_filter1d
    if mask is not None:
        prof = mask.reshape(mask.shape[0], -1).sum(1).astype(float)
    else:
        prof = vol.reshape(vol.shape[0], -1).sum(1).astype(float)
    if prof.max() <= 0:
        return vol.shape[0] // 2
    ps = gaussian_filter1d(prof, 1.5)
    zr = np.where(ps > 0.10 * ps.max())[0]          # tissue-containing depths
    if len(zr) < 3:
        return int(round(zr.mean())) if len(zr) else len(ps) // 2
    z0, z1 = int(zr[0]), int(zr[-1])
    return z0 + int(np.argmin(ps[z0:z1 + 1]))       # valley between the 2 walls


def _prep_masked(proj, foot, clahe=True):
    """Normalize a projection using only tissue pixels; background stays 0 so
    CLAHE never amplifies noise."""
    from skimage.exposure import equalize_adapthist
    out = np.zeros_like(proj, dtype=np.float32)
    px = proj[foot]
    if px.size == 0:
        return out
    lo, hi = np.percentile(px, [1, 99.5])
    n = np.clip((proj - lo) / (hi - lo + 1e-9), 0, 1).astype(np.float32)
    n[~foot] = 0
    if clahe:
        n = equalize_adapthist(n, clip_limit=0.02).astype(np.float32)
        n[~foot] = 0
    return n


def _seg2d(model, img2d, newapi, diam, cellprob=0.0):
    kw = dict(diameter=diam, channel_axis=None, cellprob_threshold=cellprob)
    if not newapi:
        kw.update(channels=[0, 0])
    return model.eval(img2d, **kw)[0]


def _preprocess(img, denoise_r=1, bg_sigma=20.0, clahe=True):
    """Prep an en-face projection for Cellpose so faint/noisy walls segment as
    well as bright ones:
      1. median denoise  -> kills speckle without blurring junction lines
      2. background subtraction (subtract a heavily-blurred copy) -> removes the
         diffuse haze that swamps the far wall, flattens uneven illumination
      3. percentile normalize + CLAHE -> local contrast so thin junctions pop
    """
    from scipy.ndimage import median_filter, gaussian_filter
    from skimage.exposure import equalize_adapthist
    x = img.astype(np.float32)
    if denoise_r > 0:
        x = median_filter(x, size=2 * denoise_r + 1)
    if bg_sigma > 0:
        x = np.clip(x - gaussian_filter(x, bg_sigma), 0, None)
    lo, hi = np.percentile(x, [1, 99.5])
    n = np.clip((x - lo) / (hi - lo + 1e-9), 0, 1)
    return equalize_adapthist(n, clip_limit=0.02) if clahe else n


def _count_masks(m, min_area):
    if min_area > 0 and m.max() > 0:
        ids, cnt = np.unique(m[m > 0], return_counts=True)
        for i in ids[cnt < min_area]:
            m[m == i] = 0
    return int(len(np.unique(m[m > 0]))), m


def _enface_count(model, vol, newapi, diam, min_area, invert, cellprob=0.0,
                  clahe=True, denoise_r=1, bg_sigma=20.0, tissue=True):
    """Split into top/bottom walls, project each en-face, segment 2D, sum.

    With tissue=True (default): mask real tissue vs background noise, split at
    the valley within the tissue depth range, and normalize per wall over tissue
    pixels only -- robust to branch points, off-center vessels, and tubes that
    sit high/low in the stack (where a single global Z-split fed noise to one
    wall). tissue=False falls back to the plain split + _preprocess."""
    if tissue:
        mask = _tissue_mask(vol)
        mvol = np.where(mask, vol, 0).astype(np.float32)
        zs = _split_z(vol, mask)
        bottom = mvol[:zs].max(0) if zs > 0 else mvol[0]
        top = mvol[zs:].max(0) if zs < vol.shape[0] else mvol[-1]
        foot_top = mask[zs:].any(0) if zs < vol.shape[0] else mask[-1]
        foot_bot = mask[:zs].any(0) if zs > 0 else mask[0]
        imgs = [_prep_masked(top, foot_top, clahe),
                _prep_masked(bottom, foot_bot, clahe)]
    else:
        zs = _split_z(vol)
        bottom = vol[:zs].max(0) if zs > 0 else vol[0]
        top = vol[zs:].max(0) if zs < vol.shape[0] else vol[-1]
        imgs = [_preprocess(top, denoise_r, bg_sigma, clahe),
                _preprocess(bottom, denoise_r, bg_sigma, clahe)]
    if invert:
        imgs = [im.max() - im for im in imgs]
    masks = [_seg2d(model, im, newapi, diam, cellprob) for im in imgs]
    counts, masks = zip(*[_count_masks(m, min_area) for m in masks])
    return sum(counts), list(counts), [imgs[0], imgs[1]], list(masks), zs


def _enface_qc(projs, masks, counts, out_png, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    for a, img, m, c, nm in zip(ax, projs, masks, counts, ("top wall", "bottom wall")):
        lo, hi = np.percentile(img, [1, 99.5])
        a.imshow(np.clip((img - lo) / (hi - lo + 1e-9), 0, 1), cmap="gray")
        b = find_boundaries(m, mode="outer")
        a.imshow(np.ma.masked_where(~b, b), cmap="autumn", alpha=0.8)
        a.set_title(f"{nm} en-face: {c} cells"); a.axis("off")
    fig.suptitle(title)
    fig.tight_layout(); fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--pattern", default="*.tif")
    ap.add_argument("--channel", default="VE-Cadherin")
    ap.add_argument("--model", default="cyto3",
                    help="Cellpose model (cyto3/cyto2/cyto). cyto3 is default.")
    ap.add_argument("--diameter", type=float, default=0.0,
                    help="cell diameter in PIXELS; 0 = Cellpose auto-estimate")
    ap.add_argument("--ds", type=int, default=2,
                    help="XY downsample factor (2 = recommended for speed/memory)")
    ap.add_argument("--do-3d", action="store_true",
                    help="true 3D segmentation (slower) instead of 2D+stitch")
    ap.add_argument("--enface", action="store_true",
                    help="RECOMMENDED for a tube: split the volume into its top "
                         "and bottom walls by depth, project each en-face, and "
                         "segment the honeycomb 2D on each (no superposition). "
                         "Total endothelial cells = top + bottom.")
    ap.add_argument("--stitch", type=float, default=0.4,
                    help="stitch_threshold linking 2D masks across Z into 3D "
                         "objects (used unless --do-3d)")
    ap.add_argument("--invert", action="store_true",
                    help="invert VE-cad so cell interiors are bright (try if the "
                         "QC overlay shows borders, not cells, being segmented)")
    ap.add_argument("--cellprob", type=float, default=0.0,
                    help="Cellpose cell-probability threshold; LOWER (e.g. -1, -2) "
                         "accepts fainter cells (helps the dim far wall)")
    ap.add_argument("--no-clahe", action="store_true",
                    help="disable per-wall CLAHE local-contrast normalization "
                         "(on by default in --enface mode)")
    ap.add_argument("--denoise", type=int, default=1,
                    help="median-filter radius for speckle denoise (0 = off)")
    ap.add_argument("--bgsub", type=float, default=20.0,
                    help="background-subtraction sigma in pixels; removes diffuse "
                         "haze so junctions stand out (0 = off). Only used with "
                         "--no-tissue.")
    ap.add_argument("--no-tissue", action="store_true",
                    help="disable tissue-aware wall split (masking + in-tissue "
                         "depth valley). On by default in --enface mode; it fixes "
                         "branch points / off-center vessels where a plain split "
                         "fed background noise to one wall.")
    ap.add_argument("--min-vox", type=int, default=50,
                    help="drop 3D objects smaller than this many voxels (debris)")
    ap.add_argument("--only", default=None,
                    help="only images whose name contains this substring")
    ap.add_argument("--zmin", type=int, default=0,
                    help="skip the first N Z-slices at read time (drop empty "
                         "space above the tissue that skews the wall split)")
    ap.add_argument("--zmax", type=int, default=0,
                    help="last Z-slice to keep (0 = to the end)")
    ap.add_argument("--no-gpu", action="store_true")
    args = ap.parse_args()

    try:
        from cellpose import models
    except ImportError:
        print("Cellpose not installed. Run:  pip install cellpose", flush=True)
        raise SystemExit(1)

    # Cellpose 3.x: models.Cellpose(model_type=...). Build once.
    try:
        model = models.Cellpose(gpu=not args.no_gpu, model_type=args.model)
        newapi = False
    except Exception:
        # Cellpose 4.x (Cellpose-SAM): single model, no model_type/channels.
        model = models.CellposeModel(gpu=not args.no_gpu)
        newapi = True

    tifs = sorted(Path(args.dir).glob(args.pattern))
    if args.only:
        tifs = [t for t in tifs if args.only in t.name]
    if not tifs:
        print(f"No images matching {args.pattern} in {args.dir}", flush=True)
        raise SystemExit(1)

    outcsv = Path(args.dir) / "endothelial_counts.csv"
    rows = []
    diam = None if args.diameter <= 0 else args.diameter
    print(f"Cellpose ({'4.x' if newapi else '3.x'}, model={args.model}, "
          f"{'3D' if args.do_3d else f'2D+stitch {args.stitch}'}, ds={args.ds}) "
          f"on {len(tifs)} image(s)\n", flush=True)

    for tif in tifs:
        print(f"  {tif.name} ...", flush=True)
        try:
            vol = _read_channel(tif, args.channel, args.ds)
            if vol is None:
                print("    no VE-Cadherin channel; skipped.", flush=True)
                continue
            if args.zmin > 0 or args.zmax > 0:
                z0, z1 = args.zmin, (args.zmax if args.zmax > 0 else vol.shape[0])
                vol = vol[z0:z1]
                print(f"    Z-cropped to [{z0}:{z1}] -> {vol.shape[0]} slices",
                      flush=True)
            saved = Path(args.dir) / f"{tif.stem}_saved"

            if args.enface:
                n_endo, counts, projs, masks2d, zs = _enface_count(
                    model, vol, newapi, diam, args.min_vox, args.invert,
                    cellprob=args.cellprob, clahe=not args.no_clahe,
                    denoise_r=args.denoise, bg_sigma=args.bgsub,
                    tissue=not args.no_tissue)
                print(f"    split at Z={zs}: top wall={counts[0]}, "
                      f"bottom wall={counts[1]} (total={n_endo})", flush=True)
                if saved.exists():
                    _enface_qc(projs, masks2d, [counts[0], counts[1]],
                               saved / f"{tif.stem}_endothelial_qc.png",
                               f"{tif.name}: {n_endo} endothelial cells (en-face)")
            else:
                img = vol.max() - vol if args.invert else vol
                # (Z, Y, X) grayscale: tell Cellpose the Z axis, no channel axis.
                eval_kw = dict(diameter=diam, z_axis=0, channel_axis=None)
                if args.do_3d:
                    eval_kw.update(do_3D=True)
                else:
                    eval_kw.update(do_3D=False, stitch_threshold=args.stitch)
                if not newapi:
                    eval_kw.update(channels=[0, 0])
                masks = model.eval(img, **eval_kw)[0]
                if args.min_vox > 0 and masks.max() > 0:
                    ids, cnt = np.unique(masks[masks > 0], return_counts=True)
                    small = ids[cnt < args.min_vox]
                    if small.size:
                        masks = np.where(np.isin(masks, small), 0, masks)
                n_endo = int(len(np.unique(masks[masks > 0])))
                if saved.exists():
                    _qc_overlay(vol, masks, saved / f"{tif.stem}_endothelial_qc.png",
                                f"{tif.name}: {n_endo} endothelial cells")

            n_tot = _total_nuclei(saved)
            n_fib = (n_tot - n_endo) if n_tot is not None else None
            pct = (100.0 * n_fib / n_tot) if (n_fib is not None and n_tot) else None
            rows.append({"image": tif.name, "n_endothelial_cellpose": n_endo,
                         "n_total_nuclei": n_tot if n_tot is not None else "",
                         "n_fibroblasts": n_fib if n_fib is not None else "",
                         "pct_fibroblast": f"{pct:.1f}" if pct is not None else ""})
            print(f"    endothelial cells = {n_endo}"
                  + (f"; total nuclei = {n_tot}; fibroblasts = {n_fib} "
                     f"({pct:.0f}%)" if n_fib is not None else
                     "; (no fibroblast_count.npz -> total nuclei unknown)"),
                  flush=True)
        except Exception as exc:
            import traceback
            print(f"    FAILED: {exc}\n{traceback.format_exc()}", flush=True)

    if rows:
        with open(outcsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\nWrote {outcsv} ({len(rows)} image(s)).", flush=True)


if __name__ == "__main__":
    main()
