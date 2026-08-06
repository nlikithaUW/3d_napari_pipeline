"""Batch max-intensity-projection PNGs of the VE-Cadherin channel.

Standalone — no pipeline imports. For every OME-TIFF in a folder it finds the
VE-Cadherin channel (by OME channel name), takes the maximum-intensity
projection over Z, contrast-stretches it, and writes a PNG into a subfolder.

Usage (from anywhere, in an env with tifffile + numpy + PIL/imageio):

    python Scripts/vecad_mip.py --dir "C:\\...\\images\\WT"

Options:
    --dir         folder of OME-TIFFs (default: the WT images folder)
    --pattern     glob for images inside --dir (default: *.tif)
    --subfolder   output subfolder name, created inside --dir (default: vecad_mip)
    --channel     OME channel name to project (default: VE-Cadherin; falls back
                  to a case-insensitive match containing "cad")
    --pmin --pmax percentile contrast stretch (default 1.0 / 99.9); set
                  --pmin 0 --pmax 100 for a plain min/max stretch
    --bits        8 (default, stretched) or 16 (raw MIP, no stretch)
    --overwrite   re-make PNGs that already exist (default: skip them)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import tifffile

DEFAULT_DIR = (Path(__file__).resolve().parent.parent / "3d_images_folder" /
               "images" / "WT")


def _save_png(path: Path, img: np.ndarray) -> None:
    """Write a 2D uint8/uint16 array as PNG (PIL, else imageio)."""
    try:
        from PIL import Image
        Image.fromarray(img).save(str(path))
        return
    except Exception:
        pass
    import imageio.v2 as imageio
    imageio.imwrite(str(path), img)


def _channel_index(xml: str, want: str) -> int | None:
    names = re.findall(r'<Channel[^/>]*\sName="([^"]+)"', xml)
    if not names:
        return None
    for i, n in enumerate(names):
        if n == want:
            return i
    for i, n in enumerate(names):          # fallback: fuzzy "cad" match
        if "cad" in n.lower():
            return i
    return None


def _mip(tif: Path, channel: str) -> tuple[np.ndarray, str] | None:
    """Max-project the wanted channel over Z, reading it PLANE-BY-PLANE.

    Reading the whole stack at once (series.asarray) fails on some files with
    a libdeflate bad-data error — a corrupt compressed chunk in *another*
    channel. We only need this one channel, so we decode just its Z-planes and
    skip any that raise; the MIP is robust to a few missing slices.
    """
    with tifffile.TiffFile(str(tif)) as tf:
        series = tf.series[0]
        axes = series.axes
        xml = tf.ome_metadata or ""
        names = re.findall(r'<Channel[^/>]*\sName="([^"]+)"', xml)
        cidx = _channel_index(xml, channel)
        if cidx is None:
            print(f"  {tif.name}: no VE-Cadherin channel in {names}; skipped.",
                  flush=True)
            return None
        cname = names[cidx]

        # Map (z, channel) -> flat page index. OME stores planes in the series'
        # axis order; for the usual ZCYX that's z*C + c, for CZYX it's c*Z + z.
        if axes == "ZCYX":
            Z, C = series.shape[0], series.shape[1]
            page_of = lambda z: z * C + cidx
        elif axes == "CZYX":
            C, Z = series.shape[0], series.shape[1]
            page_of = lambda z: cidx * Z + z
        else:                                   # uncommon layout: whole-stack read
            data = series.asarray()
            cpos = axes.find("C")
            vol = np.take(data, cidx, axis=cpos)
            zpos = [a for a in axes if a != "C"].index("Z")
            return np.moveaxis(vol, zpos, 0).max(axis=0), cname

        mip = None
        good = bad = 0
        for z in range(Z):
            try:
                plane = tf.pages[page_of(z)].asarray()
            except Exception:
                bad += 1
                continue
            mip = plane if mip is None else np.maximum(mip, plane)
            good += 1

    if mip is None:
        print(f"  {tif.name}: every {cname} plane failed to decode; skipped.",
              flush=True)
        return None
    if bad:
        print(f"  {tif.name}: {bad}/{Z} {cname} planes were corrupt and "
              f"skipped ({good} used).", flush=True)
    return mip, cname                            # MIP over Z


def _stretch8(mip: np.ndarray, pmin: float, pmax: float) -> np.ndarray:
    lo, hi = np.percentile(mip, [pmin, pmax])
    if hi <= lo:
        lo, hi = float(mip.min()), float(max(mip.max(), mip.min() + 1))
    out = (np.clip((mip.astype(np.float32) - lo) / (hi - lo), 0, 1) * 255)
    return out.astype(np.uint8)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--pattern", default="*.tif")
    ap.add_argument("--subfolder", default="vecad_mip")
    ap.add_argument("--channel", default="VE-Cadherin")
    ap.add_argument("--pmin", type=float, default=1.0)
    ap.add_argument("--pmax", type=float, default=99.9)
    ap.add_argument("--bits", type=int, choices=(8, 16), default=8)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    tifs = sorted(Path(args.dir).glob(args.pattern))
    if not tifs:
        print(f"No images matching {args.pattern} in {args.dir}", flush=True)
        sys.exit(1)
    outdir = Path(args.dir) / args.subfolder
    outdir.mkdir(exist_ok=True)
    print(f"VE-Cadherin MIP: {len(tifs)} image(s) -> {outdir}", flush=True)

    ok = 0
    for tif in tifs:
        out = outdir / f"{tif.stem}_vecad_mip.png"
        if out.exists() and not args.overwrite:
            print(f"  skip (exists): {out.name}", flush=True)
            ok += 1
            continue
        try:
            res = _mip(tif, args.channel)
            if res is None:
                continue
            mip, cname = res
            img = mip.astype(np.uint16) if args.bits == 16 else _stretch8(
                mip, args.pmin, args.pmax)
            _save_png(out, img)
            print(f"  {tif.name}  [{cname}]  {mip.shape}  {args.bits}-bit "
                  f"-> {out.name}", flush=True)
            ok += 1
        except Exception as exc:
            print(f"  FAILED {tif.name}: {exc}", flush=True)

    print(f"\nDone: {ok}/{len(tifs)} PNG(s) in {outdir}", flush=True)


if __name__ == "__main__":
    main()
