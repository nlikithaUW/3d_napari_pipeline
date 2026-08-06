"""Is the tension readout confounded by imaging depth?

Depth and curvature covary on a tube — adhesions on the far wall are both
deeper and differently curved than those on the near wall. So if the
mClover/mRuby ratio drifts with depth for optical reasons, that drift shows up
as a curvature effect. This script tests for it before any such claim is made.

    python Scripts/depth_check.py --dir 3d_images_folder/images/WT
    python Scripts/depth_check.py --dir ... --channels     # + the mechanism

Reads the newest `*_adhesion_measures_*.csv` in each `*_saved/` folder. Reports
per image first and pools only afterwards: adhesions are not independent
observations (many per cell, few images), so a pooled p-value is driven by
adhesion count, not by evidence. Read the effect sizes, not the stars.

`--channels` additionally loads the raw stack and the saved adhesion instances
to split the ratio into its numerator and denominator. That answers *why* any
drift exists: if clover and ruby attenuate together the ratio is safe, and only
a differential decay actually biases it. Slow (loads the full stack per image).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BIN_UM = 20.0
Z = "centroid_z_um"
R = "ratio_clover_ruby"


def newest_csvs(root: Path) -> list[Path]:
    """Newest measurement CSV per `*_saved` folder."""
    out = []
    for saved in sorted(root.glob("*_saved")):
        csvs = sorted(saved.glob("*_adhesion_measures_*.csv"))
        if csvs:
            out.append(max(csvs, key=lambda p: p.stat().st_mtime))
    return out


def _rho(x, y) -> tuple[float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.ptp(x[m]) == 0 or np.ptp(y[m]) == 0:
        return float("nan"), float("nan"), int(m.sum())
    rho, p = spearmanr(x[m], y[m])
    return float(rho), float(p), int(m.sum())


def quartile_contrast(z, v) -> tuple[float, float, float]:
    """Median of v in the shallowest vs deepest depth quartile, and % change.

    A contrast between the extremes is more robust than a slope: the depth
    distribution on a tube is bimodal (near wall / far wall), not uniform, so a
    fitted line is dominated by whichever wall has more adhesions.
    """
    m = np.isfinite(z) & np.isfinite(v)
    z, v = z[m], v[m]
    lo, hi = v[z <= np.percentile(z, 25)], v[z >= np.percentile(z, 75)]
    a, b = float(np.median(lo)), float(np.median(hi))
    return a, b, (b / a - 1.0) * 100.0 if a else float("nan")


def depth_bins(z, v, bin_um=BIN_UM):
    edges = np.arange(0, np.ceil(np.nanmax(z) / bin_um) * bin_um + 1, bin_um)
    mids, meds, ns = [], [], []
    for i in range(len(edges) - 1):
        m = (z >= edges[i]) & (z < edges[i + 1]) & np.isfinite(v)
        if m.sum() >= 5:
            mids.append((edges[i] + edges[i + 1]) / 2)
            meds.append(float(np.median(v[m])))
            ns.append(int(m.sum()))
    return np.array(mids), np.array(meds), np.array(ns)


def per_adhesion_channels(csv: Path):
    """Per-adhesion clover/ruby totals per voxel, background-subtracted.

    Recomputed from the raw stack + saved instances, matching what Measure does,
    so the numerator and denominator can be looked at separately.
    """
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from pipeline.io import load_params, load_stack
    from pipeline.measure_core import _bs_per_z

    saved = csv.parent
    stem = saved.name[: -len("_saved")]
    tif = saved.parent / f"{stem}.tif"
    inst = saved / "adhesion_instances.npz"
    if not (tif.exists() and inst.exists()):
        print(f"    (--channels: need {tif.name} and adhesion_instances.npz)")
        return None

    params = load_params(Path(__file__).resolve().parent.parent
                         / "3d_images_folder" / "params.json")
    st = load_stack(tif, params)
    if not {"clover", "ruby"} <= set(st.roles):
        print("    (--channels: stack has no clover/ruby roles)")
        return None
    labels = np.load(inst)["instances"].astype(np.int32)
    clover = st.data[:, st.channels.index(st.roles["clover"])]
    ruby = st.data[:, st.channels.index(st.roles["ruby"])]

    # Map isotropic labels down to native Z, exactly as Measure does.
    iso_z = 0.325 if not len(st.voxel_zyx_um) else min(st.voxel_zyx_um)
    idx = np.clip(np.rint(np.arange(clover.shape[0])
                          * (st.voxel_zyx_um[0] / iso_z)).astype(np.int64),
                  0, labels.shape[0] - 1)
    ln = labels[idx]
    ids = np.unique(labels[labels > 0])
    nmax = int(labels.max())
    cb, rb = _bs_per_z(clover), _bs_per_z(ruby)

    # The background level each plane had removed. For any voxel above the
    # mode, raw - subtracted == the mode exactly, so the per-plane max of the
    # difference recovers it without re-implementing the (frozen) mode logic.
    bg_c = float(np.median([(clover[z] - cb[z]).max() for z in range(len(cb))]))
    bg_r = float(np.median([(ruby[z] - rb[z]).max() for z in range(len(rb))]))

    ins = ln > 0
    lv = ln[ins]
    n = np.bincount(lv, minlength=nmax + 1)[ids].astype(np.float64)
    c = np.bincount(lv, weights=cb[ins], minlength=nmax + 1)[ids]
    r = np.bincount(lv, weights=rb[ins], minlength=nmax + 1)[ids]
    with np.errstate(invalid="ignore", divide="ignore"):
        df = pd.DataFrame({"label": ids, "clover_per_vox": c / n,
                           "ruby_per_vox": r / n})
    df.attrs["bg"] = (bg_c, bg_r)
    return df


def report(df: pd.DataFrame, name: str, channels: bool) -> dict:
    z, ratio = df[Z].to_numpy(), df[R].to_numpy()
    print(f"\n=== {name}  (n = {len(df)} adhesions) ===")
    print(f"  depth range {np.nanmin(z):.1f} – {np.nanmax(z):.1f} µm")

    rho, p, n = _rho(z, ratio)
    print(f"\n  ratio vs depth:   rho = {rho:+.3f}   (p = {p:.1e}, n = {n})")
    lo, hi, pct = quartile_contrast(z, ratio)
    print(f"  shallowest 25%: {lo:.2f}   deepest 25%: {hi:.2f}   "
          f"change {pct:+.0f}%")

    mids, meds, ns = depth_bins(z, ratio)
    print(f"\n  {'depth (µm)':>12s} {'n':>6s} {'median ratio':>13s}")
    for m_, md, n_ in zip(mids, meds, ns):
        print(f"  {m_ - BIN_UM/2:5.0f}–{m_ + BIN_UM/2:<5.0f} {n_:6d} {md:13.2f}")

    # Does curvature covary with depth? That is what makes depth a confound
    # rather than just a nuisance.
    print("\n  covariation with depth (why this matters):")
    for col in ("curv_mean", "curv_gaussian", "path_curv_xy", "dist_mesh_um"):
        if col in df:
            r_, p_, n_ = _rho(z, df[col].to_numpy())
            print(f"    {col:16s} rho = {r_:+.3f}  (p = {p_:.1e})")

    extra = {}
    if channels:
        ch = per_adhesion_channels(Path(df.attrs["csv"]))
        if ch is not None:
            m = df.merge(ch, on="label", how="inner")
            print("\n  mechanism — numerator vs denominator vs depth:")
            for col in ("clover_per_vox", "ruby_per_vox"):
                r_, p_, n_ = _rho(m[Z].to_numpy(), m[col].to_numpy())
                a, b, pc = quartile_contrast(m[Z].to_numpy(), m[col].to_numpy())
                print(f"    {col:16s} rho = {r_:+.3f}  "
                      f"shallow {a:8.1f} -> deep {b:8.1f}  ({pc:+.0f}%)")
            print("    (if both decay ~equally the ratio survives; only a "
                  "differential\n     decay biases it)")

            # How much signal is actually there? A denominator only a few
            # percent above background makes the ratio fragile no matter how
            # the depth question resolves.
            bg_c, bg_r = ch.attrs["bg"]
            print("\n  signal quality — contrast over background:")
            for col, bg in (("clover_per_vox", bg_c), ("ruby_per_vox", bg_r)):
                sig = float(np.nanmedian(m[col]))
                print(f"    {col:16s} {sig:8.1f} counts/voxel over bg {bg:5.0f}"
                      f"  ({100 * sig / bg:5.1f}% contrast)")
            sig_r = float(np.nanmedian(m["ruby_per_vox"]))
            if sig_r > 0:
                print(f"    fret.background_subtract quantises the background to "
                      f"bin_size=10, so\n    that is ~±5 counts on a "
                      f"{sig_r:.0f}-count signal → ~{100 * 5 / sig_r:.0f}% "
                      f"systematic error on every ratio.")
            rc, _, _ = _rho(m[R].to_numpy(), m["clover_per_vox"].to_numpy())
            rr, _, _ = _rho(m[R].to_numpy(), m["ruby_per_vox"].to_numpy())
            print(f"\n    ratio tracks clover rho={rc:+.3f} | ruby rho={rr:+.3f}"
                  f"  (ratio = c/r couples to r by\n     construction, so read "
                  f"this alongside the contrast figures above)")
            extra["merged"] = m
    return {"z": z, "ratio": ratio, "bins": (mids, meds, ns), "name": name, **extra}


def _binned_line(ax, mids, meds, ns, color, label):
    """Binned medians with marker area ∝ n.

    Depth is bimodal on a tube (near wall / far wall), so the middle bins can
    hold a few dozen adhesions while the end bins hold over a thousand. Drawn
    at equal weight those sparse bins read as structure; sized by n they read
    as what they are.
    """
    if not len(mids):
        return
    ax.plot(mids, meds, "-", color=color, lw=1.5, alpha=0.8, zorder=2)
    ax.scatter(mids, meds, s=18 + 120 * np.sqrt(ns / ns.max()), color=color,
               edgecolor="white", linewidth=0.8, zorder=3, label=label)
    for x, y, n_ in zip(mids, meds, ns):
        ax.annotate(f"n={n_}", (x, y), textcoords="offset points",
                    xytext=(0, 11), ha="center", fontsize=6.5, color="0.35")


def make_plot(results: list[dict], out: Path) -> None:
    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(6.5 * n, 9.5), squeeze=False)
    for j, res in enumerate(results):
        z, ratio = res["z"], res["ratio"]
        ax = axes[0][j]
        m = np.isfinite(z) & np.isfinite(ratio)
        ax.hexbin(z[m], ratio[m], gridsize=40, cmap="Blues", mincnt=1)
        mids, meds, ns = res["bins"]
        _binned_line(ax, mids, meds, ns, "crimson",
                     f"median / {BIN_UM:.0f} µm (marker ∝ n)")
        lo, hi, pct = quartile_contrast(z, ratio)
        ax.axhline(lo, color="0.45", lw=1, ls=":")
        ax.axhline(hi, color="0.45", lw=1, ls=":")
        ax.set_ylim(0, np.nanpercentile(ratio, 99.5) * 1.1)
        ax.set_xlabel("adhesion depth  centroid_z_um (µm)")
        ax.set_ylabel("ratio_clover_ruby")
        ax.set_title(f"{res['name']}\ntension readout vs depth  "
                     f"(shallow 25% {lo:.1f} → deep 25% {hi:.1f}, {pct:+.0f}%)",
                     fontsize=10)
        ax.legend(loc="upper right", fontsize=8)

        ax2 = axes[1][j]
        if "merged" in res:
            mm = res["merged"]
            zz = mm[Z].to_numpy()
            for col, c in (("clover_per_vox", "seagreen"),
                           ("ruby_per_vox", "indianred")):
                mi, me, nn = depth_bins(zz, mm[col].to_numpy())
                if len(mi):
                    # Normalised to the shallowest bin so two channels of very
                    # different brightness share one axis.
                    _binned_line(ax2, mi, me / me[0], nn, c, col)
            ax2.axhline(1.0, color="0.6", lw=1, ls="--")
            ax2.set_ylabel("signal / voxel, relative to shallowest bin")
            ax2.set_title("mechanism: do both channels attenuate together?",
                          fontsize=10)
        else:
            mi, me, nn = depth_bins(z, ratio)
            _binned_line(ax2, mi, me, nn, "crimson", "median ratio")
            ax2.set_ylabel("median ratio")
            ax2.set_title("binned median (run --channels for the mechanism)",
                          fontsize=10)
        ax2.legend(fontsize=8)
        ax2.set_xlabel("adhesion depth  centroid_z_um (µm)")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"\nPlot -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, type=Path,
                    help="folder of images, e.g. 3d_images_folder/images/WT")
    ap.add_argument("--channels", action="store_true",
                    help="also split the ratio into clover / ruby (slow: "
                         "reloads each raw stack)")
    ap.add_argument("--out", type=Path, default=None,
                    help="plot path (default: <dir>/depth_check.png)")
    args = ap.parse_args()

    csvs = newest_csvs(args.dir)
    if not csvs:
        print(f"No *_adhesion_measures_*.csv under {args.dir}/*_saved/")
        return
    print(f"Found {len(csvs)} image(s) with measurements.")

    results, frames = [], []
    for c in csvs:
        df = pd.read_csv(c)
        missing = {Z, R} - set(df.columns)
        if missing:
            print(f"\n{c.name}: missing {sorted(missing)} — this CSV predates "
                  f"the centroid columns; re-export it.")
            continue
        df.attrs["csv"] = str(c)
        results.append(report(df, c.parent.name[: -len("_saved")], args.channels))
        frames.append(df)

    if not results:
        return

    if len(frames) > 1:
        print("\n\n=== pooled across images ===")
        allf = pd.concat(frames, ignore_index=True)
        rho, p, n = _rho(allf[Z].to_numpy(), allf[R].to_numpy())
        print(f"  ratio vs depth: rho = {rho:+.3f} (p = {p:.1e}, n = {n})")
        signs = [np.sign(_rho(r["z"], r["ratio"])[0]) for r in results]
        agree = int(max(signs.count(1), signs.count(-1)))
        print(f"  direction consistent in {agree}/{len(signs)} images")
        print("  NB: the pooled p-value assumes independent adhesions, which "
              "these\n      are not. Consistency of direction across images is "
              "the real evidence.")
    else:
        print("\n\nOnly one image — no cross-image check possible. A single "
              "image cannot\ndistinguish a depth effect from a property of "
              "this particular field.")

    make_plot(results, args.out or args.dir / "depth_check.png")


if __name__ == "__main__":
    main()
