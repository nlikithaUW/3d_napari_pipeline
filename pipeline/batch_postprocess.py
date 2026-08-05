"""Headless re-run of the adhesion post-process at the standard settings.

Targets images that were measured *before* settings-logging existed — i.e. their
`<stem>_saved` folder has a measurements CSV but **no** `_settings.json`. It
re-labels `adhesion_instances.npz` from the saved probability map at the standard
post-process (threshold 0.5, min_size 100, split ON, split_min 1000), reusing the
pipeline's own `instance_labels_3d`.

Pure Python — no napari, no GL, no GPU — so it is safe to run in a second
terminal *alongside* the main `pipeline.batch` run. Because it requires a CSV to
be present, it skips folders the main batch is still writing (those have no CSV
yet), so the two never collide.

It does NOT run Measure/export. After it finishes, re-run Measure + export per
image in the GUI (Load all saved -> Measure -> export) to regenerate the CSV and
the settings JSON with the standardized instances.

Usage (image_analysis env, from the 3d_napari_pipeline folder):

    python -m pipeline.batch_postprocess --dir "C:\\...\\images\\WT"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .rf import instance_labels_3d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True, type=Path,
                    help="folder containing the <stem>_saved subfolders")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-size", type=int, default=100)
    ap.add_argument("--split-min", type=int, default=1000)
    ap.add_argument("--all-no-json", action="store_true",
                    help="process EVERY folder with a proba and no settings JSON "
                         "(not just those with a CSV). WARNING: this will also "
                         "hit folders your main batch is actively writing.")
    args = ap.parse_args()

    targets = []
    for saved in sorted(Path(args.dir).glob("*_saved")):
        has_proba = (saved / "adhesion_proba.npy").exists()
        has_csv = bool(list(saved.glob("*adhesion_measures*.csv")))
        has_json = bool(list(saved.glob("*_settings.json")))
        if not has_proba or has_json:
            continue
        if args.all_no_json or has_csv:
            targets.append(saved)

    if not targets:
        print("Nothing to do: no matching folders (proba present, no settings "
              "JSON, and — unless --all-no-json — a CSV present).", flush=True)
        return

    T = int(round(args.threshold * 255))
    print(f"Re-post-processing {len(targets)} folder(s) at standard settings "
          f"(threshold={args.threshold}, min_size={args.min_size}, split ON, "
          f"split_min={args.split_min}):", flush=True)
    for saved in targets:
        try:
            print(f"  {saved.name} ...", flush=True)
            proba = np.load(saved / "adhesion_proba.npy")
            binary = (proba >= T).astype(np.uint8)
            labels = instance_labels_3d(binary, args.min_size, True,
                                        args.split_min)
            n = int(labels.max())
            np.savez_compressed(saved / "adhesion_instances.npz",
                                instances=labels.astype(np.int32))
            print(f"    -> {n} adhesions; adhesion_instances.npz updated.",
                  flush=True)
            del proba, binary, labels
        except Exception as exc:
            print(f"    FAILED ({exc}); skipping.", flush=True)

    print("\nDone. Re-run Measure + export per image in the GUI (Load all saved "
          "-> Measure -> export) to regenerate each CSV + settings JSON with the "
          "standardized instances.", flush=True)


if __name__ == "__main__":
    main()
