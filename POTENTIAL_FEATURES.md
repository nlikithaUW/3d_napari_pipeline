# Potential features — backlog

Ideas considered but deferred. Not yet implemented.

---

## Progress / completion bars for long-running operations

**Want:** when a tab's compute button is clicked (e.g. Nuclei →
"Compute DAPI Otsu + components", Adhesion RF → "Run RF predict"), show how far
along it is and when it's done, instead of only terminal prints.

**Key nuance — not all steps can show a true percentage:**

- **Determinate (real %) is possible** for steps that loop over items:
  - `rf.predict_proba` already runs in N z-chunks and prints `chunk i/N`
    (rf.py ~line 436) — a natural progress source.
  - Load-time resampling loops over channels (`io.to_isotropic`) and FRET
    outputs (`fret.resample_fret_outputs`).
  - FRET compute loops per-Z.
  These would need the blocking worker converted to a generator that `yield`s
  progress (and `rf.predict_proba` would need a progress callback threaded out
  of it).
- **Indeterminate only** for atomic single library calls, which have no
  sub-progress hook: `otsu_mask` (`threshold_otsu`), `label_and_sizes`
  (`ndi.label`), `niblack_mask`, watershed splitting. Best that's honest here is
  a "working…" bar that animates and clears on completion.

**Implementation options:**

1. **Inline `magicgui.ProgressBar` per tab** — bar appears in the tab, right
   where the button is. Matches the request most closely. More work: wire each
   worker's `yielded`/`started`/`finished` signals to a per-tab bar, and make
   the determinate workers yield progress.
2. **napari's built-in activity dock** (`napari.utils.progress`) — the
   hourglass panel, bottom-right. Much less code and idiomatic, but off to the
   side rather than in the tab.

**Open decisions (unresolved):** inline vs activity dock; scope (all tabs vs
just the heavy steps vs a Nuclei-tab prototype first); whether atomic steps get
a busy/indeterminate bar or no bar at all.
