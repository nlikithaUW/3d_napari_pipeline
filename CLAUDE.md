# 3d_napari_pipeline — read this first

Refer to this file at the start of every session in this repo.

## What this project is

A clean-slate rewrite of an older `napari_pipeline` (sibling repo at
`/mnt/d/OneDrive - UW/Documents/GitHub/napari_pipeline`). The old one got
tangled handling multiple data types; this one starts from scratch and
stays narrow.

## Biology / imaging

- Cardiac cells cultured on a 3D **vessel / lumen** surface.
- Imaged as z-stacks. Voxel size **0.325 × 0.325 × 0.4 µm** (xy × z).
  Every 3D kernel must scale its z component for this anisotropy.
- Condition: `WT` (folder under `3d_images_folder/images/`).
- Channels (see `3d_images_folder/params.json`):
  - `FRET` — FRET channel
  - `mClover` — donor
  - `mRuby` — acceptor
  - `VE-Cadherin` — substrate / cell-cell junctions, **also the only
    proxy for the vessel surface** (no dedicated vessel stain)
  - `Nuclei` — DAPI
- A cell = DAPI nucleus + a VE-cadherin **ring**. The ring is embedded
  in the curved vessel surface, so the **union of rings approximates
  the vessel surface** — that's how we'll infer geometry later.

## End goal (long term)

Segment cells (DAPI + cadherin) and adhesions with **uSegment3D**,
reconstruct the vessel surface from cadherin, then relate per-cell /
per-adhesion **force** (from FRET LUT) to cell and vessel geometry.

## What exists now

- `Scripts/nd2_to_ometiff.py` — ND2 → OME-TIFF, drops junk channels, writes
  pixel size + channel names into OME metadata.
- `pipeline/io.py` — role-addressable loader (`Stack.by_role("clover")`).
- `pipeline/fret.py` — per-Z background subtraction + bleedthrough
  correction + force LUT lookup. **Math is verbatim from the old
  pipeline; do not modify without flagging.**
- `pipeline/view.py` — napari 3D viewer with raw channels + FRET
  efficiency + force overlays.
- `Fret_LUT/Clover_mRuby2_GGSGGS7_force_new.txt` — efficiency→force LUT.

## Environment

- **Dev (mine):** WSL2, conda env `image_analysis` (see `environment.yml`).
- **Run (user's):** Windows Anaconda terminal, same env name (separate env,
  same name `image_analysis`).
- **GUI stack is pinned** to `napari=0.6.0` + `vispy=0.14.3`. vispy 0.15
  breaks 3D rendering (`Cannot SIZE object N` GLIR crash). Keep dev and run
  envs on the same pinned versions.
- **Hardware:** 32/64 c/t CPU, RTX 3090. GPU is available but unused so far.
- **If I add a dependency, I must tell the user to install it in their
  Windows conda env.** I cannot install it for them.

## Working agreement (durable preferences)

- **Do not code until the user explicitly says so.** When given a new
  task, propose the plan, ask clarifying questions, wait for approval.
- Tunable parameters → **magicgui sliders in napari** for live tuning,
  not hard-coded defaults. The user wants to see and adjust.
- Keep the pipeline narrow. Resist adding abstractions for hypothetical
  future data types — that's exactly what sank the previous version.
- FRET math is load-bearing and copied from the old pipeline; treat it
  as frozen unless the user asks otherwise.

## Current focus

Noise filtering + clean 3D rendering, **display-only** (no persisted
filtered volumes yet). Target channels: DAPI and VE-cadherin. Leave
FRET / mClover / mRuby / derived eff+force raw for now. Segmentation
and vessel-surface inference come after the user is happy with what
the viewer shows.
