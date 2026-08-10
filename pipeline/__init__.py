"""3D cardiac-vessel image pipeline: load, segment, mesh, and measure.

Layout:

- `io`          — OME-TIFF loading, role-addressable channels, isotropic resampling.
- `fret`        — FRET math (background subtract, bleedthrough, force LUT). FROZEN:
                  copied verbatim from the old pipeline; do not edit.
- `gpu`         — optional cupy shim for windowed ndimage ops, CPU fallback.
- `denoise`     — thresholding (Otsu/Niblack/Sauvola/Bernsen) + size filtering.
- `rf`          — random-forest voxel classifier + instance labeling.
- `cell_mesh`   — per-cell surfaces, vertex normals, quadric surface curvature.
- `vessel_mesh` — vessel surface from the VE-cadherin mask (alpha shape).
- `centerline`  — in-plane vessel centerline + path curvature.
- `measure_core`— per-adhesion measurements, headless.
- `batch*`      — headless multi-image runners.
- `view`        — the napari app; `gui/` holds one module per tab.

Start with `python -m pipeline.view`. See README.md for configuration and
docs/ISSUES.md for what is known to be wrong.
"""
