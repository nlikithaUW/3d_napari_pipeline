"""Adhesion orientation tab — which way the adhesions point, and whether the
neighbours agree.

Everything orientation-related lives here. The maths is in
`pipeline/adhesion_orient.py` (pure functions, no Qt) so that a headless caller
runs the *same* code — CLAUDE.md's standing rule, and the one this project has
broken three times.

The tab is four collapsible groups, the same layout doctrine as the Vessel Mesh
tab: inputs above the button, display-only controls below it, one status line
per group saying what is cached and whether it was computed on the mesh
currently on screen.

  1. **Orientation** — the compute. Reads the adhesion label volume, the vessel
     mesh and the centerline graph; writes nothing to disk.
  2. **Per-adhesion axes** — one arrow per adhesion. This is the raw data and it
     is the view that *misleads*: read group 3 before drawing conclusions from
     it.
  3. **Binned average glyphs** — one arrow per patch of vessel. The reason this
     tab exists.
  4. **Local coherence** — per-adhesion agreement with its neighbours, with no
     binning choice at all, as the control on group 3.

**Nothing here is saved.** No `.npz`, no CSV column, no batch step. That is
deliberate: the orientation measurements are still being argued about
(`docs/scratch/orientation/`), and a saved artifact or an exported column is a
decision about defaults that `batch.py` would then inherit silently. When the
measurement settles, promoting it means adding an artifact to `saved_io.py` and
a step to `state.RUN_STEPS` — not copying any maths out of the core module.

Why the binned view is the headline
-----------------------------------
Measured on `04step_003` (`docs/scratch/orientation/`):

  * adhesions do **not** run along the vessel. Compared in the xy projection —
    the view you are looking at when you look down z — their axes are as often
    across the tube as along it, inside the permutation null at every threshold.
  * but neighbouring adhesions **do** agree with each other, strongly out to a
    few µm and decaying to nothing by ~100 µm, and that survives measuring it in
    the xy plane alone where the z-elongated PSF cannot produce it.

Those two together predict that local patches have coherent directions which
*wander* independently of the vessel. That is what group 3 draws, and
``axial_cos`` on each glyph is the direct test: it is the cosine between the
patch's mean axis and the patch's mean vessel axis.

The z axis is not neutral
-------------------------
The PSF is elongated in z, so nearly every measured axis leans that way. The
``plane`` control exists for this: in **xy** the blur imposes no preferred
compass direction, so an alignment seen there cannot be the microscope. The
colour-by menu marks the lab-frame-contaminated choices with a warning glyph
rather than hiding them — `into_deg` in particular is mostly a statement about
where on the vessel you are standing, not about the adhesion.
"""

from __future__ import annotations

import numpy as np
from magicgui.widgets import ComboBox, FloatSlider, Label, PushButton, Slider
from napari.qt.threading import thread_worker

from .. import adhesion_orient as ao
from .widgets import call_all, collapsible_group, vector_mask, with_after

AXES_LAYER = "Adhesion axes"
GLYPH_LAYER = "Adhesion orientation glyphs"
COH_LAYER = "Adhesion local coherence"

#: Colour-by choices for the per-adhesion arrows. The warning glyph is not
#: decoration: those three are ~0.85 rank-correlated with the axis's own angle
#: to the stack z axis, so they mostly report the lab frame.
AXIS_COLOR_BY = {
    "local coherence": "local_coh",
    "wrap_deg (0 = along tube, 90 = around)": "wrap_deg",
    "elongation": "elong",
    "⚑ into_deg (lab-frame contaminated)": "into_deg",
    "⚑ angle to stack z (the lab frame itself)": "e1_z_deg",
}
GLYPH_COLOR_BY = {
    "coherence (how much the patch agrees)": "coherence",
    "axial_cos (1 = patch runs along the vessel)": "axial_cos",
    "count (adhesions in the patch)": "count",
}
PLANES = {"xy only (immune to the z-elongated PSF)": "xy",
          "full 3D (⚑ z-dominated)": "3d"}

#: How the per-adhesion arrow is drawn. The default is NOT the full 3D axis,
#: and that is deliberate. |e1_z| has a median of 0.967 on this data, so a 3D
#: arrow points almost straight at the camera in the standard top-down view: a
#: 4 µm arrow shows a median of 1.0 µm in plane, and that visible stub is the
#: small noisy residual rather than the direction being measured. Drawn that
#: way the arrows read as unrelated to the adhesions they sit on — which is
#: how this was found. The xy component is both what is visible and what is
#: measurable, so it is what is drawn.
AXIS_DRAW = {"xy component only (what is measurable)": "xy",
             "⚑ full 3D axis (mostly points along z)": "3d"}

#: Mean |cos| between two independent uniformly-random axes — 1/2 in 3D, 2/pi
#: in a plane. Printed only as context, and explicitly NOT as the baseline: the
#: vessel runs one general direction across a field of view, so this constant is
#: beaten by any global lean in the adhesion axes with no patch-by-patch
#: relationship at all. `ao.shuffled_axial_cos` is the baseline that holds.
ISOTROPIC_AXIAL_COS = {"3d": 0.5, "xy": 2.0 / np.pi}


def _adhesion_labels(state, viewer):
    """Full-res adhesion label volume from the Adhesion RF tab's layer, or None.

    `layer.data` is a `MultiScaleData` when the layer is multiscale, and
    `np.asarray` on that silently returns the *coarsest* level — the gotcha the
    other tabs guard the same way.
    """
    if "Adhesion instances" not in viewer.layers:
        return None
    lyr = viewer.layers["Adhesion instances"]
    data = lyr.data
    return data[0] if getattr(lyr, "multiscale", False) else data


def build(state, viewer) -> list:
    orient_status_w = Label(value="", label="")
    glyph_status_w = Label(value="", label="")

    def _result():
        return state.get("adhesion_orient")

    def _refresh_status(*_):
        r = _result()
        if not r:
            orient_status_w.value = "1. Orientation — not computed"
            return
        bid = r.get("build_id")
        note = (" — computed on the current mesh"
                if bid is not None and int(bid) == int(
                    state.get("vessel_build_id", 0))
                else " — computed on a PREVIOUS mesh, recompute")
        orient_status_w.value = (
            f"1. Orientation — {r['n_gated']:,} of {r['n_total']:,} adhesions "
            f"pass the gate ({r['n_degenerate']:,} frame-degenerate)" + note)

    # ---------------------------------------------------------------- group 1
    elong_w = FloatSlider(
        value=ao.ELONG_MIN_DEFAULT, min=1.0, max=4.0, step=0.05,
        label="elong_min",
        tooltip="Minimum elongation for an adhesion's axis to be used. Below "
                "about 1.3 a blob is round and its major axis is noise, not a "
                "direction. Raising this keeps fewer, better-defined axes.")
    maxdist_w = FloatSlider(
        value=8.0, min=1.0, max=40.0, step=0.5, label="max_dist_um",
        tooltip="Keep only adhesions within this distance of the vessel "
                "surface. Further out there is no meaningful 'along' or "
                "'around' the tube — the fibroblasts are cast through the gel, "
                "so many adhesions are genuinely far from the wall.")
    tanwin_w = FloatSlider(
        value=40.0, min=5.0, max=120.0, step=5.0, label="tangent_window_um",
        tooltip="Arc-length window the centerline tangent is smoothed over "
                "(Gaussian, sigma = window/4). This sets the scale of 'along "
                "the vessel'. A short window tracks medial-axis wobble; the "
                "graph's own kappa_window_um is 40 µm.")
    compute_btn = PushButton(
        text="Compute orientation",
        tooltip="Measure every adhesion's axis and resolve it in the vessel's "
                "frame. Needs the adhesion instances, the vessel mesh (with "
                "normals) and the centerline graph.")

    def _compute(after=None):
        labels = _adhesion_labels(state, viewer)
        if labels is None:
            print("Orientation: no 'Adhesion instances' layer — run the "
                  "Adhesion RF tab (or Load all saved) first.", flush=True)
            call_all(after, label="Adhesion orientation")
            return
        cache = state.get("vessel_mesh_cache") or {}
        verts, normals = cache.get("verts"), cache.get("normals")
        if verts is None:
            print("Orientation: build or load the vessel mesh first "
                  "(Vessel mesh tab).", flush=True)
            call_all(after, label="Adhesion orientation")
            return
        if normals is None:
            # The *sign* of the normal is irrelevant here — an adhesion axis has
            # no sign either, so `into_deg` uses |e1 . n| — which is why this
            # does NOT need `nuc_side.mesh_sign_is_luminal`. But a mesh with no
            # stored normals has no wall direction at all, and re-deriving one
            # from the faces is exactly the enclosed-volume rule that fails on
            # the open Poisson surface. So: refuse, don't invent.
            print("Orientation: this mesh carries no per-vertex normals, so "
                  "there is no wall direction to resolve against. Rebuild with "
                  "a Poisson-family method.", flush=True)
            call_all(after, label="Adhesion orientation")
            return
        graph = state.get("vessel_graph")
        if not graph or not len(graph.get("points_zyx_um", [])):
            print("Orientation: compute the centerline graph first "
                  "(Vessel mesh tab → 3. Centerline graph).", flush=True)
            call_all(after, label="Adhesion orientation")
            return

        scale = state.get("scale", (1.0, 1.0, 1.0))
        elong_min = float(elong_w.value)
        max_dist = float(maxdist_w.value)
        win = float(tanwin_w.value)
        build_id = int(state.get("vessel_build_id", 0))

        @thread_worker
        def _run():
            ax = ao.adhesion_axes(np.asarray(labels), scale)
            tan = ao.centerline_tangents(graph["points_zyx_um"],
                                         graph["s_um"], graph["edge_id"], win)
            fr = ao.wall_frames(ax["centroid_zyx_um"], verts, normals,
                                graph["points_zyx_um"], tan)
            res = ao.resolve(ax["e1"], fr["n"], fr["a"], fr["c"])
            gate = (~fr["degenerate"]) & np.isfinite(ax["elong"]) \
                & (ax["elong"] >= elong_min) \
                & (fr["dist_mesh_um"] <= max_dist)
            out = {**ax, **fr, **res, "gate": gate,
                   "e1_z_deg": ao.angle_to_z(ax["e1"]),
                   "elong_min": elong_min, "max_dist_um": max_dist,
                   "tangent_window_um": win, "build_id": build_id,
                   "n_total": int(len(gate)), "n_gated": int(gate.sum()),
                   "n_degenerate": int(fr["degenerate"].sum())}
            return out

        @with_after
        def _done(out):
            state["adhesion_orient"] = out
            _refresh_status()
            g = out["gate"]
            if not g.any():
                print("Orientation: the gate keeps nothing — lower elong_min "
                      "or raise max_dist_um.", flush=True)
                return
            # The lab-frame component, printed every time rather than kept for
            # whoever thinks to look: it is the single number that says how much
            # of what follows is the microscope.
            print(f"Orientation: {int(g.sum()):,}/{len(g):,} adhesions gated "
                  f"(elong >= {elong_min:g}, within {max_dist:g} µm of wall). "
                  f"Median angle of the axes to the stack z axis: "
                  f"{np.median(out['e1_z_deg'][g]):.1f}° — the PSF is elongated "
                  f"in z, so a small number here means the axes are largely "
                  f"lab-frame and every vessel-frame angle inherits that.",
                  flush=True)
            R, pref = ao.axial_resultant(out["inwall_deg"][g],
                                         out["inwall_mag"][g])
            ctl = ao.lab_axis_control(out["n"][g], out["a"][g], out["c"][g])
            print(f"  In-wall preferred direction {pref:.1f}° "
                  f"(0 = along the tube, 90 = around it), resultant "
                  f"R = {R:.3f}.", flush=True)
            # Printed together, always. Projecting a z-leaning axis into the
            # wall leaves a residual pointing where the WALL tilts, and on a
            # tube that is circumferential — so this statistic reads
            # "circumferential" on data containing no adhesions at all. On
            # 04step_003 the control beats the measurement.
            print(f"  CONTROL — the same statistic with every axis replaced by "
                  f"a constant stack-z axis, which carries no information about "
                  f"any adhesion: R = {ctl['R']:.3f} at {ctl['pref_deg']:.1f}°. "
                  f"The in-wall angle is only a measurement where it beats "
                  f"THIS, not where it beats zero. Use group 3 in the xy plane "
                  f"instead — a constant z axis has no xy component, so that "
                  f"view is immune by construction.", flush=True)

        w = _run()
        w.returned.connect(with_after(_done, after,
                                      label="Adhesion orientation"))
        print("Measuring adhesion orientation ...", flush=True)
        w.start()

    compute_btn.changed.connect(lambda *_: _compute())
    state["run_adhesion_orient"] = _compute

    g1 = collapsible_group("1. Orientation",
                           [orient_status_w, elong_w, maxdist_w, tanwin_w,
                            compute_btn], expanded=True)

    # ---------------------------------------------------------------- group 2
    axis_len_w = FloatSlider(
        value=4.0, min=0.5, max=20.0, step=0.5, label="arrow_um",
        tooltip="Drawn length of each adhesion's axis, in µm. Display only — "
                "it does not scale with the adhesion.")
    axis_draw_w = ComboBox(
        choices=list(AXIS_DRAW.keys()),
        value="xy component only (what is measurable)", label="draw",
        tooltip="The full 3D axis points almost along z for most adhesions "
                "(median |e1_z| = 0.97), so looking down z it points at the "
                "camera and its visible stub is noise — the arrows then look "
                "unrelated to the adhesions. Drawing the xy component shows "
                "the direction that is actually measurable, and it is the same "
                "quantity the binned glyphs average.")
    min_xy_w = FloatSlider(
        value=0.20, min=0.0, max=0.9, step=0.05, label="min_xy",
        tooltip="In xy mode, skip adhesions whose axis has less of an xy "
                "component than this. Their xy direction is a small residual "
                "of a z-pointing axis, i.e. noise — NOT DRAWN rather than "
                "drawn as a confident arrow.")
    axis_color_w = ComboBox(
        choices=list(AXIS_COLOR_BY.keys()), value="local coherence",
        label="colour by",
        tooltip="⚑ marks a quantity dominated by the stack's z axis rather "
                "than by the vessel. Those are offered because hiding them "
                "would hide the confound, not because they are safe to read.")
    coh_radius_w = FloatSlider(
        value=10.0, min=1.0, max=60.0, step=1.0, label="coherence_radius_um",
        tooltip="Neighbourhood radius for 'local coherence' — the mean "
                "|axis · axis| with every adhesion within this distance. The "
                "measured alignment decays over tens of µm, so this is the "
                "scale you are asking about.")
    axes_btn = PushButton(
        text="Show per-adhesion axes",
        tooltip="One arrow per gated adhesion. This view is easy to "
                "over-read: local agreement looks like alignment with the "
                "vessel. Group 3 is the one that separates those.")

    def _plane_of(widget):
        return PLANES[widget.value]

    axis_plane_w = ComboBox(choices=list(PLANES.keys()),
                            value="xy only (immune to the z-elongated PSF)",
                            label="coherence plane",
                            tooltip="Which components the local coherence is "
                                    "measured on. In xy the PSF imposes no "
                                    "preferred direction, so an excess there "
                                    "cannot be the microscope.")

    def _show_axes():
        r = _result()
        if not r:
            print("Show axes: compute the orientation first.", flush=True)
            return
        g = r["gate"]
        if not g.any():
            print("Show axes: the gate keeps nothing.", flush=True)
            return
        P, e1 = r["centroid_zyx_um"][g], r["e1"][g]
        L = float(axis_len_w.value)
        mode = AXIS_DRAW[axis_draw_w.value]
        if mode == "xy":
            # Flattened into the viewing plane and renormalised, so the arrow
            # has a visible length and shows the component that is measurable.
            # `_project` returns NaN where nothing survives the flattening, and
            # those are dropped by `vector_mask` below rather than renormalised
            # up out of noise.
            drawn, xymag = ao._project(e1, "xy")
            drawn = np.where((xymag >= float(min_xy_w.value))[:, None],
                             drawn, np.nan)
        else:
            drawn, xymag = e1, np.linalg.norm(e1[:, 1:], axis=1)
        # Drawn centred on the adhesion, not starting at it: the axis has no
        # sign, and an arrow drawn from the centroid would assert one.
        vec = np.stack([P - 0.5 * L * drawn, L * drawn], axis=1)
        key = AXIS_COLOR_BY[axis_color_w.value]
        if key == "local_coh":
            vals = ao.local_coherence(P, e1, float(coh_radius_w.value),
                                      plane=_plane_of(axis_plane_w))["coherence"]
        else:
            vals = np.asarray(r[key], float)[g]
        # One mask for the arrow and its colour, so a dropped row cannot shift
        # the colouring by one. An adhesion with no measurable colour value —
        # local coherence with too few neighbours, say — is NOT DRAWN rather
        # than drawn in some default shade, the same rule the wall-shear layer
        # follows for un-measured wall ([P49](docs/ISSUES.md#p49)).
        ok = vector_mask(vec, vals)
        if not ok.any():
            print(f"Show axes: nothing drawable — all {len(vec):,} axes have "
                  f"no finite {axis_color_w.value}.", flush=True)
            return
        vec, vals = vec[ok], np.asarray(vals, float)[ok]
        if AXES_LAYER in viewer.layers:
            viewer.layers.remove(AXES_LAYER)
        viewer.add_vectors(vec, name=AXES_LAYER, edge_width=0.4,
                           blending="translucent", features={"value": vals},
                           edge_color="value", edge_colormap="viridis")
        dropped = int((~ok).sum())
        print(f"{AXES_LAYER}: {len(vec):,} axes drawn as "
              f"{axis_draw_w.value}, coloured by {axis_color_w.value} "
              f"(median {np.median(vals):.3g})"
              + (f"; {dropped:,} not drawn." if dropped else "."), flush=True)
        if mode == "3d":
            # Said out loud, because an arrow pointing at the camera looks like
            # a short arrow pointing somewhere arbitrary, not like a long one
            # pointing at you.
            print(f"  ⚑ median |e1_xy| = {np.median(xymag):.3f}, so looking "
                  f"down z a {L:g} µm arrow shows a median of "
                  f"{L * np.median(xymag):.2f} µm in plane and the rest points "
                  f"at the camera. Expect these NOT to line up with the "
                  f"adhesions as drawn in a top-down view — switch 'draw' to "
                  f"the xy component.", flush=True)

    axes_btn.changed.connect(lambda *_: _show_axes())
    g2 = collapsible_group("2. Per-adhesion axes",
                           [axis_draw_w, min_xy_w, axis_len_w, axis_color_w,
                            axis_plane_w, coh_radius_w, axes_btn])

    # ---------------------------------------------------------------- group 3
    bin_w = FloatSlider(
        value=20.0, min=2.0, max=100.0, step=1.0, label="bin_um",
        tooltip="Patch size. The neighbour alignment decays over tens of µm, "
                "so a much larger bin averages the structure away and a much "
                "smaller one has too few adhesions to average. There is no "
                "right value — move it and see which features survive.")
    mincount_w = Slider(
        value=5, min=2, max=60, label="min_count",
        tooltip="Patches with fewer adhesions than this are NOT DRAWN, rather "
                "than drawn from too little data. Same rule the wall-shear "
                "layer follows for un-measured wall.")
    mincoh_w = FloatSlider(
        value=0.0, min=0.0, max=1.0, step=0.02, label="min_coherence",
        tooltip="Hide patches whose adhesions do not agree on a direction. At "
                "0 every patch above min_count is drawn, including ones whose "
                "mean axis points at nothing in particular — which is worth "
                "seeing before you filter them out.")
    glyph_plane_w = ComboBox(choices=list(PLANES.keys()),
                             value="xy only (immune to the z-elongated PSF)",
                             label="plane",
                             tooltip="Which components the patch average is "
                                     "taken on. Start in xy: an alignment "
                                     "there cannot have been produced by the "
                                     "z-elongated PSF.")
    glyph_color_w = ComboBox(
        choices=list(GLYPH_COLOR_BY.keys()),
        value="axial_cos (1 = patch runs along the vessel)", label="colour by",
        tooltip="axial_cos is the direct test of 'do they run along the "
                "vessel': the cosine between the patch's mean axis and the "
                "patch's mean vessel axis. Compare it against the isotropic "
                "level printed in the log, not against zero.")
    glyph_len_w = FloatSlider(
        value=12.0, min=1.0, max=60.0, step=1.0, label="glyph_um",
        tooltip="Drawn length of each patch's mean axis, in µm. Display only.")
    glyph_btn = PushButton(
        text="Show binned glyphs",
        tooltip="One arrow per patch of vessel: the mean axis of the "
                "adhesions in it. The view this tab exists for.")

    def _show_glyphs():
        r = _result()
        if not r:
            print("Show glyphs: compute the orientation first.", flush=True)
            return
        g = r["gate"]
        if not g.any():
            print("Show glyphs: the gate keeps nothing.", flush=True)
            return
        plane = _plane_of(glyph_plane_w)
        gl = ao.bin_glyphs(r["centroid_zyx_um"][g], r["e1"][g],
                           float(bin_w.value), min_count=int(mincount_w.value),
                           plane=plane, axial=r["a"][g])
        keep = np.isfinite(gl["coherence"]) & (gl["coherence"]
                                               >= float(mincoh_w.value))
        n_hidden = int((~keep).sum())
        if not keep.any():
            print(f"Show glyphs: no patch survives (bin {bin_w.value:g} µm, "
                  f"min_count {int(mincount_w.value)}, min_coherence "
                  f"{mincoh_w.value:g}); {gl['n_bins_dropped']} patches were "
                  "already too sparse to average.", flush=True)
            if GLYPH_LAYER in viewer.layers:
                viewer.layers.remove(GLYPH_LAYER)
            _glyph_status(gl, keep, plane)
            return
        C, A = gl["center_zyx_um"][keep], gl["axis"][keep]
        L = float(glyph_len_w.value)
        vec = np.stack([C - 0.5 * L * A, L * A], axis=1)
        vals = np.asarray(gl[GLYPH_COLOR_BY[glyph_color_w.value]], float)[keep]
        # Same one-mask rule as the per-adhesion axes. `axial_cos` is NaN for a
        # patch with no usable vessel direction, and such a patch must not be
        # drawn in the colour that happens to sit at the bottom of the map.
        ok = vector_mask(vec, vals)
        if not ok.any():
            print(f"Show glyphs: {int(keep.sum())} patches, none with a finite "
                  f"{glyph_color_w.value} to colour by.", flush=True)
            if GLYPH_LAYER in viewer.layers:
                viewer.layers.remove(GLYPH_LAYER)
            return
        if GLYPH_LAYER in viewer.layers:
            viewer.layers.remove(GLYPH_LAYER)
        viewer.add_vectors(vec[ok], name=GLYPH_LAYER, edge_width=1.2,
                           blending="translucent",
                           features={"value": vals[ok]}, edge_color="value",
                           edge_colormap="magma")
        keep_drawn = keep.copy()
        keep_drawn[np.flatnonzero(keep)[~ok]] = False
        _glyph_status(gl, keep_drawn, plane,
                      n_hidden=n_hidden + int((~ok).sum()))

    def _glyph_status(gl, keep, plane, n_hidden=0):
        n = int(keep.sum())
        if n == 0:
            glyph_status_w.value = "3. Glyphs — nothing drawn"
            return
        coh = gl["coherence"][keep]
        res = ao.shuffled_axial_cos(gl["axis"][keep], gl["vessel_axis"][keep])
        med, (lo, hi) = res["observed"], res["null"]
        verdict = ("ALONG the vessel" if med > hi else
                   "ACROSS the vessel" if med < lo else
                   "NOT aligned with the vessel")
        glyph_status_w.value = (
            f"3. Glyphs — {n} patches at {gl['bin_um']:g} µm "
            f"({gl['n_bins_dropped']} too sparse, {n_hidden} not drawn), "
            f"coherence {np.median(coh):.2f}, {verdict}")
        # The comparison is the whole point of the display, so it is printed
        # rather than left for the reader to remember — and it is printed
        # against the SHUFFLED baseline, not the isotropic constant, which this
        # data beats for reasons that have nothing to do with the patches
        # (see `ao.shuffled_axial_cos`).
        print(f"{GLYPH_LAYER}: {n} patches at {gl['bin_um']:g} µm, plane "
              f"{plane}. Median coherence {np.median(coh):.3f} — how much the "
              f"adhesions inside a patch agree.", flush=True)
        print(f"  Median axial_cos {med:.3f}, against [{lo:.3f}, {hi:.3f}] for "
              f"the same patch axes paired with the WRONG patches' vessel "
              f"directions. Verdict: the patches are {verdict}.", flush=True)
        print(f"  (For context only, two uniformly random axes would give "
              f"{ISOTROPIC_AXIAL_COS[plane]:.3f} — do not read axial_cos "
              f"against that: the vessel runs one general direction across the "
              f"field, so any global lean in the axes beats it.)", flush=True)

    glyph_btn.changed.connect(lambda *_: _show_glyphs())
    g3 = collapsible_group("3. Binned average glyphs",
                           [glyph_status_w, bin_w, mincount_w, mincoh_w,
                            glyph_plane_w, glyph_color_w, glyph_len_w,
                            glyph_btn], expanded=True)

    # ---------------------------------------------------------------- group 4
    lc_btn = PushButton(
        text="Show local coherence (points)",
        tooltip="Each adhesion coloured by how well its neighbours agree with "
                "it, using the radius and plane set in group 2. No binning, so "
                "a feature that shows here and not in group 3 is a bin "
                "artifact — and one that shows in group 3 and not here is too.")

    def _show_coherence():
        r = _result()
        if not r:
            print("Local coherence: compute the orientation first.", flush=True)
            return
        g = r["gate"]
        if not g.any():
            print("Local coherence: the gate keeps nothing.", flush=True)
            return
        P, e1 = r["centroid_zyx_um"][g], r["e1"][g]
        plane = _plane_of(axis_plane_w)
        lc = ao.local_coherence(P, e1, float(coh_radius_w.value), plane=plane)
        v = lc["coherence"]
        ok = np.isfinite(v)
        if not ok.any():
            print(f"Local coherence: no adhesion has neighbours within "
                  f"{coh_radius_w.value:g} µm.", flush=True)
            return
        if COH_LAYER in viewer.layers:
            viewer.layers.remove(COH_LAYER)
        viewer.add_points(
            P[ok], name=COH_LAYER, size=3.0, features={"value": v[ok]},
            face_color="value", face_colormap="viridis",
            blending="translucent")
        # The level to compare against is NOT 0.5 in 3D — with every axis
        # leaning toward z the no-structure level is far higher — so the number
        # is printed with that warning attached rather than bare.
        extra = ("" if plane == "xy" else
                 "  NOTE: in 3D the no-structure level is well above 0.5 "
                 "because the axes are globally z-aligned; measure in xy to "
                 "remove that.")
        print(f"{COH_LAYER}: {int(ok.sum()):,} adhesions with neighbours "
              f"within {coh_radius_w.value:g} µm, plane {plane}, median "
              f"coherence {np.median(v[ok]):.3f}." + extra, flush=True)

    lc_btn.changed.connect(lambda *_: _show_coherence())
    g4 = collapsible_group("4. Local coherence", [lc_btn])

    _refresh_status()
    return [g1, g2, g3, g4]
