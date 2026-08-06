"""Region-of-interest (X-Y box) tab.

Draw a rectangle in the imaging plane to restrict Measure to adhesions inside it
— e.g. to exclude curved vessel segments near the field edges from a "straight"
image. The box is applied by X-Y only (all Z), and Measure keeps an adhesion if
its centroid falls inside. Clearing the ROI restores whole-field measurement.
"""

from __future__ import annotations

import numpy as np
from magicgui.widgets import Container, Label, PushButton

ROI_LAYER = "ROI (X-Y)"


def build(state, viewer) -> Container:
    status = Label(value="No ROI — Measure uses the whole field.")

    def _img_shape():
        for key in ("clover_raw", "dapi_raw", "cadh_raw", "ruby_raw"):
            arr = state.get(key)
            if arr is not None:
                return np.asarray(arr).shape
        return None

    def _add_box(*_):
        shape = _img_shape()
        if shape is None:
            print("ROI: load an image first.", flush=True)
            return
        nz, ny, nx = shape
        scale = state.get("scale", (1.0, 1.0, 1.0))
        try:
            zc = int(viewer.dims.current_step[0])
        except Exception:
            zc = nz // 2
        # Default box = central 60% of the X-Y field; user moves/resizes it.
        y0, y1 = int(ny * 0.2), int(ny * 0.8)
        x0, x1 = int(nx * 0.2), int(nx * 0.8)
        rect = np.array([[zc, y0, x0], [zc, y0, x1],
                         [zc, y1, x1], [zc, y1, x0]], dtype=float)
        if ROI_LAYER in viewer.layers:
            viewer.layers.remove(ROI_LAYER)
        lyr = viewer.add_shapes(
            [rect], shape_type="rectangle", name=ROI_LAYER, scale=scale,
            edge_color="#ffff00ff", face_color="#ffff0022", edge_width=3.0,
        )
        try:
            lyr.mode = "select"
        except Exception:
            pass
        print("ROI: added a default box. Switch to the 2D X-Y view, move/resize "
              "it over the region you want, then click 'Apply ROI'. For a "
              "non-rectangular region, use the polygon tool in the ROI layer's "
              "controls (top-left) to click out any outline (or draw several "
              "shapes) — Apply keeps adhesions inside any of them.", flush=True)

    def _apply(*_):
        if ROI_LAYER not in viewer.layers:
            print("ROI: add or draw a shape first.", flush=True)
            return
        data = viewer.layers[ROI_LAYER].data
        if not len(data):
            print("ROI: no shapes in the ROI layer.", flush=True)
            return
        scale = state.get("scale", (1.0, 1.0, 1.0))
        # Store each shape's outline (Y, X in µm) as a polygon. Rectangles and
        # arbitrary polygons alike become vertex lists; Measure does a
        # point-in-polygon test, so the region can be any shape.
        polys = []
        for s in data:
            s = np.asarray(s)
            yx = s[:, -2:] * np.array([float(scale[1]), float(scale[2])])
            if len(yx) >= 3:
                polys.append(yx.astype(float))
        if not polys:
            print("ROI: shapes need at least 3 vertices.", flush=True)
            return
        state["roi_polygons_yx_um"] = polys
        state["roi_yx_um"] = None
        allv = np.vstack(polys)
        status.value = (f"ROI: {len(polys)} shape(s), "
                        f"Y [{allv[:,0].min():.0f}, {allv[:,0].max():.0f}] "
                        f"X [{allv[:,1].min():.0f}, {allv[:,1].max():.0f}] µm")
        print(f"ROI applied: {len(polys)} shape(s). Measure keeps adhesions "
              "whose centroid is inside any shape (any outline, not just a box).",
              flush=True)

    def _remove_box(*_):
        # Delete the starter rectangle(s) but keep any polygons you drew.
        if ROI_LAYER not in viewer.layers:
            print("ROI: no ROI layer.", flush=True)
            return
        lyr = viewer.layers[ROI_LAYER]
        types = list(getattr(lyr, "shape_type", []))
        rect_idx = {i for i, t in enumerate(types) if t == "rectangle"}
        if not rect_idx:
            print("ROI: no rectangle to remove.", flush=True)
            return
        lyr.selected_data = rect_idx
        lyr.remove_selected()
        print(f"ROI: removed {len(rect_idx)} starter box(es); kept "
              f"{len(lyr.data)} polygon(s). Re-click 'Apply ROI'.", flush=True)

    def _clear(*_):
        state["roi_polygons_yx_um"] = None
        state["roi_yx_um"] = None
        if ROI_LAYER in viewer.layers:
            viewer.layers.remove(ROI_LAYER)   # remove every drawn shape too
        status.value = "No ROI — Measure uses the whole field."
        print("ROI cleared — all ROI shapes removed; Measure uses the whole "
              "field.", flush=True)

    add_w = PushButton(text="Add ROI box (X-Y)")
    remove_box_w = PushButton(text="Remove starter box (keep polygons)")
    apply_w = PushButton(text="Apply ROI")
    clear_w = PushButton(text="Clear ROI")
    add_w.changed.connect(_add_box)
    remove_box_w.changed.connect(_remove_box)
    apply_w.changed.connect(_apply)
    clear_w.changed.connect(_clear)
    return Container(widgets=[add_w, remove_box_w, apply_w, clear_w, status],
                     labels=False)
