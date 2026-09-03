"""Open OME-TIFF stacks in napari with raw channels + FRET efficiency/force.

Usage:
    python -m pipeline.view [--project-dir DIR] [--params P] [--lut L]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import napari
from qtpy.QtCore import Qt
from napari.settings import get_settings
from napari.utils.triangulation_backend import TriangulationBackend
from qtpy.QtWidgets import (
    QScrollArea, QTabWidget, QWidget, QVBoxLayout)

from .gui import (
    adhesion_rf,
    cadherin,
    cadherin_rf,
    cell_mesh,
    clover,
    display,
    loader,
    logpanel,
    measure,
    nuclei,
    roi,
    vessel_mesh,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROJECT_DIR = ROOT / "3d_images_folder"
DEFAULT_PARAMS = ROOT / "3d_images_folder" / "params.json"
DEFAULT_LUT = ROOT / "Fret_LUT" / "Clover_mRuby2_GGSGGS7_force_new.txt"


def _wrap(parts) -> QWidget:
    """Stack one or more widgets (magicgui Container or Qt) into a **scrolling**
    QWidget.

    The scroll area is the whole point. A dock widget is only as tall as the
    napari window, and the tabs have outgrown that — the Vessel mesh tab alone
    carries the mesh controls, the curvature controls, the graph controls and
    the shear controls, and its save/load buttons were being cut off the bottom
    with no way to reach them. Without a scroll area a Qt layout does not
    shorten its children; it simply clips whatever does not fit, silently.

    Applied to **every** tab rather than to the one that overflowed, because
    which tab overflows depends on the screen and on how many controls a tab has
    grown since anyone last looked. Horizontal scrolling stays off and the inner
    widget is made resizable, so the widgets still track the dock's width and
    only the vertical axis scrolls.
    """
    if not isinstance(parts, (list, tuple)):
        parts = [parts]
    inner = QWidget()
    lay = QVBoxLayout(inner)
    lay.setContentsMargins(4, 4, 4, 4)
    for p in parts:
        lay.addWidget(p.native if hasattr(p, "native") else p)
    lay.addStretch(1)

    area = QScrollArea()
    area.setWidget(inner)
    area.setWidgetResizable(True)
    area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    area.setFrameShape(QScrollArea.NoFrame)
    return area


def _disable_quickedit() -> None:
    """Disable Windows console QuickEdit mode.

    QuickEdit pauses the process when the user clicks the console window
    (e.g. after pressing a napari button), requiring Enter to resume.
    """
    if sys.platform != "win32":
        return
    try:
        import ctypes
        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(-10)  # STD_INPUT_HANDLE
        mode = ctypes.c_ulong()
        if kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            ENABLE_QUICK_EDIT = 0x0040
            ENABLE_EXTENDED_FLAGS = 0x0080
            mode.value = (mode.value & ~ENABLE_QUICK_EDIT) | ENABLE_EXTENDED_FLAGS
            kernel32.SetConsoleMode(handle, mode)
    except Exception:
        pass


def main() -> None:
    _disable_quickedit()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-dir", type=Path, default=DEFAULT_PROJECT_DIR)
    ap.add_argument("--params", type=Path, default=DEFAULT_PARAMS)
    ap.add_argument("--lut", type=Path, default=DEFAULT_LUT)
    args = ap.parse_args()

    # Force pure-Python shape triangulation. napari 0.6.0 defaults to
    # "fastest available", which grabs the compiled `bermuda` backend when it
    # is importable — but the bermuda build on the Windows env is API-
    # incompatible (no `triangulate_path_edge`), crashing the "Cell mesh
    # (slice)" line overlay. The dev env has no bermuda and silently falls back
    # to pure-Python; pinning it here makes both envs behave the same. These
    # overlays are tiny, so there is no performance cost.
    get_settings().experimental.triangulation_backend = TriangulationBackend.pure_python

    viewer = napari.Viewer(ndisplay=3)
    logpanel.install(viewer)  # mirror stdout/stderr into an in-napari Log dock

    def _free_gpu():
        try:
            import cupy as cp
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
            print("GPU memory freed.", flush=True)
        except Exception:
            pass

    viewer.window._qt_window.destroyed.connect(_free_gpu)
    state: dict = {
        "dapi_raw": None, "cadh_raw": None, "clover_raw": None,
        "force_vol": None, "adhesion_rf": None, "adhesion_proba": None,
        "scale": (1.0, 1.0, 1.0),
        "dapi_mask_layer": None, "cadh_mask_layer": None,
        "clover_mask_layer": None,
        "dapi_state": {"labels": None, "sizes": None, "nuclei": None,
                       "nuclei_sizes": None, "nuclei_ids": None,
                       "nuclei_is_fib": None, "nuclei_min_size_vox": None},
        "cadh_state": {
            "labels": None, "sizes": None,
            "centroids": None, "normals": None, "label_ids": None,
            "keep_mask": None,
        },
        "clover_state": {"labels": None, "sizes": None, "method": None},
        "nuc_cache": {"labels_id": None, "centroids": None, "normals": None, "lids": None},
        "mesh_cache": {"offset_centroids": None, "outward": None,
                       "inward_offset_um": 0.0, "kdtree": None},
        "vessel_mesh_cache": {"verts": None, "faces": None},
        "raw_layers": [],
        "hists": [],
    }

    tabs = QTabWidget()
    tabs.addTab(_wrap([loader.build(state, viewer, args.project_dir,
                                    args.params, args.lut)]), "Load")
    tabs.addTab(_wrap([display.build(state, viewer)]), "Display")
    tabs.addTab(_wrap(list(nuclei.build(state, viewer))), "Nuclei")
    tabs.addTab(_wrap(list(cadherin.build(state, viewer))), "Cadherin")
    tabs.addTab(_wrap([cadherin_rf.build(state, viewer)]), "VE-cad RF")
    tabs.addTab(_wrap(list(cell_mesh.build(state, viewer))), "Cell mesh")
    tabs.addTab(_wrap(list(vessel_mesh.build(state, viewer))), "Vessel mesh")
    tabs.addTab(_wrap([clover.build(state, viewer)]), "Clover")
    tabs.addTab(_wrap([adhesion_rf.build(state, viewer)]), "Adhesion RF")
    tabs.addTab(_wrap([roi.build(state, viewer)]), "ROI")
    tabs.addTab(_wrap(list(measure.build(state, viewer))), "Measure")

    # Selecting the Vessel mesh tab fills its alpha3d_um from the circumradius
    # histogram of the current cadherin mask, so the slider shows the number the
    # build will use instead of a stale one from another image. It is a
    # suggestion, not a lock: drag the slider and it stays where you put it
    # until the next image loads. No-ops when there is no mask yet, when the
    # method does not read alpha, or when the value has been set by hand.
    def _on_tab(idx: int) -> None:
        if tabs.tabText(idx) != "Vessel mesh":
            return
        fn = state.get("vessel_auto_alpha")
        if fn is None:
            return
        try:
            fn()
        except Exception as exc:
            print(f"  (auto alpha skipped: {type(exc).__name__}: {exc})",
                  flush=True)

    tabs.currentChanged.connect(_on_tab)

    viewer.window.add_dock_widget(tabs, area="right", name="Pipeline")
    napari.run()


if __name__ == "__main__":
    main()
