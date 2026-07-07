"""Open OME-TIFF stacks in napari with raw channels + FRET efficiency/force.

Usage:
    python -m pipeline.view [--project-dir DIR] [--params P] [--lut L]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import napari
from napari.settings import get_settings
from napari.utils.triangulation_backend import TriangulationBackend
from qtpy.QtWidgets import QTabWidget, QWidget, QVBoxLayout

from .gui import (
    adhesion_rf,
    cadherin,
    cell_mesh,
    clover,
    display,
    force_proxy,
    loader,
    measure,
    nuclei,
    vessel_mesh,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROJECT_DIR = ROOT / "3d_images_folder"
DEFAULT_PARAMS = ROOT / "3d_images_folder" / "params.json"
DEFAULT_LUT = ROOT / "Fret_LUT" / "Clover_mRuby2_GGSGGS7_force_new.txt"


def _wrap(parts) -> QWidget:
    """Stack one or more widgets (magicgui Container or Qt) into a QWidget."""
    if not isinstance(parts, (list, tuple)):
        parts = [parts]
    w = QWidget()
    lay = QVBoxLayout(w)
    lay.setContentsMargins(4, 4, 4, 4)
    for p in parts:
        lay.addWidget(p.native if hasattr(p, "native") else p)
    lay.addStretch(1)
    return w


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
        "dapi_state": {"labels": None, "sizes": None},
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
        "raw_by_role": {}, "force_proxy_layer": None,
    }

    tabs = QTabWidget()
    tabs.addTab(_wrap([loader.build(state, viewer, args.project_dir,
                                    args.params, args.lut)]), "Load")
    tabs.addTab(_wrap([display.build(state, viewer)]), "Display")
    tabs.addTab(_wrap([force_proxy.build(state, viewer)]), "Force proxy")
    tabs.addTab(_wrap(list(nuclei.build(state, viewer))), "Nuclei")
    tabs.addTab(_wrap(list(cadherin.build(state, viewer))), "Cadherin")
    tabs.addTab(_wrap(list(cell_mesh.build(state, viewer))), "Cell mesh")
    tabs.addTab(_wrap(list(vessel_mesh.build(state, viewer))), "Vessel mesh")
    tabs.addTab(_wrap([clover.build(state, viewer)]), "Clover")
    tabs.addTab(_wrap([adhesion_rf.build(state, viewer)]), "Adhesion RF")
    tabs.addTab(_wrap(list(measure.build(state, viewer))), "Measure")

    viewer.window.add_dock_widget(tabs, area="right", name="Pipeline")
    napari.run()


if __name__ == "__main__":
    main()
