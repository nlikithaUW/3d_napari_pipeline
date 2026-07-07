"""Display tab — batch rendering + percentile auto-contrast for raw layers.

Deliberately NOT a duplicate of napari's built-in contrast slider. It adds the
two things napari lacks here:
  - percentile auto-contrast (napari's auto is min/max of a sparse sample,
    which collapses onto the noise floor for sparse-bright data); the Auto
    button writes into each layer's contrast_limits, so napari's native slider
    still does manual fine-tuning.
  - one rendering/attenuation control that batch-applies to every raw layer at
    once (napari only exposes rendering per-layer).
"""

from __future__ import annotations

from magicgui import magicgui
from magicgui.widgets import Container

from .widgets import percentile_limits

RENDER_MODES = ["mip", "attenuated_mip", "average", "translucent", "additive", "iso"]
FORCE_LAYER = "Force (pN)"


def build(state, viewer) -> Container:
    @magicgui(
        auto_call=True,
        rendering={"widget_type": "ComboBox", "choices": RENDER_MODES},
        attenuation={"widget_type": "FloatSlider", "min": 0.0, "max": 2.0, "step": 0.05},
    )
    def render_controls(rendering: str = "mip", attenuation: float = 1.0):
        layers = state.get("raw_layers") or []
        if not layers:
            print("Load a stack first.", flush=True)
            return
        for lyr in layers:
            lyr.rendering = rendering
            if rendering == "attenuated_mip":
                lyr.attenuation = float(attenuation)

    @magicgui(
        call_button="Auto contrast (percentile)",
        low_pct={"widget_type": "FloatSlider", "min": 0.0, "max": 50.0, "step": 0.5},
        high_pct={"widget_type": "FloatSlider", "min": 90.0, "max": 100.0, "step": 0.05},
    )
    def contrast_controls(low_pct: float = 1.0, high_pct: float = 99.9):
        layers = state.get("raw_layers") or []
        if not layers:
            print("Load a stack first.", flush=True)
            return
        for lyr in layers:
            if lyr.name == FORCE_LAYER:
                continue  # force is in pN; its fixed window is meaningful
            lo, hi = percentile_limits(lyr.data, low_pct, high_pct)
            lyr.contrast_limits_range = (0.0, 65535.0)
            lyr.contrast_limits = (lo, hi)
            print(f"  {lyr.name}: contrast [{lo:.0f}, {hi:.0f}] "
                  f"({low_pct}–{high_pct} pct)", flush=True)

    return Container(widgets=[render_controls, contrast_controls], labels=False)
