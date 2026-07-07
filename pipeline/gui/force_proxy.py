"""Uncalibrated ratio "Force proxy" tab: numerator/denominator role ÷ role.

Not the calibrated FRET `Force (pN)` layer in `fret.py` — just a live
brighter-means-more-force ratio for datasets that lack real 3-channel FRET
calibration (e.g. `fret` present only as a dummy placeholder). Default
clover/ruby, but any two loaded roles can be picked.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container

from ..fret import background_subtract
from .widgets import debounce, percentile_limits, pyramid_views

DEFAULT_NUMERATOR = "clover"
DEFAULT_DENOMINATOR = "ruby"


def _background_subtract_volume(vol: np.ndarray) -> np.ndarray:
    out = np.empty(vol.shape, dtype=np.float64)

    def _do(z: int) -> None:
        out[z] = background_subtract(vol[z])

    with ThreadPoolExecutor(max_workers=min(vol.shape[0], 16)) as ex:
        list(ex.map(_do, range(vol.shape[0])))
    return out


def build(state, viewer) -> Container:

    @magicgui(
        auto_call=True,
        numerator={"widget_type": "ComboBox", "choices": [""]},
        denominator={"widget_type": "ComboBox", "choices": [""]},
        subtract_bg={"widget_type": "CheckBox"},
        floor={"widget_type": "FloatSlider", "min": 1.0, "max": 500.0},
    )
    def ratio_controls(
        numerator: str = "",
        denominator: str = "",
        subtract_bg: bool = True,
        floor: float = 50.0,
    ):
        _apply(numerator, denominator, subtract_bg, floor)

    @debounce(200)
    def _apply(numerator, denominator, subtract_bg, floor):
        raw_by_role = state.get("raw_by_role") or {}
        if not numerator or not denominator:
            print("Force proxy: select numerator and denominator roles.",
                  flush=True)
            return
        if numerator == denominator:
            print("Force proxy: numerator and denominator must differ.",
                  flush=True)
            return
        num_raw = raw_by_role.get(numerator)
        den_raw = raw_by_role.get(denominator)
        if num_raw is None or den_raw is None:
            print("Force proxy: selected role not loaded.", flush=True)
            return

        if subtract_bg:
            num = _background_subtract_volume(num_raw)
            den = _background_subtract_volume(den_raw)
        else:
            num = num_raw.astype(np.float64)
            den = den_raw.astype(np.float64)

        ratio = np.empty(num.shape, dtype=np.float32)
        np.divide(num, np.maximum(den, floor), out=ratio)

        name = f"Force proxy ({numerator}/{denominator})"
        levels = state.get("levels", 1)
        scale = state.get("scale", (1.0, 1.0, 1.0))
        clim = percentile_limits(ratio)

        layer = state.get("force_proxy_layer")
        if layer is None or layer not in viewer.layers:
            layer = viewer.add_image(
                pyramid_views(ratio, levels), name=name, colormap="turbo",
                scale=scale, contrast_limits=clim, rendering="mip",
                blending="additive", visible=True,
                multiscale=levels > 1,
            )
            state["force_proxy_layer"] = layer
        else:
            layer.data = pyramid_views(ratio, levels)
            layer.name = name
            layer.contrast_limits_range = (float(ratio.min()), float(ratio.max()) + 1.0)
            layer.contrast_limits = clim
            layer.visible = True

    def _seed_choices():
        raw_by_role = state.get("raw_by_role") or {}
        roles = sorted(raw_by_role.keys())
        ratio_controls.numerator.choices = roles
        ratio_controls.denominator.choices = roles
        if len(roles) < 2:
            print("Force proxy: fewer than 2 roles loaded, skipping.",
                  flush=True)
            return
        if DEFAULT_NUMERATOR in roles and DEFAULT_DENOMINATOR in roles:
            num, den = DEFAULT_NUMERATOR, DEFAULT_DENOMINATOR
        else:
            num, den = roles[0], roles[1]
        ratio_controls.numerator.value = num
        ratio_controls.denominator.value = den
        ratio_controls()

    state["seed_ratio_choices"] = _seed_choices

    return Container(widgets=[ratio_controls], labels=True)
