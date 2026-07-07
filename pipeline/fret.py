"""FRET math: background subtraction, bleedthrough correction, force lookup.

Per-slice (background_subtract, fret_correct) is verbatim from the old
napari_pipeline; compute_fret_volume is rewritten to take a role->C-index
map instead of the old channel_map dict.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np


def load_force_lut(path: str) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def background_subtract(image: np.ndarray, bin_size: int = 10) -> np.ndarray:
    img = image.astype(np.float64)
    positives = img[img >= 0]
    if len(positives) == 0:
        return np.clip(img, 0, None)
    binned = np.floor(positives / bin_size).astype(np.int64)
    values, counts = np.unique(binned, return_counts=True)
    mode_val = values[np.argmax(counts)] * bin_size
    return np.clip(img - mode_val, 0, None)


def fret_correct(
    idd: np.ndarray,
    dxam: np.ndarray,
    iaa: np.ndarray,
    abt: float,
    dbt: float,
    G: float,
    venus_thres: float,
    leave_neg: bool = False,
    lut_eff: np.ndarray | None = None,
    lut_force: np.ndarray | None = None,
) -> dict:
    idd = idd.astype(np.float64)
    dxam = dxam.astype(np.float64)
    iaa = iaa.astype(np.float64)

    over_mask = (idd >= 65535) | (dxam >= 65535) | (iaa >= 65535)
    under_mask = iaa < venus_thres
    exclude = over_mask | under_mask

    idd_w = np.where(exclude, 0.0, idd)
    dxam_w = np.where(exclude, 0.0, dxam)
    iaa_w = np.where(exclude, 0.0, iaa)

    fc = dxam_w - (abt * iaa_w) - (dbt * idd_w)

    fc_over_G = fc / G if G != 0 else fc
    denom = idd_w + fc_over_G
    eff = np.zeros_like(denom)
    np.divide(fc_over_G, denom, out=eff, where=(denom != 0))

    if not leave_neg:
        eff = np.clip(eff, 0, 1)
        fc = np.clip(fc, 0, None)

    force = np.zeros_like(eff)
    if lut_eff is not None and lut_force is not None:
        pos = eff > 0
        force[pos] = np.interp(eff[pos], lut_eff, lut_force)

    fc[exclude] = 0.0
    eff[exclude] = 0.0
    force[exclude] = 0.0
    return {"fc": fc, "eff": eff, "force": force}


def compute_fret_volume(
    data_zcyx: np.ndarray,
    role_to_c: dict,
    fret_params: dict,
    lut_eff: np.ndarray,
    lut_force: np.ndarray,
) -> dict:
    """Per-Z-slice background-subtract + FRET correction over a whole stack.

    role_to_c must map "clover", "fret", "ruby" -> channel index.
    """
    c_clover = role_to_c["clover"]
    c_fret = role_to_c["fret"]
    c_ruby = role_to_c["ruby"]

    Z, _, H, W = data_zcyx.shape
    fc_vol = np.empty((Z, H, W), dtype=np.float32)
    eff_vol = np.empty((Z, H, W), dtype=np.float32)
    force_vol = np.empty((Z, H, W), dtype=np.float32)

    def _do(z: int) -> None:
        idd = background_subtract(data_zcyx[z, c_clover])
        dxam = background_subtract(data_zcyx[z, c_fret])
        iaa = background_subtract(data_zcyx[z, c_ruby])
        result = fret_correct(
            idd, dxam, iaa,
            lut_eff=lut_eff,
            lut_force=lut_force,
            **fret_params,
        )
        fc_vol[z] = result["fc"].astype(np.float32)
        eff_vol[z] = result["eff"].astype(np.float32)
        force_vol[z] = result["force"].astype(np.float32)

    # Slices are independent; numpy ops release the GIL. Cap workers so the
    # float64 intermediates per slice don't blow memory on tall stacks.
    with ThreadPoolExecutor(max_workers=min(Z, 16)) as ex:
        list(ex.map(_do, range(Z)))

    return {"fc_vol": fc_vol, "eff_vol": eff_vol, "force_vol": force_vol}


def resample_fret_outputs(
    fret_dict: dict,
    src_zyx_um: tuple[float, float, float],
    target_um: float,
) -> dict:
    """Resample FRET output volumes (fc/eff/force) to an isotropic grid.

    Runs *after* compute_fret_volume so the non-linear ratio + LUT lookup
    operate on native-spacing data (no Jensen bias on the ratio). The
    resulting volumes are linearly upsampled for display alignment with
    the resampled raw stack.
    """
    from .io import resample_volume
    keys = ("fc_vol", "eff_vol", "force_vol")

    def _do(k: str):
        return k, resample_volume(fret_dict[k], src_zyx_um, target_um, order=1)

    with ThreadPoolExecutor(max_workers=3) as ex:
        return dict(ex.map(_do, keys))
