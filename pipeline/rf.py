"""Random-forest pixel classifier for 3-D adhesion segmentation.

Feature extraction runs on GPU (cupy) when available; falls back to CPU.
Each scale uploads the volume once, runs all filter calls on-device, extracts
values at the requested indices, then frees device memory before the next scale.

Label convention in annotation volumes:
  0 = unlabeled, 1 = adhesion, 2 = background
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as _sndi

from . import gpu as _gpu

if _gpu.HAS_GPU:
    import cupy as cp
    import cupyx.scipy.ndimage as _cndi
else:
    cp = _cndi = None

SCALES_UM = [0.3, 0.7, 1.0, 1.6, 3.5, 5.0]
N_PER_SCALE = 10   # Gaussian, grad_mag, LoG, DoG, 3 Hessian eigs, 3 ST eigs
N_FEATURES = 60    # 6 scales × 10

# Where the Gaussian kernel is cut off, in sigmas. scipy/cupy default to 4.0.
# Cost per axis is 2·ceil(truncate·sigma)+1 taps, so this is linear in runtime
# at the large scales (sigma 15.4 voxels at 5 µm on the 0.325 µm grid).
#
# **Changing it changes every feature value**, and the model in
# `adhesion_rf.joblib` was fitted with features computed at 4.0. Training and
# inference must agree, so this is not a free speed knob: lowering it requires
# re-measuring the proba ≥ 0.5 set (and, strictly, retraining). Kept as a named
# constant so a test can vary it without a second copy of the filter code.
TRUNCATE = 4.0

# Chunk sizing. These are the two budgets `_auto_chunk_z` fits the chunk into;
# they are defaults, i.e. decisions (CLAUDE.md), and the chosen `chunk_z` is
# always logged so a run can be reproduced.
#
# `HOST_BUDGET_GB` caps the feature matrix X, which is the binding constraint
# once the clover gate is off: X is (chunk voxels × 60) float32, so a 20-plane
# chunk of a 1024² volume is already 5.0 GB.
# `VRAM_FRACTION` is the share of *free* device memory the slab filters may use.
HOST_BUDGET_GB = 8.0
VRAM_FRACTION = 0.8
# Device memory model for one scale, **measured** rather than counted by eye
# (`docs/scratch/i09/calib.py`, 2026-08-25, 1024² volume, three chunk sizes):
#
#     work_bytes ~= 8.71 x slab_bytes  +  13.3 x (N x 4B)
#
# The first term is padding-bearing (slab = chunk_z + 2*max_pad planes), the
# second is not. Rounded up, because under-estimating costs an OOM while
# over-estimating only costs speed.
#
# The eye-counted value used until 2026-08-25 was 10 slabs with no per-voxel
# term, which **under-estimated peak by ~2x** once the eigenvalues moved to the
# device — 12.03 GB actual against 6.3 GB predicted at chunk_z=24. Batching
# `_eigh3x3_components` cut the per-voxel term from 83.6 floats to 13.3.
_LIVE_SLABS = 9.0
_STAGING_FLOATS = 14
# Never trust a device-memory probe that implies a chunk thinner than this: the
# padding is per-chunk, so chunk_z=1 costs 2*max_pad planes of filtering per
# useful plane. Better to attempt a real allocation and let `gpu.py` fall back.
_MIN_CHUNK_Z = 8


def _auto_chunk_z(shape, max_pad: int, host_budget_gb: float,
                  vram_fraction: float) -> tuple[int, str]:
    """Largest z-chunk fitting both the host and device budgets.

    Bigger is faster: the ±`max_pad` context is a fixed per-chunk overhead, so
    a thin chunk re-filters it many times. The cap is memory, from two
    independent directions, and the tighter one wins.

    Returns (chunk_z, reason) — the reason is logged, because this silently
    sets how a run is executed and must not be invisible.
    """
    z_total, ny, nx = shape
    plane = ny * nx * 4                      # float32 bytes per z-plane
    pad2 = 2 * int(max_pad)

    # Host: X is (chunk_z · ny · nx, N_FEATURES) float32.
    c_host = (host_budget_gb * 1e9) / (ny * nx * N_FEATURES * 4)

    c_vram = float("inf")
    detail = f"host {host_budget_gb:.1f} GB -> {int(c_host)}"
    if _gpu.HAS_GPU:
        try:
            import cupy as cp
            free_b, _total_b = cp.cuda.Device().mem_info
            # The volume is already resident; the slabs come out of what's left.
            usable = free_b * vram_fraction - z_total * plane
            # Per-plane cost has a padding-bearing part and a flat part; the
            # fixed padding is subtracted once. See the model above.
            per_plane = _LIVE_SLABS * plane + _STAGING_FLOATS * ny * nx * 4
            c_vram = (usable - _LIVE_SLABS * pad2 * plane) / per_plane
            detail += (f", vram {free_b/1e9:.1f} GB free -> "
                       f"{int(c_vram) if c_vram > 0 else 0}")
        except Exception as exc:                     # pragma: no cover
            detail += f", vram probe failed ({exc})"

    if c_vram < 1:
        # The floor is the padding, not the chunk: the largest scale needs a
        # (2·max_pad)-plane slab whatever chunk_z is, so shrinking it cannot
        # rescue this. ~6.5 GB is the floor for a 1024² volume at 5 µm/0.325 µm.
        # `gpu.py` falls back to CPU on a genuine OOM, so this warns rather than
        # raises — but the run will be slow and the reason should be visible.
        need = (_LIVE_SLABS * pad2 * plane + z_total * plane) / 1e9
        print(f"  WARNING: the largest scale needs ~{need:.1f} GB of device "
              f"memory before any chunk is added, which does not fit. Expect a "
              f"CPU fallback. Reduce SCALES_UM or use a larger GPU.", flush=True)

    # A transient probe reading must not collapse the chunk. Measured
    # 2026-08-25: a run launched just after another GPU process exited read
    # "0.0 GB free", took chunk_z=1, and spent 41.9 min against 24.3 min for the
    # unoptimised code -- each 1-plane chunk still carries 2*max_pad planes of
    # padding, so tiny chunks are catastrophic, not merely slow. Below the floor
    # the device reading is not believed; the host budget still applies.
    c = int(min(c_host, c_vram))
    if c < _MIN_CHUNK_Z:
        c = int(min(c_host, _MIN_CHUNK_Z))
        detail += (f"; device reading implausible (-> {int(c_vram)}), "
                   f"using {c} and letting it OOM-and-fall-back if real")
    c = max(1, min(int(z_total), c))
    return c, detail


# ---------------------------------------------------------------------------
# Sigma helpers
# ---------------------------------------------------------------------------

def _sigma(scale_um: float, voxel_um) -> list[float]:
    """Per-axis Gaussian sigma in voxel units for a given physical scale."""
    if hasattr(voxel_um, "__len__"):
        return [float(scale_um) / float(v) for v in voxel_um]
    return [float(scale_um) / float(voxel_um)] * 3


# Voxels per eigenvalue sub-batch. The math is elementwise, so batching is
# exact — it only bounds the float64 temporaries, which were measured at ~84
# floats per voxel (8.4 GB at a 25 M-voxel chunk) and otherwise force `chunk_z`
# down to pay for them ([I09](../docs/ISSUES.md#i09)).
_EIGH_BATCH = 4_000_000


def _eigh3x3_components(a00, a11, a22, a01, a02, a12, xp=np):
    """Analytical eigenvalues from the six independent components, ascending.

    Computed in `_EIGH_BATCH`-sized slices: elementwise, so the result is
    identical, but peak memory stops scaling with the caller's chunk size.

    Takes components rather than an (N, 3, 3) array because assembling one costs
    a 0.91 GB allocation and nine strided writes per call, and every caller has
    the components already ([I09](../docs/ISSUES.md#i09)).

    `xp` is numpy or cupy — **one copy of the math for both**, per CLAUDE.md's
    no-second-copy rule. float64 throughout: `q`, `p` and `r` are differences of
    similar magnitudes and cancel badly in float32. On GPU that costs little
    despite consumer cards' 1/64-rate FP64, because this is memory-bound
    (~4 GB of traffic per call) rather than FLOP-bound.
    """
    n = a00.shape[0]
    if n > _EIGH_BATCH:
        out = xp.empty((n, 3), dtype=np.float32)
        for s0 in range(0, n, _EIGH_BATCH):
            s1 = min(s0 + _EIGH_BATCH, n)
            out[s0:s1] = _eigh3x3_components(
                a00[s0:s1], a11[s0:s1], a22[s0:s1],
                a01[s0:s1], a02[s0:s1], a12[s0:s1], xp=xp)
        return out

    a00 = a00.astype(np.float64)
    a11 = a11.astype(np.float64)
    a22 = a22.astype(np.float64)
    a01 = a01.astype(np.float64)
    a02 = a02.astype(np.float64)
    a12 = a12.astype(np.float64)

    q = (a00 + a11 + a22) / 3.0
    p1 = a01 * a01 + a02 * a02 + a12 * a12
    p2 = (a00 - q) ** 2 + (a11 - q) ** 2 + (a22 - q) ** 2 + 2.0 * p1
    p = xp.sqrt(p2 / 6.0)

    p_safe = xp.where(p > 1e-12, p, 1.0)
    inv_p = 1.0 / p_safe
    b00 = (a00 - q) * inv_p;  b11 = (a11 - q) * inv_p;  b22 = (a22 - q) * inv_p
    b01 = a01 * inv_p;  b02 = a02 * inv_p;  b12 = a12 * inv_p

    r = xp.clip(
        0.5 * (b00 * (b11 * b22 - b12 * b12)
               - b01 * (b01 * b22 - b12 * b02)
               + b02 * (b01 * b12 - b11 * b02)),
        -1.0, 1.0,
    )
    phi = xp.arccos(r) / 3.0
    TWO_PI_3 = 2.0 * np.pi / 3.0
    e1 = q + 2.0 * p * xp.cos(phi)
    e2 = q + 2.0 * p * xp.cos(phi + TWO_PI_3)
    e3 = 3.0 * q - e1 - e2
    return xp.sort(xp.stack([e1, e2, e3], axis=1), axis=1).astype(np.float32)


def _eigh3x3(m: np.ndarray) -> np.ndarray:
    """Eigenvalues of (N, 3, 3) real symmetric matrices, ascending.

    Thin wrapper on `_eigh3x3_components` for callers that already hold an
    (N, 3, 3) array. Prefer the component form — see its docstring.
    """
    return _eigh3x3_components(m[:, 0, 0], m[:, 1, 1], m[:, 2, 2],
                               m[:, 0, 1], m[:, 0, 2], m[:, 1, 2], xp=np)


# ---------------------------------------------------------------------------
# GPU feature extraction (one scale at a time)
# ---------------------------------------------------------------------------

def _extract_scale_gpu(a, iz, iy, ix, sigma: list[float],
                        scale_um: float, F: np.ndarray, col: int,
                        truncate: float = TRUNCATE) -> None:
    """Compute 10 features for one sigma on GPU; write into F[:, col:col+10].

    `a` is a cupy array already on device — caller owns upload/free.

    `truncate` is where the Gaussian kernel is cut off, in units of sigma. It
    sets the cost (taps per axis = 2·ceil(truncate·sigma)+1) and it changes the
    feature values, so it must match between training and inference — see
    TRUNCATE.
    """
    def _gf(x, sig, order=None):
        return (_cndi.gaussian_filter(x, sig, truncate=truncate) if order is None
                else _cndi.gaussian_filter(x, sig, order=order, truncate=truncate))
    s2 = scale_um ** 2
    sqrt2 = float(np.sqrt(2))
    N = len(iz)

    # All ten features are built into one device-side (N, 10) block and moved
    # host-side in a single transfer, then written to F as one contiguous slab.
    # Feature-at-a-time was measured at 140 ms per column ([I09](
    # ../docs/ISSUES.md#i09)): F is C-order (N, 60), so `F[:, col] = ...` writes
    # 4 bytes per 240-byte row and touches a cache line for each one.
    B = cp.empty((N, N_PER_SCALE), dtype=cp.float32)

    # 1. Gaussian smoothed
    g = _gf(a, sigma)
    B[:, 0] = g[iz, iy, ix]

    # 2. Gradient magnitude (1st-order derivatives)
    gz = _gf(a, sigma, order=[1, 0, 0])
    gy = _gf(a, sigma, order=[0, 1, 0])
    gx = _gf(a, sigma, order=[0, 0, 1])
    czg = gz[iz, iy, ix];  cyg = gy[iz, iy, ix];  cxg = gx[iz, iy, ix]
    B[:, 1] = cp.sqrt(czg * czg + cyg * cyg + cxg * cxg)
    del czg, cyg, cxg

    # 3. LoG (scale-normalised Laplacian of Gaussian)
    # Gathered once and reused by the Hessian below — previously these were
    # gathered twice, and the sum built a full slab-sized temporary (629 MB at
    # 5 µm) that was only ever read at the gathered indices.
    gzz = _gf(a, sigma, order=[2, 0, 0]);  czz = gzz[iz, iy, ix] * s2;  del gzz
    gyy = _gf(a, sigma, order=[0, 2, 0]);  cyy = gyy[iz, iy, ix] * s2;  del gyy
    gxx = _gf(a, sigma, order=[0, 0, 2]);  cxx = gxx[iz, iy, ix] * s2;  del gxx
    B[:, 2] = czz + cyy + cxx

    # 4. DoG  (g − gaussian(σ·√2))
    sigma_dog = [s * sqrt2 for s in sigma]
    g2 = _gf(a, sigma_dog)
    B[:, 3] = (g - g2)[iz, iy, ix]
    del g, g2

    # 5. Hessian eigenvalues (scale-normalised) — on device, from components.
    gzy = _gf(a, sigma, order=[1, 1, 0]);  czy = gzy[iz, iy, ix] * s2;  del gzy
    gzx = _gf(a, sigma, order=[1, 0, 1]);  czx = gzx[iz, iy, ix] * s2;  del gzx
    gyx = _gf(a, sigma, order=[0, 1, 1]);  cyx = gyx[iz, iy, ix] * s2;  del gyx
    B[:, 4:7] = _eigh3x3_components(czz, cyy, cxx, czy, czx, cyx, xp=cp)
    del czz, cyy, cxx, czy, czx, cyx

    # 6. Structure tensor eigenvalues
    gz2 = gz * gz;  gy2 = gy * gy;  gx2 = gx * gx
    gzgy = gz * gy; gzgx = gz * gx; gygx = gy * gx
    del gz, gy, gx

    _tmp = _gf(gz2, sigma);  del gz2;  s00 = _tmp[iz, iy, ix];  del _tmp
    _tmp = _gf(gy2, sigma);  del gy2;  s11 = _tmp[iz, iy, ix];  del _tmp
    _tmp = _gf(gx2, sigma);  del gx2;  s22 = _tmp[iz, iy, ix];  del _tmp
    _tmp = _gf(gzgy, sigma); del gzgy; s01 = _tmp[iz, iy, ix];  del _tmp
    _tmp = _gf(gzgx, sigma); del gzgx; s02 = _tmp[iz, iy, ix];  del _tmp
    _tmp = _gf(gygx, sigma); del gygx; s12 = _tmp[iz, iy, ix];  del _tmp
    B[:, 7:10] = _eigh3x3_components(s00, s11, s22, s01, s02, s12, xp=cp)
    del s00, s11, s22, s01, s02, s12

    F[:, col:col + N_PER_SCALE] = B.get()
    del B


# ---------------------------------------------------------------------------
# CPU feature extraction (one scale at a time, sequential)
# ---------------------------------------------------------------------------

def _extract_scale_cpu(vol: np.ndarray, iz, iy, ix, sigma: list[float],
                        scale_um: float, F: np.ndarray, col: int,
                        truncate: float = TRUNCATE) -> None:
    """Same feature set on CPU via scipy.ndimage. See `_extract_scale_gpu`."""
    def _gf(x, sig, order=None):
        return (_sndi.gaussian_filter(x, sig, truncate=truncate) if order is None
                else _sndi.gaussian_filter(x, sig, order=order, truncate=truncate))
    s2 = scale_um ** 2
    sqrt2 = float(np.sqrt(2))
    N = len(iz)

    # Mirrors the GPU path exactly: one (N, 10) block, contiguous write into F,
    # eigenvalues from components with no (N, 3, 3) assembly.
    B = np.empty((N, N_PER_SCALE), dtype=np.float32)

    g = _gf(vol, sigma)
    B[:, 0] = g[iz, iy, ix]

    gz = _gf(vol, sigma, order=[1, 0, 0])
    gy = _gf(vol, sigma, order=[0, 1, 0])
    gx = _gf(vol, sigma, order=[0, 0, 1])
    czg = gz[iz, iy, ix];  cyg = gy[iz, iy, ix];  cxg = gx[iz, iy, ix]
    B[:, 1] = np.sqrt(czg * czg + cyg * cyg + cxg * cxg)
    del czg, cyg, cxg

    gzz = _gf(vol, sigma, order=[2, 0, 0]);  czz = gzz[iz, iy, ix] * s2;  del gzz
    gyy = _gf(vol, sigma, order=[0, 2, 0]);  cyy = gyy[iz, iy, ix] * s2;  del gyy
    gxx = _gf(vol, sigma, order=[0, 0, 2]);  cxx = gxx[iz, iy, ix] * s2;  del gxx
    B[:, 2] = czz + cyy + cxx

    sigma_dog = [s * sqrt2 for s in sigma]
    g2 = _gf(vol, sigma_dog)
    B[:, 3] = (g - g2)[iz, iy, ix]
    del g, g2

    gzy = _gf(vol, sigma, order=[1, 1, 0]);  czy = gzy[iz, iy, ix] * s2;  del gzy
    gzx = _gf(vol, sigma, order=[1, 0, 1]);  czx = gzx[iz, iy, ix] * s2;  del gzx
    gyx = _gf(vol, sigma, order=[0, 1, 1]);  cyx = gyx[iz, iy, ix] * s2;  del gyx
    B[:, 4:7] = _eigh3x3_components(czz, cyy, cxx, czy, czx, cyx, xp=np)
    del czz, cyy, cxx, czy, czx, cyx

    gz2 = gz * gz;  gy2 = gy * gy;  gx2 = gx * gx
    gzgy = gz * gy; gzgx = gz * gx; gygx = gy * gx
    del gz, gy, gx

    s00 = _gf(gz2, sigma)[iz, iy, ix];   del gz2
    s11 = _gf(gy2, sigma)[iz, iy, ix];   del gy2
    s22 = _gf(gx2, sigma)[iz, iy, ix];   del gx2
    s01 = _gf(gzgy, sigma)[iz, iy, ix];  del gzgy
    s02 = _gf(gzgx, sigma)[iz, iy, ix];  del gzgx
    s12 = _gf(gygx, sigma)[iz, iy, ix];  del gygx
    B[:, 7:10] = _eigh3x3_components(s00, s11, s22, s01, s02, s12, xp=np)
    del s00, s11, s22, s01, s02, s12

    F[:, col:col + N_PER_SCALE] = B


# ---------------------------------------------------------------------------
# Public feature-extraction entry point
# ---------------------------------------------------------------------------

def compute_features_at_indices(
    vol_f32: np.ndarray,
    indices: np.ndarray,
    voxel_um,
    _time: bool = False,
    _vol_gpu=None,        # pre-uploaded cupy array; if given, skips H2D transfer
) -> np.ndarray:
    """Compute (N, 70) float32 feature matrix at the given voxel positions.

    Parameters
    ----------
    vol_f32 : (Z, Y, X) float32 numpy volume (used for CPU path or if _vol_gpu is None)
    indices : (N, 3) int array of (z, y, x) positions
    voxel_um : scalar or (vz, vy, vx) tuple — physical voxel size in µm
    _time : if True, print per-scale timing (filter + extract breakdown)
    _vol_gpu : optional pre-uploaded cupy array; avoids repeated H2D transfers
    """
    import time as _time_mod
    indices = np.asarray(indices, dtype=np.int32)
    N = len(indices)
    F = np.empty((N, N_FEATURES), dtype=np.float32)

    iz = np.ascontiguousarray(indices[:, 0])
    iy = np.ascontiguousarray(indices[:, 1])
    ix = np.ascontiguousarray(indices[:, 2])

    if _gpu.HAS_GPU:
        t_up = _time_mod.perf_counter()
        owns_gpu = _vol_gpu is None
        a_d = cp.asarray(vol_f32) if owns_gpu else _vol_gpu
        iz_d = cp.asarray(iz)
        iy_d = cp.asarray(iy)
        ix_d = cp.asarray(ix)
        if _time:
            cp.cuda.stream.get_current_stream().synchronize()
            mb = (vol_f32.nbytes if owns_gpu else 0) / 1e6
            print(f"    upload: {_time_mod.perf_counter()-t_up:.2f}s  "
                  f"({mb:.0f} MB)", flush=True)

        for si, scale_um in enumerate(SCALES_UM):
            t0 = _time_mod.perf_counter()
            sig = _sigma(scale_um, voxel_um)
            _extract_scale_gpu(a_d, iz_d, iy_d, ix_d, sig, scale_um,
                                F, si * N_PER_SCALE)
            if _time:
                cp.cuda.stream.get_current_stream().synchronize()
                print(f"    scale {scale_um:5.1f} µm: "
                      f"{_time_mod.perf_counter()-t0:.2f}s",
                      flush=True)

        if owns_gpu:
            del a_d
        del iz_d, iy_d, ix_d
    else:
        for si, scale_um in enumerate(SCALES_UM):
            t0 = _time_mod.perf_counter()
            sig = _sigma(scale_um, voxel_um)
            _extract_scale_cpu(vol_f32, iz, iy, ix, sig, scale_um,
                                F, si * N_PER_SCALE)
            if _time:
                print(f"    scale {scale_um:5.1f} µm: "
                      f"{_time_mod.perf_counter()-t0:.2f}s",
                      flush=True)

    return F


# ---------------------------------------------------------------------------
# AdhesionRF class
# ---------------------------------------------------------------------------

class AdhesionRF:
    def __init__(self, n_estimators: int = 200):
        from sklearn.ensemble import RandomForestClassifier
        self.rf = RandomForestClassifier(
            n_estimators=n_estimators,
            n_jobs=-1,
            random_state=42,
            class_weight="balanced",
            min_samples_leaf=5,
            max_features="sqrt",
        )
        self.is_trained = False

    def train(self, vol_f32: np.ndarray, ann_vol: np.ndarray, voxel_um) -> None:
        """Fit the RF on annotated voxels.

        ann_vol: (Z, Y, X) uint8 — 1=adhesion, 2=background, 0=unlabeled
        """
        idx = np.argwhere(ann_vol > 0).astype(np.int32)
        if len(idx) == 0:
            raise ValueError("No annotations found (ann_vol is all-zero).")

        y = ann_vol[idx[:, 0], idx[:, 1], idx[:, 2]]

        unique, counts = np.unique(y, return_counts=True)
        for cls, cnt in zip(unique, counts):
            name = {1: "adhesion", 2: "background"}.get(int(cls), f"class {cls}")
            print(f"  {name}: {cnt} labeled voxels", flush=True)

        if not (1 in unique and 2 in unique):
            raise ValueError(
                "Need both adhesion (1) and background (2) annotations to train."
            )

        print(f"  Extracting features for {len(idx)} labeled voxels ...", flush=True)
        X = compute_features_at_indices(vol_f32, idx, voxel_um)

        print("  Fitting RF ...", flush=True)
        self.rf.fit(X, y)
        self.is_trained = True

    def predict_proba(
        self,
        vol_f32: np.ndarray,
        voxel_um,
        gate_indices,          # (N, 3) int32, or None = full volume
        chunk_z: int | None = None,
        host_budget_gb: float = HOST_BUDGET_GB,
        vram_fraction: float = VRAM_FRACTION,
        truncate: float = TRUNCATE,
    ) -> np.ndarray:
        """Return (Z, Y, X) float32 adhesion probability volume.

        Processes in z-chunks so the feature matrix stays within
        `host_budget_gb` regardless of volume size. gate_indices=None means
        predict everywhere; indices are generated per chunk to avoid a 4 GB
        index array.

        Each scale uses its own minimum z-padding (4σ) rather than the global
        worst-case pad, so small-scale filters process far fewer voxels.

        **`chunk_z=None` sizes the chunk automatically** — see `_auto_chunk_z`.
        A bigger chunk amortises the fixed ±4σ padding over more useful output
        (at the 5 µm scale a 20-plane chunk filters 146 planes to keep 20), so
        the largest chunk that fits is also the fastest. Pass an int to pin it.
        """
        import time as _tm

        shape = vol_f32.shape
        proba_vol = np.zeros(shape, dtype=np.uint8)

        vox_min = min(voxel_um) if hasattr(voxel_um, "__len__") else float(voxel_um)
        # Minimum padding per scale instead of one global worst-case value.
        scale_pads = [int(np.ceil(4.0 * s / vox_min)) + 1 for s in SCALES_UM]
        max_pad = max(scale_pads)

        if chunk_z is None:
            chunk_z, why = _auto_chunk_z(shape, max_pad, host_budget_gb,
                                         vram_fraction)
            print(f"  chunk_z={chunk_z} (auto: {why})", flush=True)
        else:
            chunk_z = int(chunk_z)

        adhesion_col = list(self.rf.classes_).index(1)

        full_volume = gate_indices is None
        z_total = shape[0]
        n_chunks = int(np.ceil(z_total / chunk_z))

        # Upload full volume to GPU once; each chunk takes a zero-copy slice.
        vol_gpu = None
        if _gpu.HAS_GPU:
            import cupy as cp
            t0 = _tm.perf_counter()
            vol_gpu = cp.asarray(vol_f32)
            cp.cuda.stream.get_current_stream().synchronize()
            print(f"  GPU upload: {vol_f32.nbytes/1e6:.0f} MB in "
                  f"{_tm.perf_counter()-t0:.1f}s", flush=True)

        print(f"  Predicting in {n_chunks} z-chunks "
              f"({'full volume' if full_volume else 'gated'}, max_pad={max_pad}) ...",
              flush=True)

        BATCH = 2_000_000
        ci = 0
        _Xbuf = None
        if not full_volume:
            # Safety net: never index the output volume out of bounds if the
            # gate indices don't match this volume's shape.
            Zs, Ys, Xs = shape
            ok = ((gate_indices[:, 0] >= 0) & (gate_indices[:, 0] < Zs) &
                  (gate_indices[:, 1] >= 0) & (gate_indices[:, 1] < Ys) &
                  (gate_indices[:, 2] >= 0) & (gate_indices[:, 2] < Xs))
            if not ok.all():
                print(f"  dropping {int((~ok).sum()):,} out-of-bounds gate voxels "
                      f"(volume {shape}).", flush=True)
                gate_indices = gate_indices[ok]
        gate_iz = None if full_volume else gate_indices[:, 0]

        for z_start in range(0, z_total, chunk_z):
            z_end = min(z_start + chunk_z, z_total)
            t_chunk = _tm.perf_counter()

            if full_volume:
                nz = z_end - z_start
                ny, nx = shape[1], shape[2]
                N = nz * ny * nx
                chunk_iz = np.repeat(
                    np.arange(z_start, z_end, dtype=np.int32), ny * nx)
                chunk_iy = np.tile(
                    np.repeat(np.arange(ny, dtype=np.int32), nx), nz)
                chunk_ix = np.tile(np.arange(nx, dtype=np.int32), nz * ny)
                chunk_idx = np.stack([chunk_iz, chunk_iy, chunk_ix], axis=1)
            else:
                in_chunk = (gate_iz >= z_start) & (gate_iz < z_end)
                chunk_idx = gate_indices[in_chunk]
                if len(chunk_idx) == 0:
                    continue
                N = len(chunk_idx)

            _do_time = (ci == 0)
            t0 = _tm.perf_counter()
            # Allocated once and reused: `np.empty` per chunk means first-touch
            # page faults on 6+ GB every chunk, for a buffer whose size barely
            # changes ([I09](../docs/ISSUES.md#i09)).
            if _Xbuf is None or len(_Xbuf) < N:
                _Xbuf = np.empty((N, N_FEATURES), dtype=np.float32)
            X = _Xbuf[:N]

            if vol_gpu is not None:
                # Keep y/x index arrays on GPU; z-offset differs per scale.
                iz_base = cp.asarray(chunk_idx[:, 0])
                iy_d = cp.asarray(chunk_idx[:, 1])
                ix_d = cp.asarray(chunk_idx[:, 2])

                for si, (scale_um, pad_s) in enumerate(zip(SCALES_UM, scale_pads)):
                    t_s = _tm.perf_counter()
                    z_lo_s = max(0, z_start - pad_s)
                    z_hi_s = min(shape[0], z_end + pad_s)
                    iz_d = iz_base - z_lo_s
                    sig = _sigma(scale_um, voxel_um)
                    _extract_scale_gpu(vol_gpu[z_lo_s:z_hi_s], iz_d, iy_d, ix_d,
                                       sig, scale_um, X, si * N_PER_SCALE,
                                       truncate=truncate)
                    if _do_time:
                        cp.cuda.stream.get_current_stream().synchronize()
                        print(f"    scale {scale_um:5.1f} µm: "
                              f"{_tm.perf_counter()-t_s:.2f}s", flush=True)

                del iz_base, iz_d, iy_d, ix_d
            else:
                for si, (scale_um, pad_s) in enumerate(zip(SCALES_UM, scale_pads)):
                    t_s = _tm.perf_counter()
                    z_lo_s = max(0, z_start - pad_s)
                    z_hi_s = min(shape[0], z_end + pad_s)
                    idx_local = chunk_idx.copy()
                    idx_local[:, 0] -= z_lo_s
                    sig = _sigma(scale_um, voxel_um)
                    _extract_scale_cpu(vol_f32[z_lo_s:z_hi_s],
                                       idx_local[:, 0], idx_local[:, 1],
                                       idx_local[:, 2], sig, scale_um, X,
                                       si * N_PER_SCALE, truncate=truncate)
                    if _do_time:
                        print(f"    scale {scale_um:5.1f} µm: "
                              f"{_tm.perf_counter()-t_s:.2f}s", flush=True)

            t_feat = _tm.perf_counter() - t0

            t0 = _tm.perf_counter()
            proba_chunk = np.empty(N, dtype=np.float32)
            for start in range(0, N, BATCH):
                end = min(start + BATCH, N)
                proba_chunk[start:end] = self.rf.predict_proba(
                    X[start:end])[:, adhesion_col]
            del X          # drops the view, not `_Xbuf`
            t_rf = _tm.perf_counter() - t0

            proba_vol[chunk_idx[:, 0], chunk_idx[:, 1],
                      chunk_idx[:, 2]] = (proba_chunk * 255.0 + 0.5).astype(np.uint8)
            ci += 1
            t_total = _tm.perf_counter() - t_chunk
            print(f"  chunk {ci}/{n_chunks}  z={z_start}–{z_end-1}  "
                  f"{N} vox  feat={t_feat:.1f}s  rf={t_rf:.1f}s  "
                  f"total={t_total:.1f}s", flush=True)

        if vol_gpu is not None:
            del vol_gpu

        return proba_vol

    def save(self, path: str) -> None:
        import joblib
        joblib.dump(self, path)
        print(f"RF model saved to {path}", flush=True)

    @classmethod
    def load(cls, path: str) -> "AdhesionRF":
        import joblib
        obj = joblib.load(path)
        print(f"RF model loaded from {path} (trained={obj.is_trained})", flush=True)
        return obj


# ---------------------------------------------------------------------------
# Instance labeling
# ---------------------------------------------------------------------------

def instance_labels_3d(
    binary: np.ndarray,
    min_size_vox: int = 100,
    split_merged: bool = False,
    split_min_size: int = 1000,
) -> np.ndarray:
    """Connected components → optional watershed splitting → (Z, Y, X) int32.

    Prints progress through the splitting loop: it is the slow phase and it is
    otherwise completely silent, so a long run is indistinguishable from a hang.
    """
    from scipy.ndimage import (
        label as _label,
        distance_transform_edt as _edt,
        find_objects as _find_objects,
    )

    lbl, _ = _label(binary)

    # Remove small components
    if min_size_vox > 0:
        sizes = np.bincount(lbl.ravel())
        keep = (sizes >= min_size_vox)
        keep[0] = True  # background
        remap = np.where(keep, np.arange(len(keep), dtype=np.int32), 0)
        lbl = remap[lbl]

    if split_merged and split_min_size > 0 and lbl.max() > 0:
        from skimage.feature import peak_local_max
        from skimage.segmentation import watershed

        new_max = int(lbl.max())
        comp_ids, comp_sizes = np.unique(lbl[lbl > 0], return_counts=True)
        large = comp_ids[comp_sizes > split_min_size]

        # Every component's bounding box in ONE pass. The obvious
        # `np.where(lbl == cid)` inside the loop costs a full-volume scan per
        # component (~180 ms on a 65M-voxel stack), so a few thousand candidates
        # burn minutes just locating crops. Computing this up front is safe: the
        # loop only ever relabels voxels belonging to the component it is
        # currently splitting, so the other components' boxes stay valid.
        boxes = _find_objects(lbl)

        n_large = int(len(large))
        print(f"Splitting {n_large} merged component(s) (> {split_min_size} "
              f"vox) of {len(comp_ids)} total ...", flush=True)
        every = max(1, n_large // 20)
        n_split = 0

        for k, cid in enumerate(large, 1):
            if k % every == 0 or k == n_large:
                print(f"  {k}/{n_large} examined, {n_split} divided",
                      flush=True)
            # Crop to bounding box so EDT/watershed run on the component only.
            sl = boxes[cid - 1]
            if sl is None:                      # label absent (shouldn't happen)
                continue
            # Pad wider than peak_local_max's min_distance so an object's own
            # crop edge is never within the exclusion radius of its peaks.
            pad = 6
            z0, z1 = max(0, sl[0].start - pad), min(lbl.shape[0], sl[0].stop + pad)
            y0, y1 = max(0, sl[1].start - pad), min(lbl.shape[1], sl[1].stop + pad)
            x0, x1 = max(0, sl[2].start - pad), min(lbl.shape[2], sl[2].stop + pad)

            crop_mask = (lbl[z0:z1, y0:y1, x0:x1] == cid)
            dist = _edt(crop_mask)
            # exclude_border=False as well as the pad: the pad is clipped at the
            # volume bounds, so it cannot rescue objects touching the stack
            # edge. Without both, those objects keep at most one seed and never
            # split.
            coords = peak_local_max(dist, min_distance=5, labels=crop_mask,
                                    exclude_border=False)
            if len(coords) <= 1:
                continue
            n_split += 1

            markers = np.zeros(crop_mask.shape, dtype=np.int32)
            markers[coords[:, 0], coords[:, 1], coords[:, 2]] = np.arange(
                1, len(coords) + 1, dtype=np.int32
            )
            sub = watershed(-dist, markers, mask=crop_mask)

            lbl[z0:z1, y0:y1, x0:x1][crop_mask] = 0
            for sl_id in np.unique(sub[sub > 0]):
                new_max += 1
                lbl[z0:z1, y0:y1, x0:x1][sub == sl_id] = new_max

    # Watershed splitting can carve fragments smaller than min_size_vox. The
    # size filter above ran BEFORE the split, so re-apply it here to drop those
    # post-split slivers too — otherwise sub-threshold pieces leak through and
    # min_size_vox silently stops meaning what it says once split_merged is on.
    if split_merged and min_size_vox > 0 and lbl.max() > 0:
        sizes = np.bincount(lbl.ravel())
        keep = (sizes >= min_size_vox)
        keep[0] = True  # background
        remap = np.where(keep, np.arange(len(keep), dtype=np.int32), 0)
        lbl = remap[lbl]

    # Relabel sequentially
    from skimage.segmentation import relabel_sequential
    lbl, _, _ = relabel_sequential(lbl)
    return lbl.astype(np.int32)
