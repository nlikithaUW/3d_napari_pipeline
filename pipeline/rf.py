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


# ---------------------------------------------------------------------------
# Sigma helpers
# ---------------------------------------------------------------------------

def _sigma(scale_um: float, voxel_um) -> list[float]:
    """Per-axis Gaussian sigma in voxel units for a given physical scale."""
    if hasattr(voxel_um, "__len__"):
        return [float(scale_um) / float(v) for v in voxel_um]
    return [float(scale_um) / float(voxel_um)] * 3


def _eigh3x3(m: np.ndarray) -> np.ndarray:
    """Analytical eigenvalues of (N, 3, 3) real symmetric matrices, sorted ascending.

    Vectorised over the batch dimension — no per-matrix LAPACK dispatch.
    Typically 10–50× faster than np.linalg.eigh for large N and small (3×3).
    """
    a00 = m[:, 0, 0].astype(np.float64)
    a11 = m[:, 1, 1].astype(np.float64)
    a22 = m[:, 2, 2].astype(np.float64)
    a01 = m[:, 0, 1].astype(np.float64)
    a02 = m[:, 0, 2].astype(np.float64)
    a12 = m[:, 1, 2].astype(np.float64)

    q = (a00 + a11 + a22) / 3.0
    p1 = a01 * a01 + a02 * a02 + a12 * a12
    p2 = (a00 - q) ** 2 + (a11 - q) ** 2 + (a22 - q) ** 2 + 2.0 * p1
    p = np.sqrt(p2 / 6.0)

    p_safe = np.where(p > 1e-12, p, 1.0)
    inv_p = 1.0 / p_safe
    b00 = (a00 - q) * inv_p;  b11 = (a11 - q) * inv_p;  b22 = (a22 - q) * inv_p
    b01 = a01 * inv_p;  b02 = a02 * inv_p;  b12 = a12 * inv_p

    r = np.clip(
        0.5 * (b00 * (b11 * b22 - b12 * b12)
               - b01 * (b01 * b22 - b12 * b02)
               + b02 * (b01 * b12 - b11 * b02)),
        -1.0, 1.0,
    )
    phi = np.arccos(r) / 3.0
    TWO_PI_3 = 2.0 * np.pi / 3.0
    e1 = q + 2.0 * p * np.cos(phi)
    e2 = q + 2.0 * p * np.cos(phi + TWO_PI_3)
    e3 = 3.0 * q - e1 - e2
    return np.sort(np.stack([e1, e2, e3], axis=1), axis=1).astype(np.float32)


# ---------------------------------------------------------------------------
# GPU feature extraction (one scale at a time)
# ---------------------------------------------------------------------------

def _extract_scale_gpu(a, iz, iy, ix, sigma: list[float],
                        scale_um: float, F: np.ndarray, col: int) -> None:
    """Compute 10 features for one sigma on GPU; write into F[:, col:col+10].

    `a` is a cupy array already on device — caller owns upload/free.
    """
    s2 = scale_um ** 2
    sqrt2 = float(np.sqrt(2))

    # 1. Gaussian smoothed
    g = _cndi.gaussian_filter(a, sigma)
    F[:, col] = g[iz, iy, ix].get()

    # 2. Gradient magnitude (1st-order derivatives)
    gz = _cndi.gaussian_filter(a, sigma, order=[1, 0, 0])
    gy = _cndi.gaussian_filter(a, sigma, order=[0, 1, 0])
    gx = _cndi.gaussian_filter(a, sigma, order=[0, 0, 1])
    gm = cp.sqrt(gz[iz, iy, ix] ** 2 + gy[iz, iy, ix] ** 2 + gx[iz, iy, ix] ** 2)
    F[:, col + 1] = gm.get()
    del gm

    # 3. LoG (scale-normalised Laplacian of Gaussian)
    gzz = _cndi.gaussian_filter(a, sigma, order=[2, 0, 0])
    gyy = _cndi.gaussian_filter(a, sigma, order=[0, 2, 0])
    gxx = _cndi.gaussian_filter(a, sigma, order=[0, 0, 2])
    F[:, col + 2] = (s2 * (gzz + gyy + gxx)[iz, iy, ix]).get()

    # 4. DoG  (g − gaussian(σ·√2))
    sigma_dog = [s * sqrt2 for s in sigma]
    g2 = _cndi.gaussian_filter(a, sigma_dog)
    F[:, col + 3] = (g - g2)[iz, iy, ix].get()
    del g, g2

    # 5. Hessian eigenvalues (scale-normalised)
    gzy = _cndi.gaussian_filter(a, sigma, order=[1, 1, 0])
    gzx = _cndi.gaussian_filter(a, sigma, order=[1, 0, 1])
    gyx = _cndi.gaussian_filter(a, sigma, order=[0, 1, 1])

    N = len(iz)
    H = np.empty((N, 3, 3), dtype=np.float32)
    H[:, 0, 0] = gzz[iz, iy, ix].get() * s2
    H[:, 1, 1] = gyy[iz, iy, ix].get() * s2
    H[:, 2, 2] = gxx[iz, iy, ix].get() * s2
    _v = gzy[iz, iy, ix].get() * s2;  H[:, 0, 1] = _v;  H[:, 1, 0] = _v
    _v = gzx[iz, iy, ix].get() * s2;  H[:, 0, 2] = _v;  H[:, 2, 0] = _v
    _v = gyx[iz, iy, ix].get() * s2;  H[:, 1, 2] = _v;  H[:, 2, 1] = _v
    del gzz, gyy, gxx, gzy, gzx, gyx
    F[:, col + 4:col + 7] = _eigh3x3(H)

    # 6. Structure tensor eigenvalues
    gz2 = gz * gz;  gy2 = gy * gy;  gx2 = gx * gx
    gzgy = gz * gy; gzgx = gz * gx; gygx = gy * gx
    del gz, gy, gx

    ST = np.empty((N, 3, 3), dtype=np.float32)
    _tmp = _cndi.gaussian_filter(gz2, sigma);  del gz2
    ST[:, 0, 0] = _tmp[iz, iy, ix].get();  del _tmp
    _tmp = _cndi.gaussian_filter(gy2, sigma);  del gy2
    ST[:, 1, 1] = _tmp[iz, iy, ix].get();  del _tmp
    _tmp = _cndi.gaussian_filter(gx2, sigma);  del gx2
    ST[:, 2, 2] = _tmp[iz, iy, ix].get();  del _tmp
    _tmp = _cndi.gaussian_filter(gzgy, sigma); del gzgy
    _v = _tmp[iz, iy, ix].get();  del _tmp
    ST[:, 0, 1] = _v;  ST[:, 1, 0] = _v
    _tmp = _cndi.gaussian_filter(gzgx, sigma); del gzgx
    _v = _tmp[iz, iy, ix].get();  del _tmp
    ST[:, 0, 2] = _v;  ST[:, 2, 0] = _v
    _tmp = _cndi.gaussian_filter(gygx, sigma); del gygx
    _v = _tmp[iz, iy, ix].get();  del _tmp
    ST[:, 1, 2] = _v;  ST[:, 2, 1] = _v

    F[:, col + 7:col + 10] = _eigh3x3(ST)


# ---------------------------------------------------------------------------
# CPU feature extraction (one scale at a time, sequential)
# ---------------------------------------------------------------------------

def _extract_scale_cpu(vol: np.ndarray, iz, iy, ix, sigma: list[float],
                        scale_um: float, F: np.ndarray, col: int) -> None:
    """Same feature set on CPU via scipy.ndimage."""
    s2 = scale_um ** 2
    sqrt2 = float(np.sqrt(2))
    N = len(iz)

    g = _sndi.gaussian_filter(vol, sigma)
    F[:, col] = g[iz, iy, ix]

    gz = _sndi.gaussian_filter(vol, sigma, order=[1, 0, 0])
    gy = _sndi.gaussian_filter(vol, sigma, order=[0, 1, 0])
    gx = _sndi.gaussian_filter(vol, sigma, order=[0, 0, 1])
    F[:, col + 1] = np.sqrt(gz[iz, iy, ix] ** 2 + gy[iz, iy, ix] ** 2 + gx[iz, iy, ix] ** 2)

    gzz = _sndi.gaussian_filter(vol, sigma, order=[2, 0, 0])
    gyy = _sndi.gaussian_filter(vol, sigma, order=[0, 2, 0])
    gxx = _sndi.gaussian_filter(vol, sigma, order=[0, 0, 2])
    F[:, col + 2] = s2 * (gzz + gyy + gxx)[iz, iy, ix]

    sigma_dog = [s * sqrt2 for s in sigma]
    g2 = _sndi.gaussian_filter(vol, sigma_dog)
    F[:, col + 3] = (g - g2)[iz, iy, ix]
    del g, g2

    gzy = _sndi.gaussian_filter(vol, sigma, order=[1, 1, 0])
    gzx = _sndi.gaussian_filter(vol, sigma, order=[1, 0, 1])
    gyx = _sndi.gaussian_filter(vol, sigma, order=[0, 1, 1])

    H = np.empty((N, 3, 3), dtype=np.float32)
    H[:, 0, 0] = gzz[iz, iy, ix] * s2
    H[:, 1, 1] = gyy[iz, iy, ix] * s2
    H[:, 2, 2] = gxx[iz, iy, ix] * s2
    _v = gzy[iz, iy, ix] * s2;  H[:, 0, 1] = _v;  H[:, 1, 0] = _v
    _v = gzx[iz, iy, ix] * s2;  H[:, 0, 2] = _v;  H[:, 2, 0] = _v
    _v = gyx[iz, iy, ix] * s2;  H[:, 1, 2] = _v;  H[:, 2, 1] = _v
    del gzz, gyy, gxx, gzy, gzx, gyx
    F[:, col + 4:col + 7] = _eigh3x3(H)

    gz2 = gz * gz;  gy2 = gy * gy;  gx2 = gx * gx
    gzgy = gz * gy; gzgx = gz * gx; gygx = gy * gx
    del gz, gy, gx

    ST = np.empty((N, 3, 3), dtype=np.float32)
    ST[:, 0, 0] = _sndi.gaussian_filter(gz2, sigma)[iz, iy, ix];  del gz2
    ST[:, 1, 1] = _sndi.gaussian_filter(gy2, sigma)[iz, iy, ix];  del gy2
    ST[:, 2, 2] = _sndi.gaussian_filter(gx2, sigma)[iz, iy, ix];  del gx2
    _v = _sndi.gaussian_filter(gzgy, sigma)[iz, iy, ix]; del gzgy
    ST[:, 0, 1] = _v;  ST[:, 1, 0] = _v
    _v = _sndi.gaussian_filter(gzgx, sigma)[iz, iy, ix]; del gzgx
    ST[:, 0, 2] = _v;  ST[:, 2, 0] = _v
    _v = _sndi.gaussian_filter(gygx, sigma)[iz, iy, ix]; del gygx
    ST[:, 1, 2] = _v;  ST[:, 2, 1] = _v

    F[:, col + 7:col + 10] = _eigh3x3(ST)


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
# Gate indices
# ---------------------------------------------------------------------------

def _gate_indices(
    shape: tuple,
    dapi_labels,
    mesh_cache: dict,
    d_dapi_um: float,
    gate_dapi: bool,
    gate_mesh: bool,
    d_mesh_um: float,
    voxel_um=(1.0, 1.0, 1.0),
) -> np.ndarray:
    """Return (N, 3) int32 array of voxel positions to predict at.

    If no gate is active, returns all voxel positions (full volume).
    """
    if not gate_dapi and not gate_mesh:
        return None  # sentinel: predict everywhere; indices generated per chunk

    mask = np.ones(shape, dtype=bool)

    if gate_dapi and dapi_labels is not None:
        sampling = tuple(float(v) for v in voxel_um) if hasattr(voxel_um, "__len__") \
                   else (float(voxel_um),) * 3
        dist_um = _sndi.distance_transform_edt(dapi_labels == 0, sampling=sampling)
        mask &= dist_um <= float(d_dapi_um)
    elif gate_dapi:
        print("DAPI gate requested but no DAPI labels available — skipping.", flush=True)

    if gate_mesh:
        mc = mesh_cache or {}
        if mc.get("offset_centroids") is not None and len(mc["offset_centroids"]):
            from .cell_mesh import signed_distance_to_mesh_vertices
            scale = np.asarray(voxel_um, dtype=np.float64) if hasattr(voxel_um, "__len__") \
                    else np.full(3, float(voxel_um))
            idx = np.argwhere(mask)
            if len(idx):
                query_um = idx.astype(np.float64) * scale[None, :]
                sd = signed_distance_to_mesh_vertices(
                    query_um, mc["offset_centroids"], mc["outward"],
                    tree=mc.get("kdtree"),
                )
                keep = (sd > -float(d_mesh_um)) & (sd <= float(d_mesh_um))
                drop = idx[~keep]
                if len(drop):
                    mask[drop[:, 0], drop[:, 1], drop[:, 2]] = False
        else:
            print("Mesh gate requested but no cell mesh built — skipping.", flush=True)

    return np.argwhere(mask).astype(np.int32)


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
        chunk_z: int = 20,
    ) -> np.ndarray:
        """Return (Z, Y, X) float32 adhesion probability volume.

        Processes in z-chunks so the feature matrix stays ~6 GB per chunk
        regardless of volume size. gate_indices=None means predict everywhere;
        indices are generated per chunk to avoid a 4 GB index array.

        Each scale uses its own minimum z-padding (4σ) rather than the global
        worst-case pad, so small-scale filters process far fewer voxels.
        """
        import time as _tm

        shape = vol_f32.shape
        proba_vol = np.zeros(shape, dtype=np.uint8)

        vox_min = min(voxel_um) if hasattr(voxel_um, "__len__") else float(voxel_um)
        # Minimum padding per scale instead of one global worst-case value.
        scale_pads = [int(np.ceil(4.0 * s / vox_min)) + 1 for s in SCALES_UM]
        max_pad = max(scale_pads)

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
            X = np.empty((N, N_FEATURES), dtype=np.float32)

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
                                       sig, scale_um, X, si * N_PER_SCALE)
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
                                       si * N_PER_SCALE)
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
            del X
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
    """Connected components → optional watershed splitting → (Z, Y, X) int32."""
    from scipy.ndimage import label as _label, distance_transform_edt as _edt

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

        for cid in large:
            # Crop to bounding box so EDT/watershed run on the component only.
            zz, yy, xx = np.where(lbl == cid)
            pad = 2
            z0, z1 = max(0, zz.min() - pad), min(lbl.shape[0], zz.max() + pad + 1)
            y0, y1 = max(0, yy.min() - pad), min(lbl.shape[1], yy.max() + pad + 1)
            x0, x1 = max(0, xx.min() - pad), min(lbl.shape[2], xx.max() + pad + 1)

            crop_mask = (lbl[z0:z1, y0:y1, x0:x1] == cid)
            dist = _edt(crop_mask)
            coords = peak_local_max(dist, min_distance=5, labels=crop_mask)
            if len(coords) <= 1:
                continue

            markers = np.zeros(crop_mask.shape, dtype=np.int32)
            markers[coords[:, 0], coords[:, 1], coords[:, 2]] = np.arange(
                1, len(coords) + 1, dtype=np.int32
            )
            sub = watershed(-dist, markers, mask=crop_mask)

            lbl[z0:z1, y0:y1, x0:x1][crop_mask] = 0
            for sl_id in np.unique(sub[sub > 0]):
                new_max += 1
                lbl[z0:z1, y0:y1, x0:x1][sub == sl_id] = new_max

    # Relabel sequentially
    from skimage.segmentation import relabel_sequential
    lbl, _, _ = relabel_sequential(lbl)
    return lbl.astype(np.int32)
