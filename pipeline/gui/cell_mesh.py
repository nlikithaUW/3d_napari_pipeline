"""Cell mesh tab — alpha-shape over DAPI centroids + suggest button."""

from __future__ import annotations

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, PushButton
from napari.qt.threading import thread_worker
from scipy.spatial import cKDTree

from ..cell_mesh import (
    alpha_shape_faces,
    nuclei_vertices,
    prune_by_normal,
    signed_distance_to_mesh_vertices,
    suggest_alpha_um,
    tet_circumradii,
    vertex_outward_normals,
)
from .widgets import Histogram, finite_surface, finite_vectors


def _section_segments(verts: np.ndarray, faces: np.ndarray, plane_z: float) -> np.ndarray:
    """Line segments where the mesh intersects z = plane_z. Vectorized.

    Returns (M, 2, 3) in (z, y, x) — both endpoints have z = plane_z.
    """
    if len(faces) == 0:
        return np.zeros((0, 2, 3), dtype=np.float32)
    z = verts[:, 0]
    fz = z[faces]                          # (M, 3)
    above = fz > plane_z
    n_above = above.sum(1)
    cross = (n_above == 1) | (n_above == 2)
    if not cross.any():
        return np.zeros((0, 2, 3), dtype=np.float32)
    f = faces[cross]
    sa = above[cross]
    # The "lone" vertex is the one with the unique sign — True if n_above==1,
    # False if n_above==2. The two edges from the lone vertex are the ones
    # that cross the plane.
    n1 = sa.sum(1) == 1
    lone = np.where(n1[:, None], sa, ~sa)
    idx0 = np.argmax(lone, axis=1)
    rows = np.arange(len(f))
    A = f[rows, idx0]
    B = f[rows, (idx0 + 1) % 3]
    C = f[rows, (idx0 + 2) % 3]
    pA, pB, pC = verts[A], verts[B], verts[C]
    zA = pA[:, 0]
    tAB = (plane_z - zA) / (pB[:, 0] - zA)
    tAC = (plane_z - zA) / (pC[:, 0] - zA)
    iAB = pA + tAB[:, None] * (pB - pA)
    iAC = pA + tAC[:, None] * (pC - pA)
    segs = np.stack([iAB, iAC], axis=1).astype(np.float32)
    segs[:, :, 0] = plane_z
    # Drop degenerate segments before they reach the Shapes layer: a coplanar
    # edge gives a 0/0 intersection (NaN/inf endpoints), and a lone vertex
    # sitting on the plane gives a zero-length segment. Either makes vispy's
    # line frame divide by zero and crashes the GL draw with an access
    # violation. Keep only finite, non-zero-length segments.
    finite = np.isfinite(segs).all(axis=(1, 2))
    length = np.linalg.norm(segs[:, 1] - segs[:, 0], axis=1)
    good = finite & (length > 1e-6)
    return segs[good]


def _keep_mask(sizes: np.ndarray, lids_all: np.ndarray, min_size: int) -> np.ndarray:
    """Boolean mask over lids_all: keep nucleus iff its component size >= min_size.

    Indexing trick: keep_bool[k] = True iff label k passes; index by lids_all.
    """
    keep_bool = np.zeros(sizes.size + 2, dtype=bool)
    keep_bool[1:sizes.size + 1] = sizes >= min_size
    return keep_bool[lids_all]


def build(state, viewer) -> Container:
    hist = Histogram("Tet circumradii (alpha_um)", "log10(circumradius µm)", "seagreen")
    state.setdefault("hists", []).append(hist)

    SECTION_NAME = "Cell mesh (slice)"

    def _update_section(*_):
        mc = state.get("mesh_cache") or {}
        verts = mc.get("verts")
        faces = mc.get("faces")
        if verts is None or faces is None or len(faces) == 0:
            if SECTION_NAME in viewer.layers:
                viewer.layers[SECTION_NAME].data = np.zeros((0, 2, 3), dtype=np.float32)
            return
        # World-z of current slice (mesh verts are in µm; image scale = µm/voxel).
        scale = state.get("scale", (1.0, 1.0, 1.0))
        step_z = int(viewer.dims.current_step[0])
        plane_z = step_z * float(scale[0])
        segs = _section_segments(verts, faces, plane_z)
        if SECTION_NAME in viewer.layers:
            layer = viewer.layers[SECTION_NAME]
            layer.data = segs if len(segs) else np.zeros((0, 2, 3), dtype=np.float32)
        else:
            data = segs if len(segs) else np.zeros((0, 2, 3), dtype=np.float32)
            viewer.add_shapes(
                data, shape_type="line", name=SECTION_NAME,
                edge_color="cyan", edge_width=0.4,
                blending="additive",
            )

    # Connect once; survives mesh rebuilds. Idempotent guard via state flag.
    if not state.get("_section_hooked"):
        viewer.dims.events.current_step.connect(_update_section)
        state["_section_hooked"] = True

    @magicgui(
        call_button="Build cell mesh",
        alpha_um={"widget_type": "FloatSlider", "min": 5.0, "max": 80.0, "step": 0.5},
        min_abs_dot={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0, "step": 0.05},
        inward_offset_um={"widget_type": "FloatSlider", "min": 0.0, "max": 15.0, "step": 0.25},
        show_normals={"widget_type": "CheckBox"},
    )
    def mesh_controls(
        alpha_um: float = 25.0,
        min_abs_dot: float = 0.0,
        inward_offset_um: float = 0.0,
        show_normals: bool = False,
    ):
        labels = state["dapi_state"]["labels"]
        if labels is None:
            print("Run DAPI Otsu first.", flush=True)
            return
        scale = state["scale"]
        nuc_cache = state["nuc_cache"]

        dapi_size = state["dapi_size_controls"]
        min_size = int(dapi_size.min_size.value)
        sizes = state["dapi_state"]["sizes"]
        cache_valid = nuc_cache["labels_id"] == id(labels)

        @thread_worker
        def _run():
            if cache_valid:
                centroids_all = nuc_cache["centroids"]
                normals_all = nuc_cache["normals"]
                lids_all = nuc_cache["lids"]
            else:
                centroids_all, normals_all, lids_all = nuclei_vertices(labels, scale)
                nuc_cache.update(
                    labels_id=id(labels),
                    centroids=centroids_all,
                    normals=normals_all,
                    lids=lids_all,
                )
            if len(centroids_all):
                mask = _keep_mask(sizes, lids_all, min_size)
                centroids = centroids_all[mask]
                normals = normals_all[mask]
            else:
                centroids = centroids_all
                normals = normals_all
            faces = alpha_shape_faces(centroids, alpha_um)
            faces = prune_by_normal(centroids, normals, faces, min_abs_dot)
            if len(faces):
                outward = vertex_outward_normals(centroids, faces)
            else:
                outward = np.zeros_like(centroids)
            offset_centroids = (centroids - inward_offset_um * outward
                                if inward_offset_um > 0 and len(faces)
                                else centroids)
            return centroids, normals, faces, offset_centroids, outward

        def _done(result):
            centroids, normals, faces, offset_centroids, outward = result
            tree = cKDTree(offset_centroids) if len(offset_centroids) else None
            state["mesh_cache"] = {
                "offset_centroids": offset_centroids,
                "outward": outward,
                "inward_offset_um": float(inward_offset_um),
                "kdtree": tree,
                "verts": offset_centroids,
                "faces": faces,
            }
            for name in ("Cell mesh", "Cell mesh (slice)",
                         "Nuclei normals", "Nuclei vertices"):
                if name in viewer.layers:
                    viewer.layers.remove(name)
            if len(centroids) == 0:
                print("No nuclei survived filtering.", flush=True)
                _update_section()
                return
            if len(faces) == 0:
                print(f"{len(centroids)} vertices but 0 faces — increase alpha_um.",
                      flush=True)
                viewer.add_points(
                    offset_centroids, name="Nuclei vertices", size=2.0,
                    face_color="cyan", blending="additive",
                )
                _update_section()
                return
            verts_disp = offset_centroids
            vals = np.linalg.norm(verts_disp - verts_disp.mean(0), axis=1)
            viewer.add_surface(
                finite_surface(verts_disp, faces, vals), name="Cell mesh",
                colormap="cividis", blending="translucent", opacity=0.6,
            )
            _update_section()
            if show_normals:
                tail = verts_disp
                head = verts_disp + 4.0 * outward
                vectors = finite_vectors(np.stack([tail, head - tail], axis=1))
                viewer.add_vectors(
                    vectors, name="Nuclei normals", edge_color="yellow",
                    edge_width=0.3, blending="additive",
                )
            print(
                f"Cell mesh: {len(centroids)} vertices, {len(faces)} faces "
                f"(alpha_um={alpha_um}, min_abs_dot={min_abs_dot}, "
                f"inward_offset_um={inward_offset_um}).",
                flush=True,
            )

        print(
            f"Building cell mesh (alpha_um={alpha_um}, "
            f"inward_offset_um={inward_offset_um}) ...",
            flush=True,
        )
        w = _run()
        w.returned.connect(_done)
        w.start()

    suggest_w = PushButton(text="Suggest alpha_um from histogram")

    def _suggest(*_):
        labels = state["dapi_state"]["labels"]
        if labels is None:
            print("Run DAPI Otsu first.", flush=True)
            return
        scale = state["scale"]
        nuc_cache = state["nuc_cache"]
        dapi_size = state["dapi_size_controls"]
        min_size = int(dapi_size.min_size.value)
        sizes = state["dapi_state"]["sizes"]
        cache_valid = nuc_cache["labels_id"] == id(labels)

        @thread_worker
        def _run():
            if cache_valid:
                centroids_all = nuc_cache["centroids"]
                normals_all = nuc_cache["normals"]
                lids_all = nuc_cache["lids"]
            else:
                centroids_all, normals_all, lids_all = nuclei_vertices(labels, scale)
                nuc_cache.update(
                    labels_id=id(labels),
                    centroids=centroids_all,
                    normals=normals_all,
                    lids=lids_all,
                )
            if len(centroids_all):
                mask = _keep_mask(sizes, lids_all, min_size)
                centroids = centroids_all[mask]
            else:
                centroids = centroids_all
            R = tet_circumradii(centroids)
            suggested = suggest_alpha_um(R)
            return centroids, R, suggested

        def _done(result):
            centroids, R, suggested = result
            if R.size == 0:
                print(f"Only {len(centroids)} centroids — need ≥4 to "
                      "tetrahedralize.", flush=True)
                hist.clear()
                return
            print(
                f"alpha histogram: {R.size} tets, "
                f"R range {R[np.isfinite(R)].min():.2f}–"
                f"{R[np.isfinite(R)].max():.2f} µm, "
                f"suggested alpha_um = {suggested:.2f}.", flush=True,
            )
            R_plot = R[np.isfinite(R) & (R > 0)]
            if R_plot.size:
                cap = np.percentile(R_plot, 99.0)
                R_plot = R_plot[R_plot <= cap]
            hist.update(R_plot, threshold=suggested if suggested > 0 else None,
                        threshold_label=f"alpha_um={suggested:.2f}")
            if suggested > 0:
                if suggested > mesh_controls.alpha_um.max:
                    mesh_controls.alpha_um.max = float(suggested) * 2.0
                mesh_controls.alpha_um.value = suggested

        print("Computing tet circumradii ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.start()

    suggest_w.changed.connect(_suggest)

    CADH_VERTS_NAME = "Cadherin vertices"

    def _rebuild_layers(centroids, normals, faces, offset_centroids, outward,
                        inward_offset_um, n_cadh, alpha_um, min_abs_dot,
                        show_cadh_pts, cadh_points):
        tree = cKDTree(offset_centroids) if len(offset_centroids) else None
        state["mesh_cache"] = {
            "offset_centroids": offset_centroids,
            "outward": outward,
            "inward_offset_um": float(inward_offset_um),
            "kdtree": tree,
            "verts": offset_centroids,
            "faces": faces,
        }
        for name in ("Cell mesh", "Cell mesh (slice)",
                     "Nuclei normals", "Nuclei vertices", CADH_VERTS_NAME):
            if name in viewer.layers:
                viewer.layers.remove(name)
        if len(centroids) == 0:
            print("No vertices to mesh.", flush=True)
            _update_section()
            return
        if len(faces) == 0:
            print(f"{len(centroids)} vertices but 0 faces — increase alpha_um.",
                  flush=True)
            viewer.add_points(
                offset_centroids, name="Nuclei vertices", size=2.0,
                face_color="cyan", blending="additive",
            )
            _update_section()
            return
        verts_disp = offset_centroids
        vals = np.linalg.norm(verts_disp - verts_disp.mean(0), axis=1)
        viewer.add_surface(
            (verts_disp, faces, vals), name="Cell mesh",
            colormap="cividis", blending="translucent", opacity=0.6,
        )
        _update_section()
        if show_cadh_pts and cadh_points is not None and len(cadh_points):
            viewer.add_points(
                cadh_points, name=CADH_VERTS_NAME, size=2.0,
                face_color="magenta", blending="additive",
            )
        print(
            f"Cell mesh: {len(centroids)} vertices ({n_cadh} cadherin), "
            f"{len(faces)} faces (alpha_um={alpha_um}, "
            f"min_abs_dot={min_abs_dot}, inward_offset_um={inward_offset_um}).",
            flush=True,
        )

    @magicgui(
        call_button="Refine with cadherin",
        refine_band_um={"widget_type": "FloatSlider", "min": 0.5, "max": 50.0, "step": 0.5},
        min_sep_um={"widget_type": "FloatSlider", "min": 0.0, "max": 30.0, "step": 0.25},
        show_cadherin_vertices={"widget_type": "CheckBox"},
    )
    def refine_controls(
        refine_band_um: float = 5.0,
        min_sep_um: float = 5.0,
        show_cadherin_vertices: bool = True,
    ):
        mesh = state.get("mesh_cache") or {}
        if mesh.get("kdtree") is None or mesh.get("offset_centroids") is None \
                or len(mesh["offset_centroids"]) == 0:
            print("Build DAPI cell mesh first.", flush=True)
            return
        cs = state.get("cadh_state") or {}
        cadh_cents_all = cs.get("centroids")
        cadh_norms_all = cs.get("normals")
        keep_mask = cs.get("keep_mask")
        if cadh_cents_all is None or keep_mask is None or len(cadh_cents_all) == 0:
            print("Run the Cadherin tab first.", flush=True)
            return

        labels = state["dapi_state"]["labels"]
        if labels is None:
            print("Run DAPI Otsu first.", flush=True)
            return
        scale = state["scale"]
        nuc_cache = state["nuc_cache"]
        dapi_size = state["dapi_size_controls"]
        min_size = int(dapi_size.min_size.value)
        sizes = state["dapi_state"]["sizes"]
        cache_valid = nuc_cache["labels_id"] == id(labels)

        alpha_um = float(mesh_controls.alpha_um.value)
        min_abs_dot = float(mesh_controls.min_abs_dot.value)
        inward_offset_um = float(mesh_controls.inward_offset_um.value)

        offset_centroids_prev = mesh["offset_centroids"]
        outward_prev = mesh["outward"]
        kdtree_prev = mesh["kdtree"]

        @thread_worker
        def _run():
            if cache_valid:
                centroids_all = nuc_cache["centroids"]
                normals_all = nuc_cache["normals"]
                lids_all = nuc_cache["lids"]
            else:
                centroids_all, normals_all, lids_all = nuclei_vertices(labels, scale)
                nuc_cache.update(
                    labels_id=id(labels),
                    centroids=centroids_all,
                    normals=normals_all,
                    lids=lids_all,
                )
            if len(centroids_all):
                mask = _keep_mask(sizes, lids_all, min_size)
                dapi_cents = centroids_all[mask]
                dapi_norms = normals_all[mask]
            else:
                dapi_cents = centroids_all
                dapi_norms = normals_all

            cadh_cents = cadh_cents_all[keep_mask]
            cadh_norms = cadh_norms_all[keep_mask]

            if len(cadh_cents):
                sd = signed_distance_to_mesh_vertices(
                    cadh_cents, offset_centroids_prev, outward_prev,
                    tree=kdtree_prev,
                )
                band = np.abs(sd) <= float(refine_band_um)
                cadh_cents = cadh_cents[band]
                cadh_norms = cadh_norms[band]

            if len(cadh_cents) and len(dapi_cents) and min_sep_um > 0:
                t = cKDTree(dapi_cents)
                d, _ = t.query(cadh_cents, k=1)
                far = d > float(min_sep_um)
                cadh_cents = cadh_cents[far]
                cadh_norms = cadh_norms[far]

            if len(cadh_cents):
                centroids = np.vstack([dapi_cents, cadh_cents])
                normals = np.vstack([dapi_norms, cadh_norms])
            else:
                centroids = dapi_cents
                normals = dapi_norms

            faces = alpha_shape_faces(centroids, alpha_um)
            faces = prune_by_normal(centroids, normals, faces, min_abs_dot)
            if len(faces):
                outward = vertex_outward_normals(centroids, faces)
            else:
                outward = np.zeros_like(centroids)
            offset_centroids = (centroids - inward_offset_um * outward
                                if inward_offset_um > 0 and len(faces)
                                else centroids)
            return centroids, normals, faces, offset_centroids, outward, cadh_cents

        def _done(result):
            centroids, normals, faces, offset_centroids, outward, cadh_pts = result
            n_cadh = len(cadh_pts)
            cadh_disp = (offset_centroids[-n_cadh:] if n_cadh else None)
            _rebuild_layers(
                centroids, normals, faces, offset_centroids, outward,
                inward_offset_um, n_cadh, alpha_um, min_abs_dot,
                show_cadherin_vertices, cadh_disp,
            )

        print(f"Refining cell mesh with cadherin "
              f"(refine_band_um={refine_band_um}, min_sep_um={min_sep_um}) ...",
              flush=True)
        w = _run()
        w.returned.connect(_done)
        w.start()

    return (Container(widgets=[mesh_controls, suggest_w, refine_controls], labels=False),
            hist.canvas)
