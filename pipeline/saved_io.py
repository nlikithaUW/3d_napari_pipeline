"""The single owner of what a ``<stem>_saved`` folder holds.

Filenames, member keys and dtypes for every intermediate the pipeline writes.
numpy + pathlib only: no Qt, no napari, no ``state`` dict, no printing, and **no
voxel size** — the same rules `saved_check.py` states for itself, for the same
reason. An I/O layer that knows about voxels is the first step back toward the
tangle the rewrite exists to avoid.

**Why this is a module.** Each artifact's schema — its filename and its key
names — used to be re-typed at every writer and every reader. ``("mean", "k1",
"k2", "gaussian")`` appeared in nine places across eight files; the four vessel
artifacts each had *two* writers (the Vessel Mesh tab and ``batch_mesh.run_one``)
that were the same expression typed twice, and ``adhesion_instances.npz`` had
three. Nothing made them agree. A key added to one writer and forgotten in the
other yields a folder that loads fine, measures fine, and exports a CSV with one
column quietly blank — the failure `docs/ISSUES.md` [P12] was opened for. The
only safeguard was a comment asking a future editor to remember.

``docs/scratch/persist/owner_check.py`` is what makes that structural rather than
disciplinary: it fails if an ``np.savez``/``np.load`` call or a bare ``.npz``
filename literal appears anywhere outside this module. A second writer becomes
impossible to add by accident, not merely discouraged.

**Names are a table; I/O is hand-written per artifact.** The names are data
because two consumers need them as data rather than as calls: ``batch_mesh``
composes tagged siblings, and ``saved_check.describe_saved`` matches filenames
while deliberately never calling a loader — it reads ``.npy`` headers so a
pre-flight does not inflate hundreds of MB. The *I/O* is not a table: a generic
writer would need, per artifact, compressed-or-not, a dtype cast, a member
filter, a scalar promotion, an optional member, ``allow_pickle``, and two
different fallbacks for one missing key. Seven option axes for twelve artifacts
is a table bigger than the functions, which is exactly the abstraction-for-
hypothetical-cases `CLAUDE.md` names as what sank the previous version. So the
awkward cases keep their logic, and their comments, visibly.

Loaders always return **materialised arrays, never an open ``NpzFile``**:
a ``with`` block inside here would close the archive before the caller indexes
it.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import numpy as np

# --------------------------------------------------------------------------- #
# Names. The one place a `_saved` filename is spelled.
# --------------------------------------------------------------------------- #

ADHESION_PROBA = "adhesion_proba"
ADHESION_INSTANCES = "adhesion_instances"
CLOVER_RUBY_RATIO = "clover_ruby_ratio"
DAPI_LABELS = "dapi_labels"
NUCLEI_SPLIT = "nuclei_split"
FIBROBLAST_COUNT = "fibroblast_count"
CADHERIN_LABELS = "cadherin_labels"
CADHERIN_JUNCTION_PROBA = "cadherin_junction_proba"
VESSEL_MESH = "vessel_mesh"
VESSEL_CURVATURE = "vessel_curvature"
VESSEL_GRAPH = "vessel_graph"
VESSEL_SHEAR = "vessel_shear"
VESSEL_RECON = "vessel_recon"

#: Retired: the in-plane centerline, replaced by the 3D lumen graph
#: (docs/DIRECTION.md §8B). Nothing reads it, but two messages have to name it
#: so an operator holding an old folder is told what the file is.
VESSEL_PATH_CURVATURE_RETIRED = "vessel_path_curvature.npz"

#: Default extension per artifact, so a caller names the artifact and not the
#: container. `filename(VESSEL_MESH)` is "vessel_mesh.npz" and cannot become
#: ".npy" by a typo at one of eleven call sites.
EXT = {
    ADHESION_PROBA: ".npy",
    ADHESION_INSTANCES: ".npz",
    CLOVER_RUBY_RATIO: ".npy",
    DAPI_LABELS: ".npz",
    NUCLEI_SPLIT: ".npz",
    FIBROBLAST_COUNT: ".npz",
    CADHERIN_LABELS: ".npz",
    CADHERIN_JUNCTION_PROBA: ".npy",
    VESSEL_MESH: ".npz",
    VESSEL_CURVATURE: ".npz",
    VESSEL_GRAPH: ".npz",
    VESSEL_SHEAR: ".npz",
    VESSEL_RECON: ".json",
}

#: The four per-vertex surface curvatures, in the order they are written.
#: Moved here from `saved_check.py`, which defined it and was the only module
#: that used it while eight others spelled the tuple out again.
CURV_KEYS = ("mean", "k1", "k2", "gaussian")

#: Artifacts `describe_saved` reports by presence and size only: reading them
#: means decompressing hundreds of MB and nothing in a pre-flight needs their
#: contents.
BULK = ("adhesion_proba.npy", "adhesion_instances.npz", "dapi_labels.npz",
        "cadherin_labels.npz")


def filename(base: str, tag: str | None = None, ext: str | None = None) -> str:
    """``"vessel_mesh.npz"``, or ``"vessel_mesh__ds.npz"`` for a tagged sibling.

    ``ext`` defaults to the artifact's own extension from `EXT`.
    """
    if ext is None:
        ext = EXT[base]
    return f"{base}__{tag}{ext}" if tag else f"{base}{ext}"


def path(saved, base: str, tag: str | None = None,
         ext: str | None = None) -> Path:
    """Where an artifact lives inside a ``<stem>_saved`` folder.

    Tagged siblings live **beside** the canonical files, not in a subfolder:
    `describe_saved` iterates files and `batch_measure` addresses artifacts by
    name, so anything one directory down is invisible to both.
    """
    return Path(saved) / filename(base, tag, ext)


def tag_of(name: str, base: str, ext: str | None = None) -> str | None:
    """The tag in ``name``, ``""`` for the canonical file, ``None`` if unrelated.

    Distinguishing "canonical" from "not this artifact" is what lets a caller
    scan a folder for every variant of one artifact — which is what
    `describe_saved` does for both the meshes and the recon JSONs.
    """
    if ext is None:
        ext = EXT[base]
    if not name.endswith(ext):
        return None
    stem = name[:-len(ext)]
    if stem == base:
        return ""
    prefix = f"{base}__"
    return stem[len(prefix):] if stem.startswith(prefix) else None


# --------------------------------------------------------------------------- #
# Record types. NamedTuples so an optional member is explicit in the signature
# rather than a `.get` at each of several call sites.
# --------------------------------------------------------------------------- #

class Mesh(NamedTuple):
    """A vessel surface.

    ``faces`` is ``None`` for a folder whose npz predates them, and that
    ``None`` is load-bearing: it is what makes `measure_core` warn that
    ``dist_mesh_um`` is measured to every vertex rather than restricted to the
    surface ([P09a](docs/ISSUES.md#p09)). `batch_measure` passes it through for
    exactly that reason; `batch_mesh` substitutes an empty ``(0, 3)`` array
    instead, because it is about to compute curvature and needs an array. Both
    are right for their caller, so this returns the truth and neither is
    unified away — unifying would silence the warning.

    ``normals`` is present only on the Poisson family, which is what
    `nuc_side.mesh_sign_is_luminal` keys off. Use `has_normals` to ask that
    question: it reads the archive directory and inflates nothing.
    """
    verts: np.ndarray
    faces: np.ndarray | None
    normals: np.ndarray | None


class Labels(NamedTuple):
    """A labeled volume plus its per-component voxel counts."""
    labels: np.ndarray
    sizes: np.ndarray


class CadherinLabels(NamedTuple):
    """The VE-cadherin components.

    ``labels`` and ``sizes`` are the whole live half — `batch_mesh`'s
    ``filter_by_size`` is what consumes them. The three PCA fields are no longer
    written ([P55](docs/ISSUES.md#p55)) and come back ``None``; they are still
    read because files that stored them as pickled ``None`` are on disk.
    """
    labels: np.ndarray
    sizes: np.ndarray
    centroids: np.ndarray | None
    normals: np.ndarray | None
    label_ids: np.ndarray | None


class NucleiSplit(NamedTuple):
    """Watershed-split nuclei and their per-nucleus VnTs classification.

    ``min_size_vox`` is 0 for a file written before it was recorded.
    """
    labels: np.ndarray
    sizes: np.ndarray
    ids: np.ndarray
    is_fib: np.ndarray
    min_size_vox: int


# --------------------------------------------------------------------------- #
# Adhesions
# --------------------------------------------------------------------------- #

def save_proba(dest, proba) -> None:
    """The RF probability volume. Uncompressed: it is the largest artifact on
    disk (~300 MB) and is written once per image but read on every re-measure.
    """
    np.save(dest, proba)


def load_proba(src) -> np.ndarray:
    return np.load(src)


def save_instances(dest, instances) -> None:
    np.savez_compressed(dest, instances=np.asarray(instances).astype(np.int32))


def load_instances(src) -> np.ndarray:
    with np.load(src) as d:
        return d["instances"]


def save_ratio(dest, ratio) -> None:
    np.save(dest, ratio)


def load_ratio(src) -> np.ndarray:
    return np.load(src)


# --------------------------------------------------------------------------- #
# Nuclei
# --------------------------------------------------------------------------- #

def save_dapi_labels(dest, labels, sizes) -> None:
    np.savez_compressed(dest, labels=np.asarray(labels).astype(np.int32),
                        sizes=np.asarray(sizes))


def load_dapi_labels(src) -> Labels:
    with np.load(src) as d:
        return Labels(d["labels"], d["sizes"])


def save_nuclei_split(dest, labels, sizes, ids, is_fib, min_size_vox,
                      threshold) -> None:
    """``min_size_vox`` and ``threshold`` are written through ``int``/``float``,
    not ``np.int64``/``np.float64``. That is deliberate: every file on disk was
    written this way, and "improving" the promotion here would silently change
    the stored dtype of an artifact nothing asked to change.
    """
    np.savez_compressed(dest,
                        labels=np.asarray(labels),
                        sizes=np.asarray(sizes),
                        ids=np.asarray(ids),
                        is_fib=np.asarray(is_fib),
                        min_size_vox=int(min_size_vox),
                        threshold=float(threshold))


def load_nuclei_split(src) -> NucleiSplit:
    with np.load(src) as d:
        return NucleiSplit(
            d["labels"], d["sizes"], d["ids"], d["is_fib"],
            int(d["min_size_vox"]) if "min_size_vox" in d.files else 0)


def save_fibroblast_count(dest, total, fibroblast, other, threshold) -> None:
    """Uncompressed — four scalars; ``savez_compressed`` would only add a
    deflate header per member. Matches every file already on disk.
    """
    np.savez(dest, total=int(total), fibroblast=int(fibroblast),
             other=int(other), threshold=float(threshold))


def load_fibroblast_count(src) -> tuple[dict, list[tuple[str, tuple]]]:
    """Returns ``(values, skipped)``.

    The tolerant loop rather than ``int(d[k])`` for each expected key: that
    raises ``TypeError`` on any non-scalar member, which would turn a
    slightly-different file into a hard failure of a whole image. Non-scalar
    members are reported in ``skipped`` for the caller to warn about — the
    wording lives in `saved_check.fibroblast_count_notes` so the GUI and the
    headless path say the same thing.
    """
    values: dict = {}
    skipped: list[tuple[str, tuple]] = []
    with np.load(src) as d:
        for k in d.files:
            v = d[k]
            if np.ndim(v) != 0:
                skipped.append((k, np.shape(v)))
                continue
            values[k] = float(v) if k == "threshold" else int(v)
    return values, skipped


# --------------------------------------------------------------------------- #
# VE-cadherin
# --------------------------------------------------------------------------- #

def save_cadherin_labels(dest, labels, sizes) -> None:
    """Two keys only. The three PCA fields used to be written here as
    ``np.asarray(None)``, which is a 0-d *object* array: numpy pickles it on
    save and then refuses to unpickle it on load, so every file saved without
    them was unreadable ([P55](docs/ISSUES.md#p55)).
    """
    np.savez_compressed(dest, labels=np.asarray(labels).astype(np.int32),
                        sizes=np.asarray(sizes))


def _lazy_in(d, key):
    """A field that is no longer saved: absent, size-0, or (in files written
    before P55) a pickled ``None`` all mean the same thing — not computed."""
    if key not in d.files:
        return None
    a = d[key]
    if a.dtype == object or a.ndim == 0 or a.size == 0:
        return None
    return a


def load_cadherin_labels(src) -> CadherinLabels:
    """``allow_pickle=True`` — the ONLY loader here that sets it.

    Files written before P55 stored the three PCA fields as pickled ``None``,
    and there are such files on disk (five of the six folders). They are this
    pipeline's own artifacts, and `_lazy_in` throws away anything that is not a
    plain array, so the flag buys backward compatibility and nothing else. It
    must not spread: `load_graph` and `load_shear` pass ``False`` deliberately,
    as a safety property.
    """
    with np.load(src, allow_pickle=True) as d:
        return CadherinLabels(d["labels"], d["sizes"],
                              _lazy_in(d, "centroids"), _lazy_in(d, "normals"),
                              _lazy_in(d, "label_ids"))


def save_junction_proba(dest, proba) -> None:
    np.save(dest, proba)


def load_junction_proba(src) -> np.ndarray:
    return np.load(src)


# --------------------------------------------------------------------------- #
# Vessel geometry
# --------------------------------------------------------------------------- #

def save_mesh(dest, verts, faces, normals=None) -> None:
    """Uncompressed — a multi-million-vertex float64 array, and deflate buys
    little on one. Matches every file already on disk.

    ``normals`` is an additive key: an older reader that names verts/faces is
    unaffected, and a Poisson mesh keeps the normals it cannot re-derive (it is
    open at the volume edges, so `vertex_outward_normals` has no enclosed volume
    to read the global sign from).
    """
    arrays = {"verts": np.asarray(verts), "faces": np.asarray(faces)}
    if normals is not None:
        arrays["normals"] = np.asarray(normals)
    np.savez(dest, **arrays)


def load_mesh(src) -> Mesh:
    with np.load(src) as d:
        return Mesh(d["verts"],
                    d["faces"] if "faces" in d.files else None,
                    d["normals"] if "normals" in d.files else None)


def has_normals(src) -> bool:
    """Does this mesh carry stored normals? Header-only.

    ``.files`` comes from the zip central directory, so this is O(1) in the
    payload. `nuc_side.mesh_sign_is_luminal` needs exactly this and must NOT
    call `load_mesh` — that would inflate a multi-million-vertex mesh to answer
    a yes/no. ``allow_pickle`` is irrelevant here because no member is ever
    indexed.
    """
    with np.load(src, allow_pickle=False) as d:
        return "normals" in d.files


def save_curvature(dest, curv) -> None:
    np.savez_compressed(dest, **{k: np.asarray(curv[k]) for k in CURV_KEYS})


def load_curvature(src) -> dict:
    """Guarded: a partial file yields the keys it has rather than raising.

    This is the unified behaviour. The GUI used to read all four unguarded, so a
    partial file made it load *nothing* while the headless path loaded three of
    four — a divergence over a file neither could fully use.
    `saved_check.curvature_notes` already has the per-key wording for what comes
    back short.
    """
    with np.load(src) as d:
        return {k: d[k] for k in CURV_KEYS if k in d.files}


#: The graph's Python scalars, with the dtype each is stored as. They are
#: promoted explicitly because the array filter below would otherwise drop
#: them — which is what it silently did until [P57](docs/ISSUES.md#p57).
#: Pinned dtypes, not bare ``int``/``float``: ``np.savez`` stores a Python int
#: as the platform ``int_``, so leaving them bare would make the stored dtype
#: depend on which machine wrote the file.
GRAPH_SCALARS = {"n_junctions": np.int64, "n_endpoints": np.int64,
                 "n_crop_endpoints": np.int64, "n_terminus_endpoints": np.int64,
                 "total_length_um": np.float64, "lumen_volume_um3": np.float64}


def save_graph(dest, graph) -> None:
    """Arrays only, plus the scalars in `GRAPH_SCALARS`; the notes and nested
    info dicts are provenance and go to ``vessel_recon.json``, which is where
    provenance already lives.

    The six scalars used to be dropped on every save ([P57](docs/ISSUES.md#p57))
    because the array filter does not match a Python ``int``. The GUI refilled
    them with 0 on load and the headless reader did not, so one file yielded two
    different dicts on the two paths — and `saved_check.graph_notes` records
    what that cost: "a ``graph.get("n_junctions", 0)`` default of 0 silently
    reported 'no junction' on a branched vessel loaded from disk."
    `save_shear` had solved the identical problem correctly one function away.

    A scalar the caller does not supply is simply not written, rather than
    written as a zero: this must not invent a count either.
    """
    arrays = {k: v for k, v in graph.items() if isinstance(v, np.ndarray)}
    scalars = {k: dt(graph[k]) for k, dt in GRAPH_SCALARS.items()
               if k in graph and not isinstance(graph[k], np.ndarray)}
    np.savez_compressed(dest, **arrays, **scalars)


def load_graph(src) -> dict:
    """``allow_pickle=False`` deliberately, as a safety property.

    Fabricates no defaults, including for the `GRAPH_SCALARS` — every file
    written before P57 lacks them, and a 0 standing in for "not recorded" is
    exactly the silently-wrong count that issue is about. `graph_notes` keeps
    counting from the node arrays for that reason, which is what keeps those
    older folders correct.
    """
    with np.load(src, allow_pickle=False) as d:
        return {k: d[k] for k in d.files}


def save_shear(dest, wall) -> None:
    """The ragged wall-shear cloud — one row per wall point, not per centerline
    point, which is why it is a file of its own.

    ``max_dist_um`` is promoted explicitly. The array filter alone would drop
    it, the way `save_graph` still drops its six scalars: this is the same
    problem solved correctly one function away.
    """
    np.savez_compressed(dest, **{k: v for k, v in wall.items()
                                 if isinstance(v, np.ndarray)},
                        max_dist_um=np.float64(wall["max_dist_um"]))


def load_shear(src) -> dict:
    """``allow_pickle=False`` deliberately. ``max_dist_um`` comes back as a
    Python float, which is what every consumer passes on to `shear.vertex_field`.
    """
    with np.load(src, allow_pickle=False) as d:
        wall = {k: d[k] for k in d.files}
    wall["max_dist_um"] = float(wall["max_dist_um"])
    return wall
