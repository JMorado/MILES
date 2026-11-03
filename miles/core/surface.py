"""Smooth molecular surfaces and the per-atom values used to colour them.

The surface is the level set of a sum of Gaussian densities, one per atom,
in the spirit of the "blobby molecule" models and VMD's QuickSurf. Each
Gaussian equals the iso-value exactly at the atom's van der Waals radius, so
an isolated atom is a sphere of that radius; where atoms are close the
densities add up and the surface bridges small crevices, much as a probe
sphere rolling over the atoms would. A larger probe gives a smoother, more
blobby surface.

The density is evaluated on a regular grid only around atoms: every atom
adds its Gaussian to the grid points within its cutoff, for a chunk of atoms
at a time, with vectorised NumPy. The isosurface is extracted with Surface
Nets: one vertex per grid cell that the surface crosses, placed at the mean
of the edge crossings, then projected onto the level set along the density
gradient; each grid edge crossed by the surface gives a quad between the
four cells around it. Atoms far apart are meshed on separate grids, so
distant molecules do not make one huge grid.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np

from miles.core.elements import VDW_RADII
from miles.core.molecule import AtomLike, Molecule
from miles.core.neighbors import cross_pairs

ISO_VALUE = 1.0
"""Density at an atom's van der Waals radius, and the level of the surface."""

TAIL = 0.01
"""Gaussians are cut where they fall to this fraction of the iso-value and shifted to reach zero there."""

HYDROGEN_WEIGHT = 0.5
"""Hydrogens add half a Gaussian, so the surface looks similar with or without them."""

CHUNK_EVALUATIONS = 4_000_000
"""Atoms are splatted in chunks of about this many grid evaluations."""

KYTE_DOOLITTLE = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8, "GLY": -0.4,
    "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}
"""Hydropathy of the amino acids (Kyte and Doolittle, 1982)."""

_RESIDUE_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "CYX": "CYS",
    "CYM": "CYS", "SEC": "CYS", "MSE": "MET", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "PYL": "LYS",
}

_CHARGED_ATOMS = {
    ("LYS", "NZ"): 1.0,
    ("ARG", "NE"): 1 / 3, ("ARG", "NH1"): 1 / 3, ("ARG", "NH2"): 1 / 3,
    ("HIP", "ND1"): 0.5, ("HIP", "NE2"): 0.5, ("HSP", "ND1"): 0.5, ("HSP", "NE2"): 0.5,
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
}
"""Formal charges of ionisable side chains at neutral pH, shared over the charged atoms."""

_ION_CHARGES = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Mg": 1, "Ca": 1, "Sr": 1, "Ba": 1,
                "Zn": 1, "Mn": 1, "Fe": 1, "Co": 1, "Ni": 1, "Cu": 1, "Cd": 1,
                "F": -1, "Cl": -1, "Br": -1, "I": -1}
"""Signs of free ions (divalent ions are capped at 1, the colouring range)."""

_FULL_ION_CHARGES = {"Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Zn": 2, "Mn": 2, "Fe": 2, "Co": 2, "Ni": 2,
                     "Cu": 2, "Cd": 2}
"""Charges of the divalent ions above, for :func:`ionisable_charges`."""

_HALOGENS = {"F", "Cl", "Br", "I"}
_POLAR = {"N", "O"}


class SurfaceColoring(StrEnum):
    """What the surface colour shows."""

    ELEMENT = "element"
    HYDROPHOBICITY = "hydrophobicity"
    CHARGE = "charge"
    POTENTIAL = "potential"   # electrostatic potential of partial charges (core/potential.py)
    UNIFORM = "uniform"


@dataclass
class SurfaceMesh:
    """Triangle mesh of a molecular surface.

    Attributes
    ----------
    vertices : numpy.ndarray
        Vertex positions in angstrom, shape ``(n, 3)``, float32.
    normals : numpy.ndarray
        Unit normals pointing out of the molecule, shape ``(n, 3)``, float32.
    faces : numpy.ndarray
        Triangles as vertex indices, shape ``(m, 3)``, uint32; counter-clockwise
        when seen from outside.
    atoms : numpy.ndarray
        Index into the molecule of the atom each vertex belongs to, shape ``(n,)``.
    """

    vertices: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    normals: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    faces: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.uint32))
    atoms: np.ndarray = field(default_factory=lambda: np.zeros(0, int))

    @property
    def empty(self) -> bool:
        """Whether the mesh has no triangles."""
        return len(self.faces) == 0


def default_spacing(n_atoms: int) -> float:
    """Grid spacing in angstrom that keeps a surface of ``n_atoms`` atoms fast to compute."""
    if n_atoms < 2000:
        return 0.5
    if n_atoms < 8000:
        return 0.65
    return 0.8


def blobbiness(probe: float) -> float:
    """Gaussian sharpness for a probe radius: smaller values give smoother surfaces."""
    return min(10.0, max(0.8, 2.25 / max(probe, 0.2)))


def molecular_surface(mol: Molecule, atoms: Sequence[int] | None = None, probe: float = 1.4,
                      spacing: float | None = None) -> SurfaceMesh:
    """Smooth molecular surface around some atoms.

    Parameters
    ----------
    mol : Molecule
        Structure.
    atoms : sequence of int, optional
        Atoms to wrap (default: all).
    probe : float, optional
        Radius in angstrom of the probe the surface imitates; larger values
        bridge wider crevices.
    spacing : float, optional
        Grid spacing in angstrom (default: :func:`default_spacing`).

    Returns
    -------
    SurfaceMesh
        The surface; separate molecules give separate closed pieces.
    """
    idx = np.arange(len(mol)) if atoms is None else np.unique(np.asarray(atoms, dtype=int))
    if len(idx) == 0:
        return SurfaceMesh()
    h = default_spacing(len(idx)) if spacing is None else float(spacing)
    sharpness = blobbiness(probe)
    coords = mol.coords[idx].astype(np.float64)
    numbers = mol.numbers[idx]
    radii = VDW_RADII[numbers]
    weights = np.where(numbers == 1, HYDROGEN_WEIGHT, 1.0)
    cutoffs = radii * np.sqrt(1.0 + math.log(1.0 / TAIL) / sharpness)

    pieces = []
    for cluster in _clusters(coords, 2.0 * float(cutoffs.max()) + 2.0 * h):
        piece = _cluster_surface(coords[cluster], radii[cluster], weights[cluster], cutoffs[cluster],
                                 sharpness, h)
        if piece is not None:
            pieces.append(piece)
    if not pieces:
        return SurfaceMesh()
    vertices, normals, faces = [], [], []
    offset = 0
    for v, n, f in pieces:
        vertices.append(v)
        normals.append(n)
        faces.append(f + offset)
        offset += len(v)
    vertex_array = np.concatenate(vertices).astype(np.float32)
    owner = _nearest_atoms(vertex_array, coords, radii, float(cutoffs.max()) + 2.0 * h)
    return SurfaceMesh(vertices=vertex_array, normals=np.concatenate(normals).astype(np.float32),
                       faces=np.concatenate(faces).astype(np.uint32), atoms=idx[owner])


# -- density ------------------------------------------------------------------------------

def _clusters(coords: np.ndarray, size: float) -> list[np.ndarray]:
    """Split atoms into groups whose densities may touch (cells of atoms joined to neighbouring cells)."""
    cells = np.floor((coords - coords.min(axis=0)) / size).astype(np.int64)
    unique_cells, inverse = np.unique(cells, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    keys = [tuple(k) for k in unique_cells.tolist()]
    lookup = {k: n for n, k in enumerate(keys)}
    label = -np.ones(len(keys), dtype=np.int64)
    offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    current = 0
    for start in range(len(keys)):
        if label[start] >= 0:
            continue
        label[start] = current
        stack = [start]
        while stack:
            x, y, z = keys[stack.pop()]
            for dx, dy, dz in offsets:
                other = lookup.get((x + dx, y + dy, z + dz))
                if other is not None and label[other] < 0:
                    label[other] = current
                    stack.append(other)
        current += 1
    atom_label = label[inverse]
    return [np.nonzero(atom_label == c)[0] for c in range(current)]


def _density_grid(coords: np.ndarray, radii: np.ndarray, weights: np.ndarray, cutoffs: np.ndarray,
                  sharpness: float, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian density on a grid around the atoms, and the grid origin."""
    margin = float(cutoffs.max()) + 2.0 * h
    origin = coords.min(axis=0) - margin
    shape = tuple(int(n) for n in np.ceil((coords.max(axis=0) + margin - origin) / h).astype(int) + 1)
    size = shape[0] * shape[1] * shape[2]

    reach = math.ceil(float(cutoffs.max()) / h) + 1
    steps = np.arange(-reach, reach + 1)
    stencil = np.stack(np.meshgrid(steps, steps, steps, indexing="ij"), axis=-1).reshape(-1, 3)
    stencil = stencil[(stencil ** 2).sum(axis=1) * h * h <= (float(cutoffs.max()) + h * math.sqrt(3)) ** 2]
    stencil_xyz = stencil * h
    stencil_sq = (stencil_xyz ** 2).sum(axis=1)

    density = np.zeros(size, dtype=np.float64)
    base = np.rint((coords - origin) / h).astype(np.int64)
    delta = origin + base * h - coords                      # nearest grid point relative to the atom
    chunk = max(1, CHUNK_EVALUATIONS // len(stencil))
    for start in range(0, len(coords), chunk):
        sl = slice(start, start + chunk)
        d2 = (delta[sl] ** 2).sum(axis=1)[:, None] + 2.0 * delta[sl] @ stencil_xyz.T + stencil_sq[None, :]
        r2 = (radii[sl] ** 2)[:, None]
        inside = d2 < (cutoffs[sl] ** 2)[:, None]
        rows, cols = np.nonzero(inside)
        g = np.exp(-sharpness * (d2[rows, cols] / r2[rows, 0] - 1.0))
        values = weights[sl][rows] * ISO_VALUE * (g - TAIL) / (1.0 - TAIL)
        points = base[sl][rows] + stencil[cols]
        flat = (points[:, 0] * shape[1] + points[:, 1]) * shape[2] + points[:, 2]
        density += np.bincount(flat, weights=values, minlength=size)
    return density.reshape(shape), origin


def _cluster_surface(coords: np.ndarray, radii: np.ndarray, weights: np.ndarray, cutoffs: np.ndarray,
                     sharpness: float, h: float) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    density, origin = _density_grid(coords, radii, weights, cutoffs, sharpness, h)
    return surface_nets(density, ISO_VALUE, origin, np.eye(3) * h)


# -- surface nets ---------------------------------------------------------------------------

def surface_nets(values: np.ndarray, level: float, origin: np.ndarray,
                 axes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Triangle mesh of the level set ``values == level`` on a regular grid.

    The region where the values exceed ``level`` counts as inside. Vertices
    are placed by Surface Nets and then moved onto the level set along the
    gradient of the trilinearly interpolated field.

    Parameters
    ----------
    values : numpy.ndarray
        Grid values, shape ``(nx, ny, nz)``; point ``(i, j, k)`` lies at
        ``origin + i * axes[0] + j * axes[1] + k * axes[2]``.
    level : float
        Iso-value.
    origin : numpy.ndarray
        Position of grid point ``(0, 0, 0)``, shape ``(3,)``.
    axes : numpy.ndarray
        Step vectors of the three grid axes as rows, shape ``(3, 3)``; they
        need not be orthogonal.

    Returns
    -------
    tuple of numpy.ndarray or None
        Vertices ``(n, 3)``, unit normals pointing out of the inside region
        ``(n, 3)`` and triangles ``(m, 3)`` counter-clockwise seen from
        outside; None if the level set does not cross the grid. Parts of the
        surface that reach the edge of the grid are left open.
    """
    field_values = (np.asarray(values, dtype=np.float64) - level).astype(np.float32)
    mesh = _surface_nets(field_values)
    if mesh is None:
        return None
    grid_vertices, faces = mesh
    gradient = [np.asarray(g, dtype=np.float32) for g in np.gradient(field_values)]   # per grid step
    grid_vertices = _project(grid_vertices, field_values, gradient)
    axes = np.asarray(axes, dtype=np.float64)
    grid_gradient = np.stack([_interpolate(g, grid_vertices) for g in gradient], axis=1)
    normals = -grid_gradient @ np.linalg.inv(axes).T          # chain rule: world gradient
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-12)
    if np.linalg.det(axes) < 0:                               # a left-handed grid mirrors the winding
        faces = faces[:, ::-1]
    return np.asarray(origin, dtype=np.float64) + grid_vertices @ axes, normals, faces


_CYCLIC = ((0, 1, 2), (1, 2, 0), (2, 0, 1))   # edge axis, then the two other axes, right-handed


def _surface_nets(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    """Vertices (in grid units) and triangles of the zero level set of ``values`` (positive inside)."""
    cell_shape = np.array(values.shape) - 1
    inside = values > 0
    cells, positions, edges = [], [], []
    for a, b, c in _CYCLIC:
        lower = [slice(None)] * 3
        upper = [slice(None)] * 3
        lower[a] = slice(0, -1)
        upper[a] = slice(1, None)
        crossing = inside[tuple(lower)] != inside[tuple(upper)]
        points = np.argwhere(crossing)
        if not len(points):
            continue
        v0 = values[tuple(lower)][crossing]
        v1 = values[tuple(upper)][crossing]
        t = v0 / (v0 - v1)
        position = points.astype(np.float64)
        position[:, a] += t
        for db in (0, 1):
            for dc in (0, 1):
                cell = points.copy()
                cell[:, b] -= db
                cell[:, c] -= dc
                ok = (cell >= 0).all(axis=1) & (cell < cell_shape).all(axis=1)
                cells.append(np.ravel_multi_index(cell[ok].T, cell_shape))
                positions.append(position[ok])
        edges.append((b, c, points, inside[tuple(lower)][crossing]))
    crossed = np.concatenate(cells) if cells else np.zeros(0, np.int64)
    if not len(crossed):
        return None
    # Only the cells the surface crosses are kept, so a large, mostly empty grid costs little
    used, vertex_index = np.unique(crossed, return_inverse=True)
    position = np.concatenate(positions)
    counts = np.bincount(vertex_index, minlength=len(used))
    vertices = np.stack([np.bincount(vertex_index, weights=position[:, k], minlength=len(used))
                         for k in range(3)], axis=1) / counts[:, None]

    def vertex_of_cell(flat: np.ndarray) -> np.ndarray:
        index = np.minimum(np.searchsorted(used, flat), len(used) - 1)
        return np.where(used[index] == flat, index, -1)

    quads = []
    for b, c, points, outward_positive in edges:
        interior = (points[:, b] >= 1) & (points[:, c] >= 1) & (points[:, b] <= cell_shape[b] - 1) \
            & (points[:, c] <= cell_shape[c] - 1)
        p = points[interior]
        corners = []
        for db, dc in ((1, 1), (0, 1), (0, 0), (1, 0)):   # (b-1, c-1), (b, c-1), (b, c), (b-1, c)
            cell = p.copy()
            cell[:, b] -= db
            cell[:, c] -= dc
            corners.append(vertex_of_cell(np.ravel_multi_index(cell.T, cell_shape)))
        quad = np.stack(corners, axis=1)
        flip = ~outward_positive[interior]
        quad[flip] = quad[flip][:, ::-1]
        quads.append(quad)
    quad_array = np.concatenate(quads)
    quad_array = quad_array[(quad_array >= 0).all(axis=1)]
    faces = np.concatenate([quad_array[:, [0, 1, 2]], quad_array[:, [0, 2, 3]]])
    return vertices, faces


def _interpolate(grid: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Trilinear interpolation of a grid at points given in grid units."""
    upper = np.array(grid.shape) - 2
    base = np.clip(np.floor(points).astype(np.int64), 0, upper)
    f = points - base
    result = np.zeros(len(points))
    for dx in (0, 1):
        wx = f[:, 0] if dx else 1.0 - f[:, 0]
        for dy in (0, 1):
            wy = f[:, 1] if dy else 1.0 - f[:, 1]
            for dz in (0, 1):
                wz = f[:, 2] if dz else 1.0 - f[:, 2]
                result += wx * wy * wz * grid[base[:, 0] + dx, base[:, 1] + dy, base[:, 2] + dz]
    return result


def _project(points: np.ndarray, values: np.ndarray, gradient: list[np.ndarray],
             steps: int = 2) -> np.ndarray:
    """Move points onto the zero level set with Newton steps along the gradient (grid units)."""
    for _ in range(steps):
        f = _interpolate(values, points)
        g = np.stack([_interpolate(component, points) for component in gradient], axis=1)
        g2 = np.maximum((g ** 2).sum(axis=1), 1e-12)
        move = -(f / g2)[:, None] * g
        length = np.linalg.norm(move, axis=1, keepdims=True)
        points = points + move * np.minimum(1.0, 0.5 / np.maximum(length, 1e-12))   # at most half a cell
    return points


def _nearest_atoms(vertices: np.ndarray, coords: np.ndarray, radii: np.ndarray, cutoff: float) -> np.ndarray:
    """Index (into ``coords``) of the atom whose van der Waals sphere is closest to each vertex.

    Atoms are looked for within ``cutoff`` of each vertex. A short search
    goes first: a sphere it finds less than ``reach - largest radius`` away
    is certainly the closest, since every atom further than ``reach`` is
    further than that from the vertex. Only the remaining vertices are
    searched up to ``cutoff``.
    """
    owner = np.zeros(len(vertices), dtype=np.int64)
    best = np.full(len(vertices), np.inf)
    largest = float(radii.max())
    todo = np.arange(len(vertices))
    for reach in (min(cutoff, largest + 1.5), cutoff):
        vi, aj = cross_pairs(vertices[todo], coords, reach)
        if len(vi):
            score = np.linalg.norm(vertices[todo[vi]] - coords[aj], axis=1) - radii[aj]
            order = np.lexsort((score, vi))
            first_vertex, first = np.unique(vi[order], return_index=True)
            rows = todo[first_vertex]
            owner[rows] = aj[order][first]
            best[rows] = score[order][first]
        todo = todo[best[todo] > reach - largest]
        if not len(todo):
            break
    for k in todo[np.isinf(best[todo])]:   # nothing within the cutoff
        owner[k] = int(np.argmin(np.linalg.norm(coords - vertices[k], axis=1) - radii))
    return owner


# -- colouring ------------------------------------------------------------------------------

def atom_property(mol: Molecule, coloring: SurfaceColoring) -> np.ndarray:
    """Per-atom value in [-1, 1] for colouring a surface.

    Parameters
    ----------
    mol : Molecule
        Structure.
    coloring : SurfaceColoring
        ``HYDROPHOBICITY``: positive for hydrophobic atoms, from the Kyte and
        Doolittle scale for amino acids and from the element and its
        neighbours otherwise. ``CHARGE``: formal charges of ionisable groups
        at neutral pH. ``ELEMENT``, ``POTENTIAL`` (computed per vertex, not
        per atom) and ``UNIFORM``: zeros.

    Returns
    -------
    numpy.ndarray
        One value per atom; hydrogens take the value of their parent.
    """
    if coloring == SurfaceColoring.HYDROPHOBICITY:
        values = np.array([_hydrophobicity(mol, i) for i in range(len(mol))])
    elif coloring == SurfaceColoring.CHARGE:
        values = _charges(mol)
    else:
        return np.zeros(len(mol))
    for i, a in enumerate(mol.atoms):
        if a.symbol == "H" and mol.degree(i) == 1:
            values[i] = values[mol.neighbors(i)[0]]
    return np.clip(values, -1.0, 1.0)


def ionisable_charges(mol: Molecule) -> np.ndarray:
    """Charges of ionisable groups at neutral pH, for structures without hydrogens.

    The charges of the ``CHARGE`` colouring (Lys, Arg and His+ side chains,
    Asp and Glu, chain termini, carboxylates, ammonium, phosphates and
    sulfonates of other residues), with free ions at their full charge.
    They need no hydrogens, so they suit crystal structures, where partial
    charges would be wrong.

    Parameters
    ----------
    mol : Molecule
        Structure.

    Returns
    -------
    numpy.ndarray
        Charge of each atom in elementary charges; hydrogens and neutral
        atoms are 0.
    """
    values = _charges(mol)
    for i, a in enumerate(mol.atoms):
        free = not mol.neighbors(i) and _residue(a.resname) not in KYTE_DOOLITTLE
        if free and a.symbol in _FULL_ION_CHARGES:
            values[i] = _FULL_ION_CHARGES[a.symbol]
    return values


def _residue(resname: str) -> str:
    return _RESIDUE_ALIASES.get(resname, resname)


def _hydrophobicity(mol: Molecule, i: int) -> float:
    a = mol.atoms[i]
    scale = KYTE_DOOLITTLE.get(_residue(a.resname))
    if scale is not None:
        return scale / 4.5
    neighbors = {mol.atoms[j].symbol for j in mol.neighbors(i)}
    if a.symbol in ("C", "S"):
        return 0.0 if neighbors & _POLAR else 0.6
    if a.symbol in _HALOGENS:
        return 0.6
    if a.symbol in _POLAR:
        return -0.8
    if a.symbol == "P":
        return -0.5
    return 0.0


def _charges(mol: Molecule) -> np.ndarray:
    values = np.zeros(len(mol))
    termini: dict[str, tuple[int, int]] = {}
    for i, a in enumerate(mol.atoms):
        values[i] = _CHARGED_ATOMS.get((a.resname, a.name), 0.0)   # neutral variants (ASH, LYN...) are absent
        if _residue(a.resname) in KYTE_DOOLITTLE:
            first, last = termini.get(a.chain, (a.resseq, a.resseq))
            termini[a.chain] = (min(first, a.resseq), max(last, a.resseq))
    for i, a in enumerate(mol.atoms):
        if a.chain not in termini or _residue(a.resname) not in KYTE_DOOLITTLE:
            continue
        first, last = termini[a.chain]
        if a.resseq == first and a.name == "N":
            values[i] += 1.0
        if a.resseq == last and a.name in ("O", "OXT") and _has_atom(mol, a, "OXT"):
            values[i] -= 0.5
    _ligand_charges(mol, values)
    return values


def _has_atom(mol: Molecule, residue_atom: AtomLike, name: str) -> bool:
    """Whether the residue of ``residue_atom`` has an atom called ``name``."""
    return any(b.name == name and b.resseq == residue_atom.resseq and b.chain == residue_atom.chain
               for b in mol.atoms)


def _ligand_charges(mol: Molecule, values: np.ndarray) -> None:
    """Charges of common groups outside amino acids: ions, carboxylates, ammonium, phosphates, sulfonates."""
    for i, a in enumerate(mol.atoms):
        if _residue(a.resname) in KYTE_DOOLITTLE:
            continue
        neighbors = mol.neighbors(i)
        if not neighbors:
            values[i] = _ION_CHARGES.get(a.symbol, 0)
            continue
        heavy = [j for j in neighbors if mol.atoms[j].symbol != "H"]
        bare_oxygens = [j for j in heavy if mol.atoms[j].symbol == "O" and mol.degree(j) == 1]
        if a.symbol == "C" and len(bare_oxygens) == 2:                    # carboxylate
            for j in bare_oxygens:
                values[j] = -0.5
        elif a.symbol == "N" and len(neighbors) == 4:                      # ammonium
            values[i] = 1.0
        elif a.symbol == "P" and len(bare_oxygens) >= 2:                   # phosphates
            for j in bare_oxygens:
                values[j] = -(len(bare_oxygens) - 1) / len(bare_oxygens)
        elif a.symbol == "S" and len(bare_oxygens) >= 3:                   # sulfonates, sulfates
            for j in bare_oxygens:
                values[j] = -1.0 / len(bare_oxygens)
