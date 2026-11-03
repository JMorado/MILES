"""Cartoon (ribbon) meshes for proteins and nucleic acids.

The backbone is a Catmull-Rom spline through the CA atoms (P or C4' for
nucleic acids). It is swept with a superellipse cross-section whose size
depends on the secondary structure: round coils, flat helices and arrow
strands.
"""

from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
from PySide6.QtGui import QColor

from miles.core.molecule import Molecule
from miles.core.protein import NUCLEOTIDES
from miles.render.theme import Theme

SS_COLORS = {"H": "#E4627F", "E": "#F0B43C", "": "#A9B0B9"}
CHAIN_COLORS = ["#5B8DEF", "#2BB673", "#F2994A", "#EB5757", "#9B51E0", "#56CCF2", "#F2C94C",
                "#6FCF97", "#BB6BD9", "#2D9CDB"]

SHAPES = {"": (0.24, 0.24, 1.0), "H": (1.15, 0.24, 1.0), "E": (0.95, 0.26, 0.45)}
"""Cross-section half-width, half-thickness and superellipse exponent per secondary structure."""

ARROW_WIDTH = 1.6
NUCLEIC_RADIUS = 0.42
SAMPLES = 8   # spline points per residue
SIDES = 14    # points around the cross-section


class CartoonColor(StrEnum):
    """How cartoons are coloured."""

    SECONDARY_STRUCTURE = "ss"
    CHAIN = "chain"
    RAINBOW = "rainbow"


@dataclass
class Mesh:
    """Triangle mesh of a cartoon.

    Attributes
    ----------
    vertices, normals, colors : numpy.ndarray
        Per-vertex data, shape ``(n, 3)``; ``colors`` may be tinted.
    indices : numpy.ndarray
        Triangles, shape ``(m, 3)``.
    residue_atoms : list of list of int
        Atom indices of each residue in the mesh.
    vertex_residue : numpy.ndarray
        Residue (index into ``residue_atoms``) of each vertex.
    base_colors : numpy.ndarray
        Untinted vertex colours.
    centers, center_colors : numpy.ndarray
        Spline points and their colours (the QPainter fallback draws these).
    value_range : float
        For surfaces coloured by a value, such as the electrostatic
        potential, the value that gets the full colours; 0 otherwise.
    """

    vertices: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    normals: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    colors: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    indices: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.uint32))
    residue_atoms: list[list[int]] = field(default_factory=list)
    vertex_residue: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    base_colors: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    centers: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    center_colors: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), np.float32))
    value_range: float = 0.0

    @property
    def empty(self) -> bool:
        """Whether the mesh has no triangles."""
        return len(self.indices) == 0


@dataclass
class _Residue:
    chain: str
    ss: str
    nucleic: bool
    atoms: list[int] = field(default_factory=list)
    guide: int | None = None    # atom the spline passes through
    orient: int | None = None   # carbonyl O, sets the ribbon's side direction
    chain_index: int = 0
    chain_fraction: float = 0.0


def _rgb(hex_color: str) -> np.ndarray:
    c = QColor(hex_color)
    return np.array([c.redF(), c.greenF(), c.blueF()], np.float32)


def _unit_rows(v: np.ndarray) -> np.ndarray:
    return v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-9)


def cartoon_residues(mol: Molecule, mask: np.ndarray) -> list[_Residue]:
    """Return the residues whose atoms are in ``mask``, with the atoms their backbone trace runs through.

    This depends on the atoms and residues, not on their positions, so it
    can be kept while atoms only move (pass it to :func:`cartoon_mesh`).
    """
    residues: dict[tuple[str, int, str], _Residue] = {}
    idx = np.flatnonzero(mask & (mol.column("resname") != ""))
    fields = ("chain", "resseq", "resname", "ss", "name")
    rows = zip(idx.tolist(), *(mol.column(c)[idx].tolist() for c in fields), strict=True)
    for i, chain, resseq, resname, ss, name in rows:
        key = (chain, resseq, resname)
        r = residues.get(key)
        if r is None:
            r = residues[key] = _Residue(chain, ss, resname in NUCLEOTIDES)
        r.atoms.append(i)
        if r.nucleic:
            if name == "P" or (name == "C4'" and r.guide is None):
                r.guide = i
        elif name == "CA":
            r.guide = i
        elif name == "O":
            r.orient = i
    return list(residues.values())


def _segments(mol: Molecule, residues: list[_Residue]) -> list[list[_Residue]]:
    """Unbroken backbone traces: chains split where consecutive guide atoms are too far apart."""
    segments: list[list[_Residue]] = []
    current: list[_Residue] = []
    last: _Residue | None = None
    for r in residues:
        if r.guide is None:
            continue
        if current and last is not None and last.guide is not None:
            limit = 8.0 if r.nucleic else 4.4
            gap = np.linalg.norm(mol.coords[r.guide] - mol.coords[last.guide])
            if r.chain != last.chain or gap > limit:
                segments.append(current)
                current = []
        current.append(r)
        last = r
    if current:
        segments.append(current)
    return [s for s in segments if len(s) >= 2]


def _number_chains(segments: list[list[_Residue]]) -> None:
    """Set each residue's chain index and its position along the chain (0 to 1)."""
    length: dict[str, int] = {}
    for r in (r for seg in segments for r in seg):
        length[r.chain] = length.get(r.chain, 0) + 1
    index: dict[str, int] = {}
    position: dict[str, int] = {}
    for r in (r for seg in segments for r in seg):
        r.chain_index = index.setdefault(r.chain, len(index))
        pos = position.get(r.chain, 0)
        r.chain_fraction = pos / max(1, length[r.chain] - 1)
        position[r.chain] = pos + 1


def _residue_color(r: _Residue, coloring: str, theme: Theme | None) -> np.ndarray:
    if coloring == CartoonColor.CHAIN:
        return _rgb(CHAIN_COLORS[r.chain_index % len(CHAIN_COLORS)])
    if coloring == CartoonColor.RAINBOW:
        c = QColor.fromHsvF(max(0.0, 0.66 * (1 - r.chain_fraction)), 0.62, 0.92)
        return np.array([c.redF(), c.greenF(), c.blueF()], np.float32)
    key = r.ss if r.ss in SS_COLORS else ""
    color = _rgb(SS_COLORS[key])
    if not key and theme is not None and not theme.dark:
        color = color * np.float32(0.86)
    return color


def _catmull_rom(P: np.ndarray, samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Points and tangents of a Catmull-Rom spline through the rows of ``P``."""
    ext = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    t = np.linspace(0.0, 1.0, samples, endpoint=False)
    t2, t3 = t * t, t * t * t
    pts, tans = [], []
    for k in range(len(P) - 1):
        p0, p1, p2, p3 = ext[k], ext[k + 1], ext[k + 2], ext[k + 3]
        a = -p0 + p2
        b = 2 * p0 - 5 * p1 + 4 * p2 - p3
        c = -p0 + 3 * p1 - 3 * p2 + p3
        pts.append(0.5 * (2 * p1 + np.outer(t, a) + np.outer(t2, b) + np.outer(t3, c)))
        tans.append(0.5 * (a + np.outer(2 * t, b) + np.outer(3 * t2, c)))
    p0, p1, p2, p3 = ext[-4], ext[-3], ext[-2], ext[-1]
    pts.append(P[-1][None, :])
    end_tangent = 0.5 * ((-p0 + p2) + 2 * (2 * p0 - 5 * p1 + 4 * p2 - p3) + 3 * (-p0 + 3 * p1 - 3 * p2 + p3))
    tans.append(end_tangent[None, :])
    return np.vstack(pts), np.vstack(tans)


def _side_vectors(mol: Molecule, seg: list[_Residue], P: np.ndarray) -> np.ndarray:
    """Compute the unit vector across the ribbon at each residue, flipped to avoid twists."""
    n = len(seg)
    sides = np.zeros((n, 3))
    for k, r in enumerate(seg):
        if r.orient is not None:
            sides[k] = mol.coords[r.orient] - P[k]
        else:
            prev, nxt = P[max(k - 1, 0)], P[min(k + 1, n - 1)]
            sides[k] = (prev + nxt) / 2 - P[k]
            if np.linalg.norm(sides[k]) < 1e-3:
                sides[k] = np.cross(nxt - prev, [0.3, 0.5, 0.8])
    sides = _unit_rows(sides)
    for k in range(1, n):
        if sides[k] @ sides[k - 1] < 0:
            sides[k] = -sides[k]
    return sides


def _profile(seg: list[_Residue], ridx: np.ndarray, frac: np.ndarray) -> tuple[np.ndarray, ...]:
    """Half-width, half-thickness and exponent of the cross-section at each spline point."""
    n, K = len(seg), len(ridx)
    if seg[0].nucleic:
        return np.full(K, NUCLEIC_RADIUS), np.full(K, NUCLEIC_RADIUS), np.ones(K)
    ss = [r.ss if r.ss in SHAPES else "" for r in seg]
    a, b, p = np.empty(K), np.empty(K), np.empty(K)
    for k in range(K):
        r, f = ridx[k], frac[k]
        s0, s1 = ss[r], ss[min(r + 1, n - 1)]
        s2 = ss[r + 2] if r + 2 < n else ""
        if s0 == "E" and s1 == "E" and s2 != "E":   # arrow head at the end of a strand
            a[k] = ARROW_WIDTH * (1 - f) + SHAPES[""][0] * f
            b[k], p[k] = SHAPES["E"][1], SHAPES["E"][2]
        elif s0 == s1:
            a[k], b[k], p[k] = SHAPES[s0]
        else:
            g = 1.0 if s0 == "E" else f * f * (3 - 2 * f)
            lo, hi = np.array(SHAPES[s0]), np.array(SHAPES[s1])
            a[k], b[k], p[k] = lo + (hi - lo) * g
    return a, b, p


def cartoon_mesh(mol: Molecule, mask: np.ndarray, coloring: str = CartoonColor.SECONDARY_STRUCTURE,
                 theme: Theme | None = None, residues: list[_Residue] | None = None) -> Mesh:
    """Ribbon mesh for the residues whose atoms are flagged in ``mask``.

    Parameters
    ----------
    mol : Molecule
        Structure with residue information.
    mask : numpy.ndarray
        Boolean mask of the atoms drawn as cartoon.
    coloring : CartoonColor, optional
        Colouring scheme.
    theme : Theme, optional
        Coils are darkened slightly on light themes.
    residues : list, optional
        :func:`cartoon_residues` of the same atoms, when it is at hand.

    Returns
    -------
    Mesh
        The mesh (empty if no residue has a backbone atom).
    """
    segments = _segments(mol, cartoon_residues(mol, mask) if residues is None else residues)
    if not segments:
        return Mesh()
    _number_chains(segments)
    theta = np.linspace(0, 2 * np.pi, SIDES, endpoint=False)
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    verts, norms, cols, tris, vres, centers, ccols = [], [], [], [], [], [], []
    residue_atoms: list[list[int]] = []
    offset = 0
    for seg in segments:
        n = len(seg)
        first_residue = len(residue_atoms)
        residue_atoms.extend(r.atoms for r in seg)
        seg_colors = np.array([_residue_color(r, coloring, theme) for r in seg], np.float32)
        P = mol.coords[np.array([r.guide for r in seg], dtype=int)]
        pts, tans = _catmull_rom(P, SAMPLES)
        T = _unit_rows(tans)
        K = len(pts)
        rpos = np.minimum(np.arange(K) / SAMPLES, n - 1)
        ridx = np.minimum(np.floor(rpos).astype(int), n - 1)
        frac = rpos - ridx

        sides = _side_vectors(mol, seg, P)
        nxt = np.minimum(ridx + 1, n - 1)
        W = sides[ridx] * (1 - frac)[:, None] + sides[nxt] * frac[:, None]
        W = _unit_rows(W - np.sum(W * T, axis=1)[:, None] * T)
        N = _unit_rows(np.cross(T, W))
        a, b, p = _profile(seg, ridx, frac)

        cx = np.sign(cos_t)[None, :] * np.abs(cos_t)[None, :] ** p[:, None]
        cy = np.sign(sin_t)[None, :] * np.abs(sin_t)[None, :] ** p[:, None]
        ring = (pts[:, None, :] + W[:, None, :] * (a[:, None] * cx)[:, :, None]
                + N[:, None, :] * (b[:, None] * cy)[:, :, None])
        nx = np.sign(cos_t)[None, :] * np.abs(cos_t)[None, :] ** (2 - p[:, None]) / a[:, None]
        ny = np.sign(sin_t)[None, :] * np.abs(sin_t)[None, :] ** (2 - p[:, None]) / b[:, None]
        nrm = W[:, None, :] * nx[:, :, None] + N[:, None, :] * ny[:, :, None]
        nrm = nrm / np.maximum(np.linalg.norm(nrm, axis=2, keepdims=True), 1e-9)

        nearest = np.minimum(np.round(rpos).astype(int), n - 1)
        col = seg_colors[nearest]
        verts.append(ring.reshape(-1, 3))
        norms.append(nrm.reshape(-1, 3))
        cols.append(np.repeat(col, SIDES, axis=0))
        vres.append(np.repeat(first_residue + nearest, SIDES))

        # tube surface: two triangles per quad between consecutive rings
        j = np.arange(SIDES)
        k = np.arange(K - 1)
        q0 = (k[:, None] * SIDES + j[None, :]).ravel()
        q1 = (k[:, None] * SIDES + (j[None, :] + 1) % SIDES).ravel()
        ring_start = offset
        tris.append(np.stack([q0, q1, q1 + SIDES], axis=1) + ring_start)
        tris.append(np.stack([q0, q1 + SIDES, q0 + SIDES], axis=1) + ring_start)
        offset += K * SIDES

        # end caps: a fan around the spline end point
        for end, sign in ((0, -1.0), (K - 1, 1.0)):
            center = offset
            verts.append(pts[end][None, :])
            norms.append((T[end] * sign)[None, :])
            cols.append(col[end][None, :])
            vres.append(np.array([first_residue + nearest[end]]))
            offset += 1
            j0 = ring_start + end * SIDES + j
            j1 = ring_start + end * SIDES + (j + 1) % SIDES
            first, second = (j0, j1) if sign > 0 else (j1, j0)
            tris.append(np.stack([np.full(SIDES, center), first, second], axis=1))
        centers.append(pts)
        ccols.append(col)

    base = np.concatenate(cols).astype(np.float32)
    return Mesh(
        vertices=np.concatenate(verts).astype(np.float32),
        normals=np.concatenate(norms).astype(np.float32),
        colors=base.copy(),
        indices=np.concatenate(tris).astype(np.uint32),
        residue_atoms=residue_atoms,
        vertex_residue=np.concatenate(vres).astype(np.int32),
        base_colors=base,
        centers=np.concatenate(centers).astype(np.float32),
        center_colors=np.concatenate(ccols).astype(np.float32),
    )


def tint_selection(mesh: Mesh, selected: set[int], color: QColor) -> np.ndarray:
    """Vertex colours with the residues that contain selected atoms tinted."""
    if mesh.empty or not selected:
        return mesh.base_colors
    hit = np.array([bool(selected.intersection(atoms)) for atoms in mesh.residue_atoms])
    if not hit.any():
        return mesh.base_colors
    out = mesh.base_colors.copy()
    mask = hit[mesh.vertex_residue]
    accent = np.array([color.redF(), color.greenF(), color.blueF()], np.float32)
    out[mask] = out[mask] * 0.45 + accent * 0.55
    return out
