"""A small valence force field used by Clean to tidy up built structures.

Terms: harmonic bonds, cosine-harmonic angles taken from each atom's
coordination geometry, planarity torsions about bonds between trigonal
atoms, a weak staggering torsion on sp3-sp3 single bonds and a soft
repulsion between non-bonded atoms of the same fragment. The aim is a
sensible starting geometry, not an accurate structure.

All terms are evaluated with vectorised NumPy over index arrays, and the
geometry is minimised with L-BFGS.
"""

import itertools
import math
from collections.abc import Callable, Iterable, Sequence
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from miles.core.chem import minimize as minimize_mmff
from miles.core.chem import reembed
from miles.core.elements import VDW_RADII, bond_length
from miles.core.geometry import IDEAL_ANGLES, Geometry
from miles.core.molecule import Molecule
from miles.core.neighbors import near_mask, pairs_within
from miles.core.perception import infer_atom_types

K_BOND = 400.0
K_ANGLE = 120.0
K_PLANAR = 8.0
K_STAGGER = 0.6
K_REPULSE = 60.0

CONTACT_CUTOFF = 6.0
"""Contacts are searched within this distance for structures over 150 atoms."""

GENTLE_RESTRAINT = 10.0
"""Spring constant tying experimental atoms to their positions in a gentle clean."""

ENVIRONMENT = 6.0
"""Atoms this close to a cleaned selection are held fixed around it."""

MMFF_ELEMENTS = frozenset({"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"})
"""Elements MMFF94 covers; structures with other elements use the built-in force field."""

MMFF_ATOM_LIMIT = 400
"""Larger structures use the built-in force field, which is much faster."""

_LONE_PAIR_ANGLES = {("O", 2): 104.5, ("S", 2): 96.0, ("Se", 2): 92.0, ("N", 3): 107.0, ("P", 3): 98.0}

EnergyGradient = Callable[[np.ndarray], tuple[float, np.ndarray]]


class CleanMode(StrEnum):
    """What :func:`clean` did."""

    MMFF = "mmff94"    # relaxed with MMFF94
    FULL = "full"      # relaxed with the built-in force field
    GENTLE = "gentle"  # experimental atoms restrained to their positions


def _pair_keys(a: np.ndarray, b: np.ndarray, n: int) -> np.ndarray:
    """One integer per unordered atom pair, for fast set membership tests."""
    return np.minimum(a, b).astype(np.int64) * n + np.maximum(a, b)


class ForceField:
    """Energy and gradient of a structure under the Clean force field.

    The terms are set up once from the bonds and atom geometries; the
    energy can then be evaluated for any coordinates.

    Parameters
    ----------
    mol : Molecule
        Structure to set the terms up for (not modified).
    soft : bool, optional
        Leave out torsions and contact repulsion, so a structure far from its
        target (e.g. a chain just closed into a ring) can fold.
    restraint : array_like, optional
        Per-atom spring constants tying atoms to their starting positions.
    """

    def __init__(self, mol: Molecule, soft: bool = False, restraint: ArrayLike | None = None) -> None:
        self.soft = soft
        self.x0 = mol.coords.copy()
        self.restraint = None if restraint is None else np.asarray(restraint, dtype=float)
        work = mol.copy()
        infer_atom_types(work)
        self.n = len(work)
        self._setup_bonds(work)
        self._setup_angles(work)
        self._setup_torsions(work)
        self._setup_contacts(work)

    # -- setup ---------------------------------------------------------------------

    def _setup_bonds(self, mol: Molecule) -> None:
        bonds = list(mol.bonds.items())
        self.b_idx = np.array([key for key, _ in bonds], dtype=int).reshape(-1, 2)
        self.b_len = np.array([bond_length(mol.atoms[i].symbol, mol.atoms[j].symbol, order)
                               for (i, j), order in bonds])

    def _setup_angles(self, mol: Molecule) -> None:
        """Every angle i-j-k gets up to three ideal cosines; the closest one applies."""
        idx, cosines = [], []
        for j, atom in enumerate(mol.atoms):
            nbrs = mol.neighbors(j)
            ideals = IDEAL_ANGLES.get(atom.geom or "")
            if len(nbrs) < 2 or not ideals:
                continue
            lone_pairs = _LONE_PAIR_ANGLES.get((atom.symbol, len(nbrs)))
            if lone_pairs and atom.geom == Geometry.TETRAHEDRAL and (atom.nsites or 0) == len(nbrs):
                ideals = (lone_pairs,)
            for i, k in itertools.combinations(nbrs, 2):
                targets = ideals
                if mol.bond_order(i, k) > 0:  # three-membered ring
                    targets = (60.0,)
                elif set(mol.neighbors(i)) & set(mol.neighbors(k)) - {j}:  # four-membered ring
                    targets = (90.0,)
                cos = [math.cos(math.radians(t)) for t in targets]
                idx.append((i, j, k))
                cosines.append(cos + [cos[-1]] * (3 - len(cos)))
        self.a_idx = np.array(idx, dtype=int).reshape(-1, 3)
        self.a_cos = np.array(cosines, dtype=float).reshape(-1, 3)

    def _setup_torsions(self, mol: Molecule) -> None:
        """Planar torsions about trigonal-trigonal bonds, staggering about sp3-sp3 single bonds."""
        idx, planar_flags = [], []
        for (j, k), order in mol.bonds.items():
            gj, gk = mol.atoms[j].geom, mol.atoms[k].geom
            planar = gj == Geometry.TRIGONAL and gk == Geometry.TRIGONAL
            stagger = order == 1.0 and gj == Geometry.TETRAHEDRAL and gk == Geometry.TETRAHEDRAL
            if not (planar or stagger):
                continue
            for i in mol.neighbors(j):
                if i == k:
                    continue
                for m in mol.neighbors(k):
                    if m not in (j, i):
                        idx.append((i, j, k, m))
                        planar_flags.append(planar)
        self.t_idx = np.array(idx, dtype=int).reshape(-1, 4)
        self.t_planar = np.array(planar_flags, dtype=bool)

    def _setup_contacts(self, mol: Molecule) -> None:
        """Repulsion between atoms of the same molecule that are neither 1-2 nor 1-3.

        Large structures use a neighbour search, so the cost stays linear.
        """
        n = self.n
        excluded = _pair_keys(self.b_idx[:, 0], self.b_idx[:, 1], n)
        excluded = np.concatenate([excluded, _pair_keys(self.a_idx[:, 0], self.a_idx[:, 2], n)])
        one_four = _pair_keys(self.t_idx[:, 0], self.t_idx[:, 3], n)
        component = np.zeros(n, dtype=int)
        for c, members in enumerate(mol.components()):
            component[list(members)] = c
        if n > 150:
            pairs = pairs_within(mol.coords, CONTACT_CUTOFF)
            ii, jj = pairs[:, 0], pairs[:, 1]
        elif n > 1:
            ii, jj = np.triu_indices(n, 1)
        else:
            ii = jj = np.zeros(0, int)
        same = component[ii] == component[jj]
        ii, jj = ii[same], jj[same]
        keys = _pair_keys(ii, jj, n)
        keep = ~np.isin(keys, excluded)
        ii, jj, keys = ii[keep], jj[keep], keys[keep]
        vdw = VDW_RADII[mol.numbers]
        factor = np.where(np.isin(keys, one_four), 0.65, 0.8)
        self.p_idx = np.stack([ii, jj], axis=1).astype(int).reshape(-1, 2)
        self.p_min = (factor * (vdw[ii] + vdw[jj])).astype(float)

    # -- evaluation --------------------------------------------------------------------

    def energy_gradient(self, x: np.ndarray) -> tuple[float, np.ndarray]:
        """Evaluate the energy and its gradient.

        Parameters
        ----------
        x : numpy.ndarray
            Coordinates, shape ``(n, 3)`` or flattened.

        Returns
        -------
        energy : float
            Total energy (arbitrary units).
        gradient : numpy.ndarray
            Derivative with respect to the coordinates, shape ``(n, 3)``.
        """
        x = x.reshape(-1, 3)
        g = np.zeros_like(x)
        e = self._bonds(x, g)
        e += self._angles(x, g)
        if not self.soft:
            e += self._torsions(x, g)
            e += self._contacts(x, g)
        if self.restraint is not None:
            d = x - self.x0
            e += float((self.restraint[:, None] * d * d).sum())
            g += 2 * self.restraint[:, None] * d
        return e, g

    def _bonds(self, x: np.ndarray, g: np.ndarray) -> float:
        if not len(self.b_idx):
            return 0.0
        i, j = self.b_idx[:, 0], self.b_idx[:, 1]
        d = x[i] - x[j]
        r = np.maximum(np.linalg.norm(d, axis=1), 1e-6)
        dr = r - self.b_len
        f = (2 * K_BOND * dr / r)[:, None] * d
        np.add.at(g, i, f)
        np.add.at(g, j, -f)
        return K_BOND * float((dr * dr).sum())

    def _angles(self, x: np.ndarray, g: np.ndarray) -> float:
        if not len(self.a_idx):
            return 0.0
        i, j, k = self.a_idx[:, 0], self.a_idx[:, 1], self.a_idx[:, 2]
        u = x[i] - x[j]
        v = x[k] - x[j]
        lu = np.maximum(np.linalg.norm(u, axis=1), 1e-6)
        lv = np.maximum(np.linalg.norm(v, axis=1), 1e-6)
        c = (u * v).sum(axis=1) / (lu * lv)
        nearest = np.abs(self.a_cos - c[:, None]).argmin(axis=1)
        dc = c - self.a_cos[np.arange(len(c)), nearest]
        dedc = (2 * K_ANGLE * dc)[:, None]
        gu = dedc * (v / (lu * lv)[:, None] - c[:, None] * u / (lu * lu)[:, None])
        gv = dedc * (u / (lu * lv)[:, None] - c[:, None] * v / (lv * lv)[:, None])
        np.add.at(g, i, gu)
        np.add.at(g, k, gv)
        np.add.at(g, j, -gu - gv)
        return K_ANGLE * float((dc * dc).sum())

    def _torsions(self, x: np.ndarray, g: np.ndarray) -> float:
        """Torsion energy, with the gradient formulation of Blondel and Karplus (1996)."""
        if not len(self.t_idx):
            return 0.0
        p1, p2, p3, p4 = (x[self.t_idx[:, m]] for m in range(4))
        F = p1 - p2
        G = p2 - p3
        H = p4 - p3
        A = np.cross(F, G)
        B = np.cross(H, G)
        a2 = np.maximum((A * A).sum(axis=1), 1e-10)
        b2 = np.maximum((B * B).sum(axis=1), 1e-10)
        gl = np.maximum(np.linalg.norm(G, axis=1), 1e-8)
        cos_phi = (A * B).sum(axis=1) / np.sqrt(a2 * b2)
        sin_phi = (np.cross(B, A) * G).sum(axis=1) / (np.sqrt(a2 * b2) * gl)
        phi = np.arctan2(sin_phi, cos_phi)
        energy = np.where(self.t_planar, K_PLANAR * np.sin(phi) ** 2, K_STAGGER * (1 + np.cos(3 * phi)))
        dedphi = np.where(self.t_planar, K_PLANAR * np.sin(2 * phi), -3 * K_STAGGER * np.sin(3 * phi))
        fg = (F * G).sum(axis=1)
        hg = (H * G).sum(axis=1)
        d1 = -(gl / a2)[:, None] * A
        d4 = (gl / b2)[:, None] * B
        d2 = (gl / a2)[:, None] * A + (fg / (a2 * gl))[:, None] * A - (hg / (b2 * gl))[:, None] * B
        d3 = -(gl / b2)[:, None] * B - (fg / (a2 * gl))[:, None] * A + (hg / (b2 * gl))[:, None] * B
        for m, dm in enumerate((d1, d2, d3, d4)):
            np.add.at(g, self.t_idx[:, m], dedphi[:, None] * dm)
        return float(energy.sum())

    def _contacts(self, x: np.ndarray, g: np.ndarray) -> float:
        if not len(self.p_idx):
            return 0.0
        i, j = self.p_idx[:, 0], self.p_idx[:, 1]
        d = x[i] - x[j]
        r = np.maximum(np.linalg.norm(d, axis=1), 1e-6)
        over = np.maximum(self.p_min - r, 0.0)
        f = (-2 * K_REPULSE * over / r)[:, None] * d
        np.add.at(g, i, f)
        np.add.at(g, j, -f)
        return K_REPULSE * float((over * over).sum())


def _lbfgs_direction(g: np.ndarray, s_hist: list[np.ndarray], y_hist: list[np.ndarray]) -> np.ndarray:
    """L-BFGS two-loop recursion: approximate ``-H^-1 g`` from recent steps."""
    q = g.copy()
    alphas = []
    for s, y in zip(reversed(s_hist), reversed(y_hist), strict=True):
        rho = 1.0 / float(y @ s)
        a = rho * float(s @ q)
        q -= a * y
        alphas.append((rho, a))
    if s_hist:
        gamma = float(s_hist[-1] @ y_hist[-1]) / float(y_hist[-1] @ y_hist[-1])
    else:
        gamma = 0.05 / max(1.0, float(np.abs(g).max()))
    r = gamma * q
    for s, y, (rho, a) in zip(s_hist, y_hist, reversed(alphas), strict=True):
        r += s * (a - rho * float(y @ r))
    return -r


def minimize(ff: ForceField, x0: np.ndarray, max_steps: int = 800, gtol: float = 1e-3,
             max_move: float = 0.3, frozen: ArrayLike | None = None, memory: int = 10) -> np.ndarray:
    """Minimise the force-field energy with L-BFGS and a backtracking line search.

    Parameters
    ----------
    ff : ForceField
        Energy function.
    x0 : numpy.ndarray
        Starting coordinates, shape ``(n, 3)``.
    max_steps : int, optional
        Iteration limit.
    gtol : float, optional
        Stop when no gradient component exceeds this.
    max_move : float, optional
        Largest displacement of any coordinate in one step, in angstrom.
    frozen : array_like of bool, optional
        Atoms that stay in place.
    memory : int, optional
        Number of past steps kept by L-BFGS.

    Returns
    -------
    numpy.ndarray
        Minimised coordinates, shape ``(n, 3)``.
    """
    energy_gradient: EnergyGradient = ff.energy_gradient
    if frozen is not None and np.any(frozen):
        movable = np.repeat(~np.asarray(frozen, dtype=bool), 3).astype(float)

        def energy_gradient(y: np.ndarray) -> tuple[float, np.ndarray]:
            e, g = ff.energy_gradient(y)
            return e, g.reshape(-1) * movable

    x = x0.reshape(-1).astype(float).copy()
    e, g = energy_gradient(x)
    g = g.reshape(-1)
    s_hist: list[np.ndarray] = []
    y_hist: list[np.ndarray] = []
    for _ in range(max_steps):
        if np.abs(g).max() < gtol:
            break
        d = _lbfgs_direction(g, s_hist, y_hist)
        if float(d @ g) >= 0:  # not a descent direction: restart along the gradient
            d = -g * 0.05 / max(1.0, float(np.abs(g).max()))
            s_hist.clear()
            y_hist.clear()
        step = np.abs(d.reshape(-1, 3)).max()
        if step > max_move:
            d *= max_move / step
        t = 1.0
        slope = float(d @ g)
        for _ in range(30):
            x_new = x + t * d
            e_new, g_new = energy_gradient(x_new)
            if e_new <= e + 1e-4 * t * slope:
                break
            t *= 0.5
        else:
            break
        g_new = g_new.reshape(-1)
        s, y = x_new - x, g_new - g
        if float(s @ y) > 1e-10:
            s_hist.append(s)
            y_hist.append(y)
            if len(s_hist) > memory:
                s_hist.pop(0)
                y_hist.pop(0)
        x, e, g = x_new, e_new, g_new
    return x.reshape(-1, 3)


def _mmff_coordinates(mol: Molecule, fixed: Sequence[int] | None = None) -> np.ndarray | None:
    """MMFF94 minimum of an organic structure, or None if MMFF94 should not or cannot be used.

    Delocalised bonds outside rings (nitro groups, carboxylates drawn with
    order 1.5) have no Kekule form MMFF94 types correctly, so structures with
    them use the built-in force field.
    """
    if (len(mol) > MMFF_ATOM_LIMIT or not mol.mask("symbol", MMFF_ELEMENTS).all()
            or (mol.column("resname") != "").any()):
        return None
    ring = mol.ring_bonds()
    if any(order == 1.5 and key not in ring for key, order in mol.bonds.items()):
        return None
    return minimize_mmff(mol, fixed, force_fields=("MMFF94",))


def clean_geometry(mol: Molecule, max_steps: int = 800) -> CleanMode:
    """Relax the geometry in place.

    Organic structures of up to :data:`MMFF_ATOM_LIMIT` atoms are minimised
    with MMFF94. Everything else (metals, unusual bonding, large systems) uses
    the built-in force field. Atoms with residue data (from PDB or mmCIF
    files) are tied to their starting positions, so a protein is only tidied,
    not rebuilt.

    Parameters
    ----------
    mol : Molecule
        Structure to relax.
    max_steps : int, optional
        Iteration limit (lower for large or experimental structures).

    Returns
    -------
    CleanMode
        Whether a full or a gentle clean was done.
    """
    if len(mol) < 2:
        return CleanMode.FULL
    experimental = np.asarray(mol.column("resname") != "", dtype=bool)
    if experimental.any():
        ff = ForceField(mol, restraint=np.where(experimental, GENTLE_RESTRAINT, 0.0))
        mol.coords = minimize(ff, mol.coords, max_steps=min(max_steps, 300), gtol=0.05)
        return CleanMode.GENTLE
    coords = _mmff_coordinates(mol)
    if coords is not None:
        mol.coords = coords
        return CleanMode.MMFF
    steps = max_steps if len(mol) < 1000 else 300
    mol.coords = minimize(ForceField(mol), mol.coords, max_steps=steps)
    return CleanMode.FULL


def clean_selection(mol: Molecule, atoms: Iterable[int]) -> None:
    """Fully relax the selected atoms in place.

    Atoms within :data:`ENVIRONMENT` of the selection are included but held
    fixed, so the selection settles into its surroundings.

    Parameters
    ----------
    mol : Molecule
        Structure.
    atoms : iterable of int
        Atoms to relax.
    """
    selected = sorted({int(i) for i in atoms})
    if not selected:
        return
    near = near_mask(mol.coords, mol.coords[selected], ENVIRONMENT)
    near[selected] = True
    region = np.nonzero(near)[0].tolist()
    chosen = set(selected)
    sub = mol.subset(region)
    frozen = np.array([i not in chosen for i in region], dtype=bool)
    coords = _mmff_coordinates(sub, fixed=np.nonzero(frozen)[0].tolist())
    if coords is None:
        coords = minimize(ForceField(sub), sub.coords, max_steps=800, frozen=frozen)
    sub.coords = coords
    mol.set_positions(selected, sub.coords[[k for k, i in enumerate(region) if i in chosen]])


def clean_atoms(mol: Molecule, indices: Iterable[int], fold: bool = False) -> None:
    """Relax only the given atoms (normally one connected molecule) in place.

    Parameters
    ----------
    mol : Molecule
        Structure.
    indices : iterable of int
        Atoms to relax; the rest are ignored.
    fold : bool, optional
        Rebuild the coordinates from the bonds first (falling back to a soft
        relaxation), which lets a freshly closed ring reach a sensible shape.
    """
    idx = sorted(set(indices))
    if len(idx) < 2:
        return
    sub = mol.subset(idx)
    if fold:
        coords = reembed(sub)
        if coords is not None:
            mol.set_positions(idx, coords)
            return
        sub.coords = minimize(ForceField(sub, soft=True), sub.coords, max_steps=1500)
    clean_geometry(sub)
    mol.set_positions(idx, sub.coords)
