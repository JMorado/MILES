"""Hydrogen bonds within a structure: the one definition the built-in plugins share.

View > Hydrogen Bonds draws what :class:`HydrogenBondFinder` finds, and the
binding pocket counts the hydrogen bonds it finds between a ligand and its
pocket, with the same criteria and the same cutoffs.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from miles import AMINO_ACIDS, NUCLEOTIDES, Molecule, is_water, pairs_within, residue_ids

HBOND_MAX = 3.5          # donor-acceptor distance, angstrom
HBOND_MIN = 2.4          # anything closer is a clash or a bond
HBOND_ANGLE = 120.0      # minimum D-H...A angle, degrees

DONORS = {"N", "O", "S"}
ACCEPTORS = {"N", "O"}


class _HydrogenModel:
    """Whether hydrogens were modelled around an atom: any hydrogen in its residue (or its molecule)."""

    def __init__(self, mol: Molecule) -> None:
        ids = residue_ids(mol)
        hydrogens = (mol.numbers == 1).astype(float)
        counts = np.bincount(ids, weights=hydrogens, minlength=int(ids.max(initial=0)) + 1)
        self._modelled = counts[ids] > 0 if len(ids) else np.zeros(0, dtype=bool)

    def __call__(self, i: int) -> bool:
        return bool(self._modelled[i])

HBOND_SEPARATION = 3     # donor and acceptor must be more than this many bonds apart
HBOND_HEAVY_ANGLE = 90.0  # minimum X-D...A and D...A-Y angles with the atoms' heavy neighbours, degrees
NH_LENGTH = 1.01          # N-H bond length for hydrogens placed on planar NH groups, angstrom
VALENCE = {"N": 3, "O": 2, "S": 2}

# Roles of amino acid atoms, used when a residue has no hydrogens to show its donors.
_BACKBONE_DONORS = {"N"}
_BACKBONE_ACCEPTORS = {"O", "OXT"}
_SIDE_CHAIN_DONORS = {
    "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"}, "HIS": {"ND1", "NE2"}, "LYS": {"NZ"},
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "TRP": {"NE1"}, "CYS": {"SG"},
    "ASH": {"OD2"}, "GLH": {"OE2"}, "LYN": {"NZ"},
}
_SIDE_CHAIN_ACCEPTORS = {
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"}, "GLN": {"OE1"}, "HIS": {"ND1", "NE2"},
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "ASH": {"OD1", "OD2"}, "GLH": {"OE1", "OE2"},
    "LYN": {"NZ"},
}
_RESIDUE_ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
                    "CYX": "CYS", "CYM": "CYS", "MSE": "MET"}
# Roles of nucleobase atoms (ring NH groups cannot be told from pyridine-type nitrogens by their bonds).
_BASE_DONORS = {"A": {"N6"}, "G": {"N1", "N2"}, "C": {"N4"}, "T": {"N3"}, "U": {"N3"}, "I": {"N1"}}
_BASE_ACCEPTORS = {"A": {"N1", "N3", "N7"}, "G": {"O6", "N3", "N7"}, "C": {"O2", "N3"}, "T": {"O2", "O4"},
                   "U": {"O2", "O4"}, "I": {"O6", "N3", "N7"}}


@dataclass(frozen=True, slots=True)
class HydrogenBond:
    """A hydrogen bond within a structure.

    Attributes
    ----------
    donor : int
        Atom carrying the hydrogen (N, O or S).
    hydrogen : int or None
        The hydrogen, or None when the donor's hydrogens are not in the model.
    acceptor : int
        Atom accepting the hydrogen bond (N or O).
    distance : float
        Donor-acceptor distance, in angstrom.
    """

    donor: int
    hydrogen: int | None
    acceptor: int
    distance: float


def _implicit_hydrogens(mol: Molecule, i: int) -> int:
    """Hydrogens an N, O or S atom would carry to fill its usual valence."""
    return max(0, round(VALENCE.get(mol.atoms[i].symbol, 0) - mol.order_sum(i)))


def _conjugated(mol: Molecule, i: int) -> bool:
    """Whether atom ``i`` has a double, triple or aromatic bond."""
    return any(mol.bond_order(i, j) > 1.0 for j in mol.neighbors(i))


def _accepts_by_bonds(mol: Molecule, i: int, implicit: int) -> bool:
    """Whether an N or O atom has a lone pair free to accept a hydrogen bond, judged from its bonds.

    Oxygen accepts unless it has three substituents. Nitrogen accepts with
    one or two substituents (pyridine, imines, nitriles) and as an amine;
    with four (ammonium), or three when it is conjugated with its neighbours
    (amides, anilines, pyrrole-type ring nitrogens, nitro groups), its lone
    pair is taken.
    """
    symbol = mol.atoms[i].symbol
    substituents = mol.degree(i) + implicit
    if symbol == "O":
        return substituents <= 2
    if symbol != "N" or substituents >= 4:
        return False
    if substituents <= 2:
        return True
    return not (_conjugated(mol, i) or any(_conjugated(mol, j) for j in mol.neighbors(i)
                                           if mol.atoms[j].symbol != "H"))


def _base_roles(mol: Molecule, i: int) -> tuple[bool, bool] | None:
    """Whether a nucleobase atom donates (when hydrogens are missing) and accepts; None for other atoms."""
    a = mol.atoms[i]
    base = a.resname[-1] if a.resname in NUCLEOTIDES else ""
    if a.name not in _BASE_DONORS.get(base, set()) | _BASE_ACCEPTORS.get(base, set()):
        return None
    return a.name in _BASE_DONORS[base], a.name in _BASE_ACCEPTORS[base]


def _amino_acid_roles(mol: Molecule, i: int) -> tuple[bool, bool]:
    """Whether an atom of a standard residue donates (when hydrogens are missing) and accepts."""
    a = mol.atoms[i]
    resname = _RESIDUE_ALIASES.get(a.resname, a.resname)
    donor = (a.name in _BACKBONE_DONORS and resname != "PRO") or a.name in _SIDE_CHAIN_DONORS.get(resname, ())
    acceptor = a.name in _BACKBONE_ACCEPTORS or a.name in _SIDE_CHAIN_ACCEPTORS.get(resname, ())
    return donor, acceptor


class HydrogenBondFinder:
    """Finds the hydrogen bonds of one structure, for any positions of its atoms.

    Everything that depends only on the bonds is worked out once: which atoms
    can donate (with their hydrogens) and accept, and which pairs are too
    close in the bond graph. :meth:`find` then only searches distances and
    angles, so the hydrogen bonds of every frame of a trajectory are found
    quickly.

    Criteria: donor and acceptor (N, O; S
    only donates) 2.4 to 3.5 angstrom apart, more than three bonds apart (so
    neighbours, 1-3 and 1-4 pairs and atoms of the same five- or six-membered
    ring are left out), a D-H...A angle above 120 degrees, and angles above 90
    degrees between the hydrogen bond and the bonds of donor and acceptor to
    their heavy neighbours (so nothing is accepted from behind).

    When a residue or molecule has no hydrogens at all, donors are taken from
    atom names for amino acids and nucleobases and from unfilled valences elsewhere
    (a hydroxyl or water oxygen, an NH or NH2 nitrogen). A nitrogen with two
    heavy neighbours and one missing hydrogen (peptide NH, ring NH) is planar,
    so its hydrogen is placed on the bisector, as DSSP does, and the
    D-H...A angle is checked as usual; for other donors, whose hydrogens can
    rotate, the heavy-atom angles decide. Acceptors are oxygens and nitrogens
    with a free lone pair: amino acid and nucleobase atoms by name, other
    nitrogens by their bonds (amines and pyridine-type nitrogens accept;
    amides, anilines, pyrrole-type nitrogens and ammonium do not). Aromatic
    NH groups of other molecules without hydrogens cannot be told from
    pyridine-type nitrogens and count as acceptors only.

    Parameters
    ----------
    mol : Molecule
        Structure; its bonds and atom names are read once.
    atoms : sequence of int, optional
        Atoms that may take part (default: all). Hydrogens are found through
        their donors, so they need not be included.
    """

    def __init__(self, mol: Molecule, atoms: Sequence[int] | None = None) -> None:
        allowed = range(len(mol)) if atoms is None else sorted({int(i) for i in atoms})
        modelled = _HydrogenModel(mol)
        polar, donates, accepts, hydrogens, heavy = [], [], [], [], []
        for i in allowed:
            a = mol.atoms[i]
            if a.symbol not in DONORS:
                continue
            has_model = modelled(i)
            implicit = 0 if has_model else _implicit_hydrogens(mol, i)
            base = _base_roles(mol, i)
            if a.resname in AMINO_ACIDS or base is not None:
                named_donor, acceptor = base if base is not None else _amino_acid_roles(mol, i)
                bare_donor = named_donor and not has_model
            elif is_water(a):
                bare_donor, acceptor = not has_model, True
            else:
                bare_donor = implicit > 0
                acceptor = a.symbol in ACCEPTORS and _accepts_by_bonds(mol, i, implicit)
            hs = mol.hydrogens_on(i) if has_model else []
            if not (bare_donor or hs or acceptor):
                continue
            polar.append(i)
            donates.append(bare_donor)
            accepts.append(acceptor)
            hydrogens.append(hs[:4])
            heavy.append([j for j in mol.neighbors(i) if mol.atoms[j].symbol != "H"][:4])
        n = len(polar)
        self.polar = np.array(polar, dtype=np.int64)
        self._bare_donor = np.array(donates, dtype=bool)
        self._acceptor = np.array(accepts, dtype=bool)
        self._hydrogens = np.full((n, 4), -1, dtype=np.int64)
        self._heavy = np.full((n, 4), -1, dtype=np.int64)
        for k in range(n):
            self._hydrogens[k, :len(hydrogens[k])] = hydrogens[k]
            self._heavy[k, :len(heavy[k])] = heavy[k]
        # planar NH donors without modelled hydrogens get one on the bisector of their two bonds
        self._placed = self._bare_donor & np.array(
            [mol.atoms[i].symbol == "N" and len(heavy[k]) == 2 for k, i in enumerate(polar)], dtype=bool)
        self._excluded = self._close_in_graph(mol)

    def _close_in_graph(self, mol: Molecule) -> np.ndarray:
        """Sorted keys ``a * n + b`` (local indices, a < b) of polar pairs too close in the bond graph."""
        local = {int(i): k for k, i in enumerate(self.polar.tolist())}
        n = len(local)
        keys: list[int] = []
        for i, k in local.items():
            seen = {i}
            frontier = [i]
            for _ in range(HBOND_SEPARATION):
                frontier = [j for v in frontier for j in mol.neighbors(v) if j not in seen]
                seen.update(frontier)
            keys.extend(min(k, local[j]) * n + max(k, local[j]) for j in seen if j in local and j != i)
        return np.unique(np.array(keys, dtype=np.int64))

    def find(self, coords: np.ndarray, max_distance: float = HBOND_MAX,
             min_angle: float = HBOND_ANGLE) -> list[HydrogenBond]:
        """Hydrogen bonds for one set of atom positions.

        Parameters
        ----------
        coords : numpy.ndarray
            Positions of all atoms of the structure, shape ``(n_atoms, 3)``.
        max_distance : float, optional
            Largest donor-acceptor distance, angstrom.
        min_angle : float, optional
            Smallest D-H...A angle, degrees.

        Returns
        -------
        list of HydrogenBond
            Each donor-acceptor pair once, in the direction with the straighter
            D-H...A angle when both atoms could donate.
        """
        if len(self.polar) < 2:
            return []
        x = np.asarray(coords, dtype=float)
        pairs = pairs_within(x[self.polar], max_distance)
        if not len(pairs):
            return []
        a, b = pairs[:, 0], pairs[:, 1]
        d = np.linalg.norm(x[self.polar[a]] - x[self.polar[b]], axis=1)
        keep = (d > HBOND_MIN) & ~np.isin(a * len(self.polar) + b, self._excluded)
        a, b, d = a[keep], b[keep], d[keep]
        cos_ab, h_ab = self._best_hydrogen(x, a, b)
        cos_ba, h_ba = self._best_hydrogen(x, b, a)
        limit = math.cos(math.radians(min_angle))
        open_ab = self._open(x, a, b) & self._open(x, b, a) & self._acceptor[b]   # a donates to b
        open_ba = self._open(x, b, a) & self._open(x, a, b) & self._acceptor[a]
        explicit_ab = (cos_ab < limit) & open_ab
        explicit_ba = (cos_ba < limit) & open_ba
        bare_ab = self._bare_donor[a] & ~self._placed[a] & open_ab
        bare_ba = self._bare_donor[b] & ~self._placed[b] & open_ba
        found = []
        for k in np.nonzero(explicit_ab | explicit_ba | bare_ab | bare_ba)[0].tolist():
            if explicit_ab[k] or explicit_ba[k]:
                forward = explicit_ab[k] and (not explicit_ba[k] or cos_ab[k] <= cos_ba[k])
                donor, acceptor = (a[k], b[k]) if forward else (b[k], a[k])
                h = int(h_ab[k] if forward else h_ba[k])
                hydrogen: int | None = h if h >= 0 else None   # -1: placed on a planar NH
            else:
                donor, acceptor = (a[k], b[k]) if bare_ab[k] else (b[k], a[k])
                hydrogen = None
            found.append(HydrogenBond(int(self.polar[donor]), hydrogen, int(self.polar[acceptor]),
                                      float(d[k])))
        return found

    def _open(self, x: np.ndarray, atoms: np.ndarray, partners: np.ndarray) -> np.ndarray:
        """Whether each partner lies more than 90 degrees from every bond of the atom to a heavy neighbour."""
        heavy = self._heavy[atoms]                                    # (k, 4), -1 for none
        valid = heavy >= 0
        centre = x[self.polar[atoms]][:, None, :]
        u = x[np.where(valid, heavy, 0)] - centre
        v = x[self.polar[partners]][:, None, :] - centre
        norms = np.linalg.norm(u, axis=2) * np.linalg.norm(v, axis=2)
        cos = (u * v).sum(axis=2) / np.where(norms > 0, norms, 1.0)
        limit = math.cos(math.radians(HBOND_HEAVY_ANGLE))
        return np.all(~valid | (cos < limit + 1e-9), axis=1)

    def _best_hydrogen(self, x: np.ndarray, donors: np.ndarray, acceptors: np.ndarray
                       ) -> tuple[np.ndarray, np.ndarray]:
        """Cosine of the straightest D-H...A angle over each donor's hydrogens, and that hydrogen.

        Placed hydrogens of planar NH groups count with index -1. Pairs whose
        donor has no hydrogens get a cosine of 2 (never a hydrogen bond).
        """
        h = self._hydrogens[donors]                                   # (k, 4), -1 for none
        valid = h >= 0
        hx = x[np.where(valid, h, 0)]                                 # (k, 4, 3)
        placed = self._placed[donors]
        if placed.any():
            d = x[self.polar[donors[placed]]]
            ends = x[self._heavy[donors[placed], :2]]                 # (m, 2, 3)
            out = d[:, None, :] - ends
            out /= np.linalg.norm(out, axis=2, keepdims=True)
            bisector = out.sum(axis=1)
            bisector /= np.maximum(np.linalg.norm(bisector, axis=1, keepdims=True), 1e-9)
            hx[placed, 0] = d + NH_LENGTH * bisector
            valid = valid.copy()
            valid[placed, 0] = True
        u = x[self.polar[donors]][:, None, :] - hx
        v = x[self.polar[acceptors]][:, None, :] - hx
        norms = np.linalg.norm(u, axis=2) * np.linalg.norm(v, axis=2)
        cos = np.where(valid & (norms > 0), (u * v).sum(axis=2) / np.where(norms > 0, norms, 1.0), 2.0)
        best = cos.argmin(axis=1)
        rows = np.arange(len(donors))
        return cos[rows, best], h[rows, best]


def hydrogen_bonds(mol: Molecule, atoms: Sequence[int] | None = None) -> list[HydrogenBond]:
    """All hydrogen bonds within a structure, with the criteria of :class:`HydrogenBondFinder`.

    Parameters
    ----------
    mol : Molecule
        Structure.
    atoms : sequence of int, optional
        Atoms that may take part (default: all).

    Returns
    -------
    list of HydrogenBond
        Each bond once.
    """
    return HydrogenBondFinder(mol, atoms).find(mol.coords)
