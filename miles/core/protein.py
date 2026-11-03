"""Residue classes, groups of atoms, default representations and secondary structure."""

from enum import StrEnum

import numpy as np

from miles.core.molecule import RESIDUE_KEY, AtomLike, Molecule, Representation

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "MSE", "HID", "HIE", "HIP",
    "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN",
}
"""Residue names of amino acids, including common protonation variants."""

NUCLEOTIDES = {"A", "C", "G", "U", "T", "I", "DA", "DC", "DG", "DT", "DU", "DI"}
"""Residue names of RNA and DNA nucleotides."""

WATERS = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"}
"""Residue names used for water."""

POLYMER_RESIDUES = AMINO_ACIDS | NUCLEOTIDES
"""Residue names of amino acids and nucleotides."""

ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C", "CYM": "C",
    "ASH": "D", "GLH": "E", "LYN": "K",
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "T", "I": "I",
    "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "DI": "I",
}
"""One-letter codes of amino acid and nucleotide residues."""


def is_polymer(atom: AtomLike) -> bool:
    """Whether an atom belongs to an amino acid or nucleotide residue."""
    return atom.resname in AMINO_ACIDS or atom.resname in NUCLEOTIDES


def is_water(atom: AtomLike) -> bool:
    """Whether an atom belongs to a water molecule."""
    return atom.resname in WATERS


def polymer_mask(mol: Molecule) -> np.ndarray:
    """Boolean mask of the atoms in amino acid and nucleotide residues."""
    return mol.mask("resname", POLYMER_RESIDUES)


def water_mask(mol: Molecule) -> np.ndarray:
    """Boolean mask of the atoms in water molecules."""
    return mol.mask("resname", WATERS)


def has_polymer(mol: Molecule) -> bool:
    """Whether a structure contains a protein or nucleic acid chain."""
    names = mol.column("name")
    alpha, phosphorus = names == "CA", names == "P"
    if not (alpha.any() or phosphorus.any()):
        return False
    return bool(((alpha & mol.mask("resname", AMINO_ACIDS))
                 | (phosphorus & mol.mask("resname", NUCLEOTIDES))).any())


def default_representations(mol: Molecule) -> np.ndarray | None:
    """Representations a structure opens with: the polymer as a cartoon, water hidden, the rest as usual.

    Chains without secondary structure records get an estimate from their
    CA trace, written to the molecule's ``ss`` column.

    Parameters
    ----------
    mol : Molecule
        Structure.

    Returns
    -------
    numpy.ndarray or None
        Representation of each atom (``""`` for the view's style), or None
        when the structure has no polymer.
    """
    if not has_polymer(mol):
        return None
    polymer = polymer_mask(mol)
    reps = np.full(len(mol), str(Representation.DEFAULT), dtype=object)
    reps[polymer] = str(Representation.CARTOON)
    reps[water_mask(mol) & ~polymer] = str(Representation.HIDDEN)
    estimate_secondary_structure(mol)
    return reps


class Group(StrEnum):
    """Groups of atoms that can be selected together."""

    RESIDUE = "residue"
    CHAIN = "chain"
    POLYMER = "polymer"   # protein and nucleic acids
    LIGANDS = "ligands"
    WATER = "water"
    ALL = "all"


def residues(mol: Molecule) -> list[list[int]]:
    """Atoms of every amino acid and nucleotide residue, in file order.

    Parameters
    ----------
    mol : Molecule
        Structure with residue data.

    Returns
    -------
    list of list of int
        One list of atom indices per polymer residue.
    """
    idx = np.flatnonzero(polymer_mask(mol))
    keys = zip(*(mol.column(c)[idx].tolist() for c in RESIDUE_KEY), strict=True)
    found: dict[tuple[str, int, str, str], list[int]] = {}
    for i, key in zip(idx.tolist(), keys, strict=True):
        found.setdefault(key, []).append(i)
    return list(found.values())


def group_atoms(mol: Molecule, group: Group, around: int | None = None) -> list[int]:
    """Return the indices of the atoms in a group.

    Parameters
    ----------
    mol : Molecule
        Structure.
    group : Group
        Which group.
    around : int, optional
        Atom whose residue or chain is meant (needed for ``RESIDUE`` and ``CHAIN``).

    Returns
    -------
    list of int
        Atom indices; everything for ``ALL`` or when ``around`` is missing.
    """
    if group == Group.RESIDUE and around is not None:
        return mol.residue_atoms(around)
    if group == Group.CHAIN and around is not None:
        return mol.chain_atoms(around)
    if group == Group.POLYMER:
        return np.flatnonzero(polymer_mask(mol)).tolist()
    if group == Group.LIGANDS:
        from miles.core.selection import find  # the selection language defines "ligand" once

        return find(mol, "ligand").tolist()
    if group == Group.WATER:
        return np.flatnonzero(water_mask(mol)).tolist()
    return list(range(len(mol)))


# -- secondary structure ---------------------------------------------------------------

def _alpha_carbons(mol: Molecule) -> dict[str, list[tuple[tuple[int, str], int]]]:
    """Residue number with insertion code, and CA index, of each residue of each protein chain, in order."""
    idx = np.flatnonzero((mol.column("name") == "CA") & mol.mask("resname", AMINO_ACIDS))
    chains: dict[str, list[tuple[tuple[int, str], int]]] = {}
    seen = set()
    keys = zip(*(mol.column(c)[idx].tolist() for c in RESIDUE_KEY), strict=True)
    for i, key in zip(idx.tolist(), keys, strict=True):
        if key not in seen:
            seen.add(key)
            chain, resseq, icode, _ = key
            chains.setdefault(chain, []).append(((resseq, icode), i))
    return chains


def _drop_short_runs(flags: np.ndarray, minimum: int) -> np.ndarray:
    """Clear runs of True values shorter than ``minimum``."""
    out = flags.copy()
    start = None
    for i, flag in enumerate([*flags.tolist(), False]):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            if i - start < minimum:
                out[start:i] = False
            start = None
    return out


def estimate_secondary_structure(mol: Molecule) -> None:
    """Assign helices and strands from CA distances.

    Uses distance criteria between residues ``i`` and ``i + 2``, ``i + 3``
    and ``i + 4`` in the spirit of P-SEA. Chains that already carry
    secondary structure (from the file) are left alone.

    Parameters
    ----------
    mol : Molecule
        Structure whose ``Atom.ss`` fields are set in place.
    """
    ss = mol.column("ss")
    for chain, residues in _alpha_carbons(mol).items():
        idx = [ca for _, ca in residues]
        if len(idx) < 5 or any(ss[idx]):
            continue
        ca = mol.coords[idx]
        n = len(ca)
        helix = np.zeros(n, bool)
        strand = np.zeros(n, bool)
        for i in range(n - 4):
            d2, d3, d4 = (float(np.linalg.norm(ca[i] - ca[i + k])) for k in (2, 3, 4))
            if abs(d2 - 5.5) < 0.5 and abs(d3 - 5.3) < 0.5 and abs(d4 - 6.4) < 0.6:
                helix[i:i + 5] = True
            if abs(d2 - 6.7) < 0.6 and abs(d3 - 9.9) < 0.9 and abs(d4 - 12.4) < 1.1:
                strand[i:i + 3] = True
        helix = _drop_short_runs(helix, 5)
        strand = _drop_short_runs(strand & ~helix, 3)
        codes = np.where(helix, "H", np.where(strand, "E", ""))
        lookup = {number: str(code) for (number, _), code in zip(residues, codes, strict=True)}
        members = np.flatnonzero((mol.column("chain") == chain) & mol.mask("resname", AMINO_ACIDS))
        old = mol.column("ss")[members].tolist()
        resseq, icode = mol.column("resseq")[members].tolist(), mol.column("icode")[members].tolist()
        numbers = zip(resseq, icode, strict=True)
        new = [lookup.get(r, s) for r, s in zip(numbers, old, strict=True)]
        mol.set_column("ss", new, members)
