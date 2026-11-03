"""Bridge to Open Babel: conversion, SMILES, 3D building and force fields.

Open Babel numbers atoms from 1; everything returned here uses the 0-based
indices of :class:`~miles.core.molecule.Molecule`.
"""

import contextlib
import os
import sys
from collections.abc import Iterator, Sequence

import numpy as np
from openbabel import openbabel as ob

from miles.core.elements import BY_NUMBER, SYMBOLS
from miles.core.geometry import kabsch
from miles.core.molecule import RESIDUE_KEY, UNSET_CHARGE, Molecule
from miles.core.perception import infer_atom_types

ob.obErrorLog.SetOutputLevel(0)

RESIDUE_FORMATS = {"pdb", "ent", "pqr", "pdbqt", "mmcif", "cif", "mcif", "gro"}
"""Formats whose residue records are meaningful (others get a dummy residue)."""

FORCE_FIELDS = ("MMFF94", "UFF")
"""Force fields to try, in order of preference."""

ISOLATED_IONS = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Zn": 2,
                 "F": -1, "Cl": -1, "Br": -1, "I": -1}
"""Formal charge of an atom with no bonds, for the common monatomic ions."""


class ChemError(ValueError):
    """A structure or SMILES string Open Babel cannot handle."""


# -- conversion -----------------------------------------------------------------

def from_openbabel(obmol: "ob.OBMol", residues: bool = False, charges: bool = False) -> Molecule:
    """Convert an Open Babel molecule.

    Parameters
    ----------
    obmol : openbabel.OBMol
        Source molecule.
    residues : bool, optional
        Keep residue names, numbers, insertion codes, chains and atom names.
    charges : bool, optional
        Keep the formal charges, for sources that state them (SMILES, MDL
        molfiles); otherwise a charge of 0 cannot be told from no charge
        given, and the charges are left to be worked out from the bonds.

    Returns
    -------
    Molecule
        Atoms, coordinates and bonds (aromatic bonds get order 1.5), with
        isotopes and atom types inferred from the bonds.
    """
    numbers, isotopes, formal = [], [], []
    names, resnames, resseqs, icodes, chains, hetero = [], [], [], [], [], []
    for a in ob.OBMolAtomIter(obmol):
        z = a.GetAtomicNum()
        numbers.append(z if z in BY_NUMBER else 0)
        isotopes.append(a.GetIsotope())
        formal.append(a.GetFormalCharge() if charges else UNSET_CHARGE)
        res = a.GetResidue() if residues else None
        if res is not None:
            names.append(res.GetAtomID(a).strip())
            resnames.append(res.GetName().strip())
            resseqs.append(res.GetNum())
            icodes.append(res.GetInsertionCode().strip("\x00 "))
            chains.append(res.GetChain().strip())
            hetero.append(bool(res.IsHetAtom(a)))
    columns: dict[str, list] = {"number": numbers, "isotope": isotopes, "formal_charge": formal}
    if residues and names:
        columns |= {"name": names, "resname": resnames, "resseq": resseqs, "icode": icodes,
                    "chain": chains, "hetero": hetero}
    mol = Molecule.from_columns(openbabel_coordinates(obmol), **columns)
    bonds = [(b.GetBeginAtomIdx() - 1, b.GetEndAtomIdx() - 1,
              1.5 if b.IsAromatic() else float(b.GetBondOrder())) for b in ob.OBMolBondIter(obmol)]
    if bonds:
        table = np.array(bonds)
        mol.set_bond_arrays(table[:, :2].astype(np.intp), table[:, 2])
    infer_atom_types(mol, only_missing=False)
    return mol


def to_openbabel(mol: Molecule, title: str = "") -> "ob.OBMol":
    """Convert to an Open Babel molecule.

    Parameters
    ----------
    mol : Molecule
        Source structure; aromatic bonds (order 1.5) are kekulized.
    title : str, optional
        Molecule title.

    Returns
    -------
    openbabel.OBMol
        Molecule with the same atoms, coordinates, bonds and residues.
    """
    obmol = ob.OBMol()
    obmol.BeginModify()
    residues: dict[tuple[str, int, str, str], ob.OBResidue] = {}
    numbers, isotopes = mol.numbers.tolist(), mol.column("isotope").tolist()
    keys = list(zip(*(mol.column(c).tolist() for c in RESIDUE_KEY), strict=True))
    names, hetero = mol.column("name").tolist(), mol.column("hetero").tolist()
    for k, (x, y, z) in enumerate(mol.coords.tolist()):
        atom = obmol.NewAtom()
        atom.SetAtomicNum(numbers[k])
        atom.SetVector(x, y, z)
        if isotopes[k]:
            atom.SetIsotope(isotopes[k])
        key = keys[k]
        if key[3]:   # a residue name: the atom belongs to a residue
            res = residues.get(key)
            if res is None:
                chain, resseq, icode, resname = key
                res = residues[key] = obmol.NewResidue()
                res.SetName(resname)
                res.SetNum(int(resseq))
                res.SetChain(chain[:1] or " ")
                if icode:
                    res.SetInsertionCode(icode[:1])
            res.AddAtom(atom)
            res.SetAtomID(atom, names[k] or str(SYMBOLS[numbers[k]]))
            res.SetHetAtom(atom, bool(hetero[k]))
    aromatic = False
    for (i, j), order in zip(mol.bond_pairs.tolist(), mol.bond_orders.tolist(), strict=True):
        if order == 1.5:
            aromatic = True
            obmol.AddBond(i + 1, j + 1, 1)
            obmol.GetBond(i + 1, j + 1).SetAromatic()
            obmol.GetAtom(i + 1).SetAromatic()
            obmol.GetAtom(j + 1).SetAromatic()
        else:
            obmol.AddBond(i + 1, j + 1, max(1, round(order)))
    obmol.EndModify()
    if residues:
        obmol.SetChainsPerceived()   # the residues above are the real ones: writers must not renumber them
    if aromatic:
        obmol.SetAromaticPerceived()
        ob.OBKekulize(obmol)
    _set_formal_charges(mol, obmol)
    obmol.SetTitle(title)
    return obmol


def _bond_charge(number: int, valence: int, symbol: str) -> int:
    """Formal charge from the element and the sum of (kekulized) bond orders."""
    if valence == 0:
        return ISOLATED_IONS.get(symbol, 0)
    if number in (7, 15, 33):          # N, P, As
        if valence == 4:
            return 1
        if number == 7 and valence == 2:
            return -1
    elif number in (8, 16, 34):        # O, S, Se
        if valence == 1:
            return -1
        if valence == 3:
            return 1
    elif number == 5 and valence == 4:  # B
        return -1
    return 0


def _set_formal_charges(mol: Molecule, obmol: "ob.OBMol") -> None:
    """Give the Open Babel atoms their formal charges: as stated, or worked out from the bonds.

    A charge is worked out only where the bonds can tell: for isolated
    atoms (monatomic ions), and in molecules that have their hydrogens or
    were built in MILES. Heavy-atom-only molecules read from a file (most
    PDB proteins) are left neutral, since a hydroxyl oxygen without its
    hydrogen would otherwise look like an anion.
    """
    stated = mol.column("formal_charge")
    unset = stated == UNSET_CHARGE
    if not len(mol):
        return
    if unset.any():
        labels = mol.component_labels()
        complete = np.zeros(int(labels.max()) + 1, dtype=bool)
        complete[labels[mol.numbers == 1]] = True
        complete[labels[~mol.column("inferred")]] = True
        work_out = unset & (complete[labels] | (mol.degrees() == 0))
    else:
        work_out = unset
    symbols = mol.symbols
    numbers = mol.numbers.tolist()
    for k, atom in enumerate(ob.OBMolAtomIter(obmol)):
        if not unset[k]:
            atom.SetFormalCharge(int(stated[k]))
        elif work_out[k]:
            atom.SetFormalCharge(_bond_charge(numbers[k], atom.GetExplicitValence(), symbols[k]))


def openbabel_coordinates(obmol: "ob.OBMol") -> np.ndarray:
    """Return the coordinates of an Open Babel molecule (current conformer), shape ``(n, 3)``."""
    return np.array([(a.GetX(), a.GetY(), a.GetZ()) for a in ob.OBMolAtomIter(obmol)],
                    dtype=float).reshape(-1, 3)


def _align(new: np.ndarray, old: np.ndarray) -> np.ndarray:
    """Superpose ``new`` onto ``old`` (same atoms) by a rigid motion."""
    a, b = new - new.mean(axis=0), old - old.mean(axis=0)
    return a @ kabsch(a, b).T + old.mean(axis=0)


# -- SMILES ------------------------------------------------------------------------

def _read_smiles(text: str) -> "ob.OBMol | None":
    text = text.strip()
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    obmol = ob.OBMol()
    if not text or " " in text or not conv.ReadString(obmol, text) or not obmol.NumAtoms():
        return None
    return obmol


def is_smiles(text: str) -> bool:
    """Whether ``text`` is a valid SMILES string."""
    return _read_smiles(text) is not None


def from_smiles(text: str) -> tuple[Molecule, int | None]:
    """Build a 3D structure from SMILES.

    Parameters
    ----------
    text : str
        SMILES string. A ``*`` marks an attachment point, as in ``*C(=O)OC``.

    Returns
    -------
    mol : Molecule
        Structure with hydrogens, relaxed with a force field.
    attach : int or None
        Index of the hydrogen that replaced the first ``*``.

    Raises
    ------
    ChemError
        If the string is not valid SMILES.
    """
    obmol = _read_smiles(text)
    if obmol is None:
        raise ChemError(f"Not a valid SMILES string: {text}")
    attach = None
    for a in ob.OBMolAtomIter(obmol):
        if a.GetAtomicNum() == 0:
            a.SetAtomicNum(1)
            if attach is None:
                attach = a.GetIdx() - 1
    obmol.AddHydrogens()
    ob.OBBuilder().Build(obmol)
    _relax(obmol, steps=500)
    return from_openbabel(obmol, charges=True), attach


def to_smiles(mol: Molecule) -> str:
    """Write canonical SMILES, with stereochemistry taken from the 3D geometry.

    Raises
    ------
    ChemError
        If the molecule is empty.
    """
    if not len(mol):
        raise ChemError("No atoms")
    obmol = to_openbabel(mol)
    ob.StereoFrom3D(obmol)
    conv = ob.OBConversion()
    conv.SetOutFormat("can")
    return conv.WriteString(obmol).split()[0]


# -- force fields --------------------------------------------------------------------

def _forcefield(obmol: "ob.OBMol", constraints: "ob.OBFFConstraints | None" = None,
                names: Sequence[str] = FORCE_FIELDS) -> "ob.OBForceField | None":
    """First of the named force fields with parameters for every atom."""
    for name in names:
        ff = ob.OBForceField.FindForceField(name)
        if ff is not None and (ff.Setup(obmol) if constraints is None else ff.Setup(obmol, constraints)):
            return ff
    return None


def _relax(obmol: "ob.OBMol", steps: int = 500, constraints: "ob.OBFFConstraints | None" = None,
           names: Sequence[str] = FORCE_FIELDS) -> float | None:
    """Minimise ``obmol`` in place; return the energy, or None without a force field."""
    ff = _forcefield(obmol, constraints, names)
    if ff is None:
        return None
    ff.ConjugateGradients(steps, 1e-6)
    ff.GetCoordinates(obmol)
    return ff.Energy()


def minimize(mol: Molecule, fixed: Sequence[int] | None = None, steps: int = 800,
             force_fields: Sequence[str] = FORCE_FIELDS) -> np.ndarray | None:
    """Minimise with MMFF94, or UFF when MMFF94 has no parameters.

    Parameters
    ----------
    mol : Molecule
        Structure (not modified).
    fixed : sequence of int, optional
        Atoms that stay in place.
    steps : int, optional
        Conjugate-gradient steps.
    force_fields : sequence of str, optional
        Open Babel force fields to try, in order.

    Returns
    -------
    numpy.ndarray or None
        New coordinates, or None if no force field handles the structure.
    """
    if len(mol) < 2:
        return None
    obmol = to_openbabel(mol)
    constraints = None
    if fixed is not None and len(fixed):
        constraints = ob.OBFFConstraints()
        for i in fixed:
            constraints.AddAtomConstraint(int(i) + 1)
    if _relax(obmol, steps, constraints, force_fields) is None:
        return None
    return openbabel_coordinates(obmol)


def reembed(mol: Molecule) -> np.ndarray | None:
    """Build fresh 3D coordinates from the bonds alone.

    Used after closing a ring, when the old geometry is far from a sensible
    one. The result is relaxed and superposed onto the old coordinates.

    Parameters
    ----------
    mol : Molecule
        Structure (not modified).

    Returns
    -------
    numpy.ndarray or None
        New coordinates, or None if the structure could not be built.
    """
    if len(mol) < 3:
        return None
    obmol = to_openbabel(mol)
    for a in ob.OBMolAtomIter(obmol):
        a.SetVector(0.0, 0.0, 0.0)
    if not ob.OBBuilder().Build(obmol):
        return None
    _relax(obmol, 800)
    return _align(openbabel_coordinates(obmol), mol.coords)


@contextlib.contextmanager
def silence_openbabel() -> Iterator[None]:
    """Discard what Open Babel's C++ code prints to the standard output and error streams."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in (*saved, devnull):
            os.close(fd)


def conformers(mol: Molecule, count: int = 30, rmsd: float = 0.25) -> list[tuple[float, np.ndarray]]:
    """Search for low-energy conformers.

    Parameters
    ----------
    mol : Molecule
        Structure (not modified).
    count : int, optional
        Conformers generated by the genetic search.
    rmsd : float, optional
        Conformers closer than this (angstrom) to a lower-energy one are dropped.

    Returns
    -------
    list of tuple of (float, numpy.ndarray)
        Energy in kcal/mol and coordinates, lowest energy first, each
        superposed onto the input structure.

    Raises
    ------
    ChemError
        If no force field has parameters for the structure.
    """
    obmol = to_openbabel(mol)
    ff = _forcefield(obmol)
    if ff is None:
        raise ChemError("No force field parameters for this structure")
    if obmol.NumRotors() == 0:
        ff.ConjugateGradients(500, 1e-6)
        ff.GetCoordinates(obmol)
        return [(ff.Energy(), _align(openbabel_coordinates(obmol), mol.coords))]
    search = ob.OBConformerSearch()
    with silence_openbabel():
        search.Setup(obmol, count, 5, 5, 25)
        search.Search()
        search.GetConformers(obmol)
    results = []
    for k in range(obmol.NumConformers()):
        obmol.SetConformer(k)
        ff.SetCoordinates(obmol)
        ff.ConjugateGradients(300, 1e-6)
        ff.GetCoordinates(obmol)
        results.append((ff.Energy(), _align(openbabel_coordinates(obmol), mol.coords)))
    results.sort(key=lambda t: t[0])
    unique: list[tuple[float, np.ndarray]] = []
    for energy, xyz in results:
        if all(np.sqrt(((xyz - kept) ** 2).sum(axis=1).mean()) > rmsd for _, kept in unique):
            unique.append((energy, xyz))
    return unique
