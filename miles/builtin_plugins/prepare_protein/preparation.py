"""Preparing experimental protein structures: missing side chains, water and hydrogens.

Structures from the Protein Data Bank often lack atoms at the ends of
flexible side chains and never have hydrogens. :func:`prepare` completes the
side chains from the amino acid templates, optionally removes water, and
adds hydrogens for a chosen pH. The standard amino acids get their
hydrogens from their atom names: each ionisable group is protonated when the
pH is below its typical pKa (:data:`PKA`), histidine takes the tautomer that
suits its surroundings, and cysteines in disulfide bridges or bound to
metals keep no hydrogen. Everything else (ligands, nucleic acids, unusual
residues) is protonated with Open Babel's pH model.
"""

import functools
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from openbabel import openbabel as ob

from miles import (
    Atom,
    AtomLike,
    Molecule,
    add_hydrogens,
    amino_acid_codes,
    amino_acid_template,
    element,
    infer_atom_types,
    infer_geometry,
    openbabel_coordinates,
    residue_ids,
    to_openbabel,
    water_mask,
)
from miles.builtin_plugins.shared.mutation import BACKBONE, mutate_residue

DEFAULT_PH = 7.4

PKA = {"ASP": 3.9, "GLU": 4.1, "HIS": 6.0, "CYS": 8.3, "TYR": 10.1, "LYS": 10.5, "ARG": 12.5,
       "N-terminus": 8.0, "C-terminus": 3.1}
"""Typical pKa values of the ionisable groups of proteins."""

# Hydrogens on each heavy side-chain atom of the standard amino acids, with
# acids protonated, bases neutral and arginine charged; ionisable groups are
# adjusted for the pH by _amino_acid_hydrogens.
SIDE_CHAIN_HYDROGENS: dict[str, dict[str, int]] = {
    "ALA": {"CB": 3},
    "ARG": {"CB": 2, "CG": 2, "CD": 2, "NE": 1, "CZ": 0, "NH1": 2, "NH2": 2},
    "ASN": {"CB": 2, "CG": 0, "OD1": 0, "ND2": 2},
    "ASP": {"CB": 2, "CG": 0, "OD1": 0, "OD2": 1},
    "CYS": {"CB": 2, "SG": 1},
    "GLN": {"CB": 2, "CG": 2, "CD": 0, "OE1": 0, "NE2": 2},
    "GLU": {"CB": 2, "CG": 2, "CD": 0, "OE1": 0, "OE2": 1},
    "GLY": {},
    "HIS": {"CB": 2, "CG": 0, "ND1": 0, "CD2": 1, "CE1": 1, "NE2": 0},
    "ILE": {"CB": 1, "CG1": 2, "CG2": 3, "CD1": 3},
    "LEU": {"CB": 2, "CG": 1, "CD1": 3, "CD2": 3},
    "LYS": {"CB": 2, "CG": 2, "CD": 2, "CE": 2, "NZ": 2},
    "MET": {"CB": 2, "CG": 2, "SD": 0, "CE": 3},
    "MSE": {"CB": 2, "CG": 2, "SE": 0, "CE": 3},
    "PHE": {"CB": 2, "CG": 0, "CD1": 1, "CD2": 1, "CE1": 1, "CE2": 1, "CZ": 1},
    "PRO": {"CB": 2, "CG": 2, "CD": 2},
    "SER": {"CB": 2, "OG": 1},
    "THR": {"CB": 1, "OG1": 1, "CG2": 3},
    "TRP": {"CB": 2, "CG": 0, "CD1": 1, "CD2": 0, "NE1": 1, "CE2": 0, "CE3": 1, "CZ2": 1, "CZ3": 1,
            "CH2": 1},
    "TYR": {"CB": 2, "CG": 0, "CD1": 1, "CD2": 1, "CE1": 1, "CE2": 1, "CZ": 0, "OH": 1},
    "VAL": {"CB": 1, "CG1": 3, "CG2": 3},
}

# Residue names that fix a protonation state (Amber and CHARMM conventions):
# the standard residue and the state, which the pH does not change.
VARIANTS = {"HID": ("HIS", "ND1"), "HSD": ("HIS", "ND1"), "HIE": ("HIS", "NE2"), "HSE": ("HIS", "NE2"),
            "HIP": ("HIS", "both"), "HSP": ("HIS", "both"), "ASH": ("ASP", "protonated"),
            "GLH": ("GLU", "protonated"), "LYN": ("LYS", "neutral"), "CYX": ("CYS", "bridged"),
            "CYM": ("CYS", "deprotonated")}
STANDARD_RESIDUES = frozenset(SIDE_CHAIN_HYDROGENS) | frozenset(VARIANTS)
HBOND_REACH = 3.2   # angstrom; an N...O/N distance up to this counts as a possible hydrogen bond

ResidueKey = tuple[str, int, str]


@dataclass
class PreparationReport:
    """What :func:`prepare` changed.

    Attributes
    ----------
    side_chains_completed : int
        Residues whose side chains were rebuilt.
    hydrogens_added : int
        Hydrogens placed (all earlier hydrogens are replaced).
    waters_removed : int
        Water molecules deleted.
    """

    side_chains_completed: int = 0
    hydrogens_added: int = 0
    waters_removed: int = 0


@functools.cache
def _side_chain_names(code: str) -> frozenset[str]:
    """Names of the heavy side-chain atoms of an amino acid."""
    template = amino_acid_template(code).mol
    return frozenset(a.name for a in template.atoms if a.symbol != "H" and a.name not in BACKBONE)


def _residues(mol: Molecule, keys: Iterable[ResidueKey] | None = None) -> dict[ResidueKey, list[int]]:
    """Atoms of every residue (or of the given ones), keyed by chain, number and name."""
    wanted = None if keys is None else set(keys)
    found: dict[ResidueKey, list[int]] = {}
    for i, a in enumerate(mol.atoms):
        key = (a.chain, a.resseq, a.resname)
        if a.resname and (wanted is None or key in wanted):
            found.setdefault(key, []).append(i)
    return found


def _incomplete_keys(mol: Molecule) -> list[ResidueKey]:
    codes = set(amino_acid_codes())
    keys = []
    for key, atoms in _residues(mol).items():
        if key[2] not in codes:
            continue
        names = {mol.atoms[i].name for i in atoms if mol.atoms[i].symbol != "H"}
        if {"N", "CA", "C"} <= names and not _side_chain_names(key[2]) <= names:
            keys.append(key)
    return keys


def incomplete_residues(mol: Molecule) -> list[int]:
    """Amino acids with missing side-chain atoms.

    Parameters
    ----------
    mol : Molecule
        Structure with residue data and PDB atom names.

    Returns
    -------
    list of int
        The CA atom of each standard amino acid whose side chain lacks heavy
        atoms that its template has. Residues without N, CA and C are
        skipped, since their side chains cannot be placed.
    """
    residues = _residues(mol)
    return [next(i for i in residues[key] if mol.atoms[i].name == "CA") for key in _incomplete_keys(mol)]


def complete_side_chains(mol: Molecule) -> int:
    """Rebuild the incomplete side chains in place; return how many were rebuilt."""
    keys = _incomplete_keys(mol)
    for key in keys:
        atoms = _residues(mol, [key])[key]
        mutate_residue(mol, next(i for i in atoms if mol.atoms[i].name == "CA"), key[2])
    return len(keys)


def water_count(mol: Molecule) -> int:
    """Count the water molecules of a structure."""
    waters = water_mask(mol)
    return len(np.unique(residue_ids(mol)[waters])) if waters.any() else 0


def remove_water(mol: Molecule) -> int:
    """Delete every water molecule in place; return how many were deleted."""
    count = water_count(mol)
    mol.remove_atoms(np.flatnonzero(water_mask(mol)).tolist())
    return count


def _hydrogen_for(mol: Molecule, parent: int, name: str) -> Atom:
    """Make a hydrogen record in its parent's residue, so it is drawn and grouped like the parent."""
    p = mol.atoms[parent]
    return Atom("H", "linear", 1, 1, name=name, resname=p.resname, resseq=p.resseq, icode=p.icode,
                chain=p.chain, hetero=p.hetero)


def _hydrogen_names(parent_name: str, parent_symbol: str, count: int) -> list[str]:
    """PDB-style names for the hydrogens on one atom, e.g. HB2 and HB3, or H1 to H3 on a terminal N."""
    if parent_name == "N":
        return ["H"] if count == 1 else [f"H{k + 1}" for k in range(count)]
    base = "H" + parent_name[1:]
    if count == 1:
        return [base]
    first = 2 if count == 2 and parent_symbol == "C" else 1
    return [f"{base}{first + k}" for k in range(count)]


def _is_metal(atom: AtomLike) -> bool:
    return element(atom.symbol).is_metal


def _metal_bound(mol: Molecule, i: int) -> bool:
    return any(_is_metal(mol.atoms[j]) for j in mol.neighbors(i))


def _histidine_tautomer(mol: Molecule, nd1: int | None, ne2: int | None) -> str:
    """``"ND1"`` or ``"NE2"``: the ring nitrogen that carries the hydrogen of neutral histidine.

    A nitrogen bound to a metal cannot carry it. Otherwise the hydrogen goes
    where it can donate a hydrogen bond: to the nitrogen with an oxygen or
    nitrogen of another residue nearby when only one has one, else to NE2,
    the more common tautomer.
    """
    bound = {name: i is not None and _metal_bound(mol, i) for name, i in (("ND1", nd1), ("NE2", ne2))}
    if bound["NE2"] and not bound["ND1"]:
        return "ND1"
    if bound["ND1"] and not bound["NE2"]:
        return "NE2"
    near_nd1 = nd1 is not None and _polar_partner(mol, nd1)
    near_ne2 = ne2 is not None and _polar_partner(mol, ne2)
    return "ND1" if near_nd1 and not near_ne2 else "NE2"


def _polar_partner(mol: Molecule, i: int) -> bool:
    """Whether an N or O of another residue lies within :data:`HBOND_REACH` of atom ``i``."""
    a = mol.atoms[i]
    close = ((mol.coords - mol.coords[i]) ** 2).sum(axis=1) <= HBOND_REACH ** 2
    return any(mol.atoms[j].symbol in ("N", "O")
               and (mol.atoms[j].chain, mol.atoms[j].resseq) != (a.chain, a.resseq)
               for j in close.nonzero()[0])


def _amino_acid_hydrogens(mol: Molecule, ph: float) -> dict[int, int]:
    """Count the hydrogens of every heavy atom of the standard amino acids at a pH.

    Atoms of other residues, and atoms with names a residue does not have,
    are left out, so they can be protonated by other means.
    """
    counts: dict[int, int] = {}
    for (_, _, resname), atoms in _residues(mol).items():
        name, state = VARIANTS.get(resname, (resname, ""))
        table = SIDE_CHAIN_HYDROGENS.get(name)
        if table is None:
            continue
        index = {mol.atoms[i].name: i for i in atoms if mol.atoms[i].symbol != "H"}
        here: dict[int, int] = {}
        for atom_name, i in index.items():
            if atom_name in table:
                here[i] = table[atom_name]
            elif atom_name == "CA":
                here[i] = 2 if name == "GLY" else 1
            elif atom_name in ("C", "O", "OXT"):
                here[i] = 0
            elif atom_name == "N":
                here[i] = 0 if name == "PRO" else 1
        if "N" in index and not any(mol.atoms[j].name == "C" for j in mol.neighbors(index["N"])):
            here[index["N"]] += 2 if ph < PKA["N-terminus"] else 1   # free amine at the chain start
        if "OXT" in index and ph < PKA["C-terminus"]:
            here[index["OXT"]] = 1
        if name in ("ASP", "GLU"):
            acid = index.get("OD2" if name == "ASP" else "OE2")
            if acid is not None and state != "protonated" and (state or ph >= PKA[name]):
                here[acid] = 0
        elif name == "LYS" and "NZ" in index and (state != "neutral" and ph < PKA["LYS"]):
            here[index["NZ"]] = 3
        elif name == "ARG" and "NH2" in index and ph >= PKA["ARG"]:
            here[index["NH2"]] = 1
        elif name == "TYR" and "OH" in index and ph >= PKA["TYR"]:
            here[index["OH"]] = 0
        elif name == "CYS" and "SG" in index:
            sg = index["SG"]
            bridged = _metal_bound(mol, sg) or any(mol.atoms[j].symbol == "S" for j in mol.neighbors(sg))
            if bridged or state in ("bridged", "deprotonated") or (not state and ph >= PKA["CYS"]):
                here[sg] = 0
        elif name == "HIS":
            nd1, ne2 = index.get("ND1"), index.get("NE2")
            both = state == "both" or (not state and ph < PKA["HIS"])
            chosen = state if state in ("ND1", "NE2") else _histidine_tautomer(mol, nd1, ne2)
            for ring_name, ring_atom in (("ND1", nd1), ("NE2", ne2)):
                if ring_atom is not None and not _metal_bound(mol, ring_atom):
                    here[ring_atom] = 1 if both or ring_name == chosen else 0
        for i in list(here):   # each bond to a metal or a ligand (a covalent inhibitor) replaces a hydrogen
            unusual = sum(1 for j in mol.neighbors(i)
                          if _is_metal(mol.atoms[j]) or mol.atoms[j].resname not in STANDARD_RESIDUES)
            here[i] = max(0, here[i] - unusual)
        counts.update(here)
    return counts


def _place_counted_hydrogens(mol: Molecule, counts: dict[int, int]) -> int:
    """Add exactly ``counts[i]`` hydrogens to each atom, with PDB-style names; return how many."""
    added = 0
    for i, n in counts.items():
        if n <= 0:
            continue
        a = mol.atoms[i]
        a.nsites = mol.degree(i) + n
        a.valence = mol.order_sum(i) + n
        a.geom = infer_geometry(mol, i, a.nsites)
        a.inferred = False
        first = len(mol)
        made = add_hydrogens(mol, [i])
        names = _hydrogen_names(a.name, a.symbol, made) if a.name else [""] * made
        for k, name in zip(range(first, first + made), names, strict=True):
            mol.atoms[k].name = name
        added += made
    return added


def add_hydrogens_at_ph(mol: Molecule, ph: float = DEFAULT_PH) -> int:
    """Replace all hydrogens by the hydrogens a structure carries at a given pH.

    Every heavy atom keeps its record (residue data, representation,
    visibility, secondary structure) and its order. The standard amino acids
    are protonated from their atom names and the pKa values in :data:`PKA`:
    at pH 7.4, for example, Lys and Arg side chains and the N-terminus are
    protonated and Asp, Glu and the C-terminus deprotonated. Everything else
    gets its usual number of hydrogens from Open Babel, adjusted by its pH
    model. New hydrogens join their parent's residue with PDB-style names,
    and the bonding pattern of every heavy atom is updated so that later
    commands (such as Add Hydrogens) take the new state as complete.

    Parameters
    ----------
    mol : Molecule
        Structure, changed in place.
    ph : float, optional
        pH for the protonation states.

    Returns
    -------
    int
        Number of hydrogens placed. The hydrogens present before are removed
        first, so this is the total number of hydrogens afterwards.
    """
    mol.remove_atoms([i for i, a in enumerate(mol.atoms) if a.symbol == "H"])
    heavy = len(mol)
    if not heavy:
        return 0
    counts = _amino_acid_hydrogens(mol, ph)
    by_parent: dict[int, list[int]] = {}
    if len(counts) < heavy:
        obmol = to_openbabel(mol)
        for atom in ob.OBMolAtomIter(obmol):
            ob.OBAtomAssignTypicalImplicitHydrogens(atom)
        obmol.AddHydrogens(False, True, ph)
        xyz = openbabel_coordinates(obmol)
        for k in range(heavy, obmol.NumAtoms()):
            neighbor = next(iter(ob.OBAtomAtomIter(obmol.GetAtom(k + 1))), None)
            if neighbor is not None and neighbor.GetIdx() - 1 not in counts:
                by_parent.setdefault(neighbor.GetIdx() - 1, []).append(k)
    added = 0
    for parent, hydrogens in by_parent.items():
        p = mol.atoms[parent]
        names = _hydrogen_names(p.name, p.symbol, len(hydrogens)) if p.name else [""] * len(hydrogens)
        for k, name in zip(hydrogens, names, strict=True):
            mol.add_bond(parent, mol.add_atom(_hydrogen_for(mol, parent, name), xyz[k]), 1.0)
            added += 1
    added += _place_counted_hydrogens(mol, counts)
    infer_atom_types(mol, only_missing=False)
    mol.set_column("inferred", np.zeros(len(mol), dtype=bool))
    return added


def prepare(mol: Molecule, complete: bool = True, hydrogens: bool = True, ph: float = DEFAULT_PH,
            waters: bool = False) -> PreparationReport:
    """Prepare a protein structure in place.

    Parameters
    ----------
    mol : Molecule
        Structure with residue data.
    complete : bool, optional
        Rebuild side chains with missing atoms.
    hydrogens : bool, optional
        Replace the hydrogens by those for pH ``ph``.
    ph : float, optional
        pH for the protonation states.
    waters : bool, optional
        Delete water molecules (before anything else).

    Returns
    -------
    PreparationReport
        What was changed.
    """
    report = PreparationReport()
    if waters:
        report.waters_removed = remove_water(mol)
    if complete:
        report.side_chains_completed = complete_side_chains(mol)
    if hydrogens:
        report.hydrogens_added = add_hydrogens_at_ph(mol, ph)
    return report
