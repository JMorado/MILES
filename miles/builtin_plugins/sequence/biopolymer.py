"""Building peptides and nucleic acids from their sequences.

Peptides are built here: a poly-glycine backbone is laid out from internal
coordinates with the backbone torsions of an alpha helix or a beta strand,
and every residue is then given its side chain with
:func:`~miles.builtin_plugins.shared.mutation.mutate_residue`, which picks the least
crowded rotamer. DNA (single or double stranded) and RNA come from Open
Babel's FASTA reader, which builds B-form helices.

Open Babel can build peptides too, but its helices put carbonyl oxygens
within 1 angstrom of the next residue's nitrogen, so they are not used. Its
RNA has O3'-P bonds of under 1 angstrom and 2'-hydroxyl groups inside the
next base, so RNA is relaxed with MMFF94 a few nucleotides at a time.
"""

import math
import re
from enum import StrEnum

import numpy as np
from openbabel import openbabel as ob

from miles import (
    Atom,
    Molecule,
    clean_selection,
    from_openbabel,
    infer_atom_types,
    minimize,
    pairs_within,
    place_atom,
    unit_vector,
)
from miles.builtin_plugins.shared.mutation import mutate_residue

MAX_RESIDUES = 2000
"""Longest sequence accepted."""

AMINO_ACID_CODES = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
    "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
    "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
"""Three-letter code of each one-letter amino acid code."""

# Backbone geometry (Engh and Huber): bond lengths in angstrom, angles in degrees.
N_CA, CA_C, C_N, C_O, C_OXT = 1.458, 1.525, 1.329, 1.231, 1.31
N_H, CA_H, O_H = 1.01, 1.09, 0.97
ANGLE_N_CA_C, ANGLE_CA_C_N, ANGLE_C_N_CA, ANGLE_CA_C_O = 111.2, 116.2, 121.7, 120.5
OMEGA = 180.0
TETRAHEDRAL_HALF = math.radians(109.47 / 2)

HELIX_PHI_PSI = (-57.0, -47.0)
STRAND_PHI_PSI = (-120.0, 130.0)

BACKBONE_NAMES = frozenset({"N", "CA", "C", "O", "H", "H1", "H2", "HA", "HA2", "HA3", "OXT", "HXT"})
CLASH_HEAVY = 2.2   # non-bonded heavy atoms closer than this clash (angstrom)
CLASH_HYDROGEN = 1.5
RNA_WINDOW = 4      # nucleotides relaxed together
RNA_STEPS = 200


class SequenceKind(StrEnum):
    """What a sequence describes and how it is built."""

    HELIX = "helix"              # peptide as an alpha helix
    STRAND = "strand"            # peptide as a beta strand
    DNA = "dna"                  # double-stranded B-DNA
    DNA_SINGLE = "dna single"    # one DNA strand
    RNA = "rna"                  # one RNA strand


_ALPHABETS = {
    SequenceKind.HELIX: set(AMINO_ACID_CODES),
    SequenceKind.STRAND: set(AMINO_ACID_CODES),
    SequenceKind.DNA: set("ACGT"),
    SequenceKind.DNA_SINGLE: set("ACGT"),
    SequenceKind.RNA: set("ACGU"),
}


def clean_sequence(text: str, kind: SequenceKind) -> str:
    """Turn typed or pasted text into a sequence of one-letter codes.

    FASTA header lines (starting with ``>``), whitespace, digits and the
    ``-`` and ``*`` characters of alignments are removed and the letters
    upper-cased.

    Parameters
    ----------
    text : str
        Sequence as typed or pasted.
    kind : SequenceKind
        Sets which letters are allowed.

    Returns
    -------
    str
        The sequence.

    Raises
    ------
    ValueError
        If the sequence is empty, too long or contains letters that are not
        valid for its kind.
    """
    lines = [line for line in text.splitlines() if not line.lstrip().startswith(">")]
    sequence = re.sub(r"[\s\d\-*]", "", "".join(lines)).upper()
    if not sequence:
        raise ValueError("The sequence is empty")
    bad = sorted(set(sequence) - _ALPHABETS[kind])
    if bad:
        what = "amino acid" if kind in (SequenceKind.HELIX, SequenceKind.STRAND) else "nucleotide"
        raise ValueError(f"Not {what} codes: {', '.join(bad)}")
    if len(sequence) > MAX_RESIDUES:
        raise ValueError(f"The sequence has {len(sequence)} residues; at most {MAX_RESIDUES} can be built")
    return sequence


def build_sequence(text: str, kind: SequenceKind) -> Molecule:
    """Build a 3D structure from a sequence.

    Parameters
    ----------
    text : str
        Sequence of one-letter codes (cleaned with :func:`clean_sequence`).
    kind : SequenceKind
        Peptide helix or strand, double- or single-stranded DNA, or RNA.

    Returns
    -------
    Molecule
        The structure with hydrogens and residue data (chain ``A``, and
        ``B`` for the second DNA strand).

    Raises
    ------
    ValueError
        If the sequence is not valid.
    """
    sequence = clean_sequence(text, kind)
    if kind == SequenceKind.HELIX:
        return build_peptide(sequence, *HELIX_PHI_PSI, ss="H")
    if kind == SequenceKind.STRAND:
        return build_peptide(sequence, *STRAND_PHI_PSI, ss="E" if len(sequence) >= 3 else "")
    return _build_nucleic_acid(sequence, kind)


# -- peptides -----------------------------------------------------------------------------------

def _tetrahedral_pair(center: np.ndarray, first: np.ndarray, second: np.ndarray,
                      length: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the two free tetrahedral positions on an atom with two bonds."""
    u, v = unit_vector(first - center), unit_vector(second - center)
    bisector = -unit_vector(u + v)
    normal = unit_vector(np.cross(u, v))
    a, b = math.cos(TETRAHEDRAL_HALF), math.sin(TETRAHEDRAL_HALF)
    return center + (bisector * a + normal * b) * length, center + (bisector * a - normal * b) * length


def _backbone(n: int, phi: float, psi: float) -> list[dict[str, np.ndarray]]:
    """Backbone atom positions (N, CA, C, O) of ``n`` residues with the given torsions."""
    residues: list[dict[str, np.ndarray]] = []
    n_pos = np.zeros(3)
    ca_pos = np.array([N_CA, 0.0, 0.0])
    c_pos = place_atom(ca_pos, n_pos, None, CA_C, ANGLE_N_CA_C, 0.0)
    for k in range(n):
        if k:
            prev = residues[-1]
            n_pos = place_atom(prev["C"], prev["CA"], prev["N"], C_N, ANGLE_CA_C_N, psi)
            ca_pos = place_atom(n_pos, prev["C"], prev["CA"], N_CA, ANGLE_C_N_CA, OMEGA)
            c_pos = place_atom(ca_pos, n_pos, prev["C"], CA_C, ANGLE_N_CA_C, phi)
        o_pos = place_atom(c_pos, ca_pos, n_pos, C_O, ANGLE_CA_C_O, psi + 180.0)
        residues.append({"N": n_pos, "CA": ca_pos, "C": c_pos, "O": o_pos})
    return residues


def build_peptide(sequence: str, phi: float, psi: float, ss: str = "") -> Molecule:
    """Build a neutral peptide with every residue at the same backbone torsions.

    Parameters
    ----------
    sequence : str
        One-letter amino acid codes.
    phi, psi : float
        Backbone torsions in degrees.
    ss : str, optional
        Secondary structure code stored on every atom (``"H"`` or ``"E"``).

    Returns
    -------
    Molecule
        Peptide with hydrogens, an amino group at the N-terminus and a
        carboxylic acid at the C-terminus.
    """
    backbone = _backbone(len(sequence), phi, psi)
    mol = Molecule()

    def add(name: str, symbol: str, position: np.ndarray, resseq: int) -> int:
        atom = Atom(symbol, name=name, resname="GLY", resseq=resseq, chain="A", ss=ss)
        return mol.add_atom(atom, position)

    previous_c: int | None = None
    for k, positions in enumerate(backbone):
        resseq = k + 1
        n = add("N", "N", positions["N"], resseq)
        ca = add("CA", "C", positions["CA"], resseq)
        c = add("C", "C", positions["C"], resseq)
        o = add("O", "O", positions["O"], resseq)
        mol.add_bond(n, ca, 1.0)
        mol.add_bond(ca, c, 1.0)
        mol.add_bond(c, o, 2.0)
        if previous_c is None:
            # N-terminal amino group: two hydrogens on a pyramidal nitrogen
            for name, torsion in (("H1", 180.0), ("H2", 60.0)):
                position = place_atom(positions["N"], positions["CA"], positions["C"], N_H, 109.5, torsion)
                mol.add_bond(n, add(name, "H", position, resseq), 1.0)
        else:
            mol.add_bond(previous_c, n, 1.0)
            # amide hydrogen in the peptide plane, opposite the bisector of C-N-CA
            n_pos = positions["N"]
            direction = unit_vector(unit_vector(n_pos - mol.coords[previous_c])
                                    + unit_vector(n_pos - positions["CA"]))
            mol.add_bond(n, add("H", "H", n_pos + direction * N_H, resseq), 1.0)
        for name, position in zip(("HA2", "HA3"), _tetrahedral_pair(positions["CA"], positions["N"],
                                                                    positions["C"], CA_H), strict=True):
            h = add(name, "H", position, resseq)
            mol.add_bond(ca, h, 1.0)
        previous_c = c

    # C-terminal carboxylic acid: OXT where the next nitrogen would be, with its hydrogen
    last = backbone[-1]
    resseq = len(backbone)
    oxt = add("OXT", "O", place_atom(last["C"], last["CA"], last["N"], C_OXT, ANGLE_CA_C_N, psi), resseq)
    assert previous_c is not None
    mol.add_bond(previous_c, oxt, 1.0)
    hxt = add("HXT", "H", place_atom(mol.coords[oxt], last["C"], last["O"], O_H, 107.0, 180.0), resseq)
    mol.add_bond(oxt, hxt, 1.0)
    infer_atom_types(mol, only_missing=False)

    for k, letter in enumerate(sequence):
        code = AMINO_ACID_CODES[letter]
        if code != "GLY":
            ca = next(i for i, a in enumerate(mol.atoms) if a.resseq == k + 1 and a.name == "CA")
            mutate_residue(mol, ca, code)
    for a in mol.atoms:
        a.ss = ss
    _relieve_clashes(mol)
    return mol


def _clashes(mol: Molecule) -> list[tuple[int, int]]:
    """Pairs of atoms that are neither bonded nor share a neighbour, yet are too close."""
    near: set[tuple[int, int]] = set(mol.bonds)
    for i, j in mol.bonds:
        near.update((min(i, k), max(i, k)) for k in mol.neighbors(j))
        near.update((min(j, k), max(j, k)) for k in mol.neighbors(i))
    pairs = []
    for i, j in pairs_within(mol.coords, CLASH_HEAVY).tolist():
        if (i, j) in near:
            continue
        limit = CLASH_HYDROGEN if "H" in (mol.atoms[i].symbol, mol.atoms[j].symbol) else CLASH_HEAVY
        if np.linalg.norm(mol.coords[i] - mol.coords[j]) < limit:
            pairs.append((i, j))
    return pairs


def _relieve_clashes(mol: Molecule) -> None:
    """Relax the side chains of residues that clash, with everything else held in place.

    Rigid side-chain templates cannot always avoid every contact; a proline
    in a helix, for example, puts a ring hydrogen against the carbonyl
    oxygen four residues back.
    """
    clashing = {(mol.atoms[k].chain, mol.atoms[k].resseq) for pair in _clashes(mol) for k in pair}
    if clashing:
        side_chains = [k for k, a in enumerate(mol.atoms)
                       if (a.chain, a.resseq) in clashing and a.name not in BACKBONE_NAMES]
        clean_selection(mol, side_chains)


# -- nucleic acids ------------------------------------------------------------------------------

def _build_nucleic_acid(sequence: str, kind: SequenceKind) -> Molecule:
    """DNA or RNA from Open Babel's FASTA reader (B-form helices)."""
    conversion = ob.OBConversion()
    conversion.SetInFormat("fasta")
    if kind == SequenceKind.DNA_SINGLE:
        conversion.AddOption("1", ob.OBConversion.INOPTIONS)
    title = "RNA" if kind == SequenceKind.RNA else "DNA"
    obmol = ob.OBMol()
    if not conversion.ReadString(obmol, f">{title}\n{sequence}\n") or obmol.NumAtoms() == 0:
        raise ValueError("Open Babel could not build this sequence")
    mol = from_openbabel(obmol, residues=True)
    if kind == SequenceKind.RNA:
        _relax_in_windows(mol)
    return mol


def _relax_in_windows(mol: Molecule) -> None:
    """Relax a strand with MMFF94, :data:`RNA_WINDOW` residues at a time.

    Each window is minimised together with the residues on either side of
    it, which are held fixed, so the cost grows linearly with the length
    instead of quadratically.
    """
    keys = list(dict.fromkeys((a.chain, a.resseq) for a in mol.atoms))
    members: dict[tuple[str, int], list[int]] = {}
    for i, a in enumerate(mol.atoms):
        members.setdefault((a.chain, a.resseq), []).append(i)
    for start in range(0, len(keys), RNA_WINDOW):
        window = keys[start:start + RNA_WINDOW]
        context = keys[max(0, start - 1):start] + keys[start + RNA_WINDOW:start + RNA_WINDOW + 1]
        moving = [i for key in window for i in members[key]]
        fixed = {i for key in context for i in members[key]}
        region = sorted(set(moving) | fixed)
        coords = minimize(mol.subset(region), fixed=[k for k, i in enumerate(region) if i in fixed],
                          steps=RNA_STEPS, force_fields=("MMFF94",))
        if coords is not None:
            position = {i: k for k, i in enumerate(region)}
            mol.set_positions(moving, coords[[position[i] for i in moving]])
