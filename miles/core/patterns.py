"""Finding substructures with SMARTS patterns.

A SMARTS pattern describes a group of atoms by elements, bonds and
properties such as charge, aromaticity or the number of attached
hydrogens; ``[CX3]=[OX1]`` is a carbonyl group. Open Babel does the
matching on a copy of the structure.

How hydrogens are counted depends on the structure. In fragments that
contain hydrogens (a ligand built in MILES, a structure from a
quantum chemistry program), the hydrogens are taken as complete: an oxygen
with one bond and no hydrogen is a negatively charged oxygen, as in a
carboxylate. Fragments without any hydrogen (proteins and nucleic acids
from most PDB files) get the usual number of hydrogens of each atom
implicitly, so ``[OX2H]`` still finds the hydroxyl groups of serines.
"""

from dataclasses import dataclass

from openbabel import openbabel as ob

from miles.core.chem import silence_openbabel, to_openbabel
from miles.core.molecule import Molecule


class PatternError(ValueError):
    """Raised for text that is not a valid SMARTS pattern."""


@dataclass(frozen=True)
class PatternMatch:
    """Where a pattern occurs in a structure.

    Attributes
    ----------
    matches : list of tuple of int
        Each occurrence as the atom indices matched by the pattern's atoms,
        in pattern order; occurrences that cover the same atoms are listed
        once.
    """

    matches: list[tuple[int, ...]]

    @property
    def atoms(self) -> list[int]:
        """Every matched atom, sorted, each once."""
        return sorted({i for match in self.matches for i in match})


def _fragments(mol: Molecule) -> list[list[int]]:
    """Split the atoms into connected groups."""
    seen = [False] * len(mol)
    groups: list[list[int]] = []
    neighbors = [mol.neighbors(i) for i in range(len(mol))]
    for start in range(len(mol)):
        if seen[start]:
            continue
        seen[start] = True
        group, stack = [], [start]
        while stack:
            i = stack.pop()
            group.append(i)
            for j in neighbors[i]:
                if not seen[j]:
                    seen[j] = True
                    stack.append(j)
        groups.append(group)
    return groups


def _prepared(mol: Molecule) -> "ob.OBMol":
    """Open Babel copy with formal charges, and implicit hydrogens where there are none (see above).

    :func:`~miles.core.chem.to_openbabel` gives every atom its formal charge
    (see :func:`~miles.core.charges.formal_charges`).
    """
    obmol = to_openbabel(mol)
    symbols = mol.symbols
    for group in _fragments(mol):
        if not any(symbols[i] == "H" for i in group):
            for i in group:
                ob.OBAtomAssignTypicalImplicitHydrogens(obmol.GetAtom(i + 1))
    return obmol


def compile_pattern(smarts: str) -> "ob.OBSmartsPattern":
    """Parse a SMARTS pattern.

    Parameters
    ----------
    smarts : str
        The pattern, e.g. ``"[CX3]=[OX1]"``.

    Returns
    -------
    openbabel.OBSmartsPattern
        The parsed pattern, ready to match.

    Raises
    ------
    PatternError
        If the text is empty or not valid SMARTS.
    """
    text = smarts.strip()
    pattern = ob.OBSmartsPattern()
    with silence_openbabel():
        valid = bool(text) and pattern.Init(text)
    if not valid:
        raise PatternError(f"Not a valid SMARTS pattern: {text}" if text else "Type a SMARTS pattern")
    return pattern


def match_pattern(mol: Molecule, smarts: str) -> PatternMatch:
    """Find every occurrence of a SMARTS pattern.

    Parameters
    ----------
    mol : Molecule
        Structure to search.
    smarts : str
        The pattern, e.g. ``"[OX2H]"`` for hydroxyl oxygens.

    Returns
    -------
    PatternMatch
        The occurrences, as indices of ``mol``'s atoms.

    Raises
    ------
    PatternError
        If ``smarts`` is not a valid pattern.
    """
    pattern = compile_pattern(smarts)
    if not len(mol):
        return PatternMatch([])
    obmol = _prepared(mol)
    with silence_openbabel():
        pattern.Match(obmol)
        found = [tuple(int(i) - 1 for i in match) for match in pattern.GetUMapList()]
    return PatternMatch(found)
