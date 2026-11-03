"""Whether a structure is a small molecule and whether it lacks hydrogens.

Properties and 2D drawings both check these first.
"""

import numpy as np

from miles import AMINO_ACIDS, NUCLEOTIDES, Molecule

MAX_HEAVY_ATOMS = 250
"""Structures with more heavy atoms get their composition only, and no 2D drawing."""

MIN_POLYMER_RESIDUES = 3
"""A structure with this many amino acids or nucleotides is treated as a polymer."""


def is_small_molecule(mol: Molecule) -> bool:
    """Whether a structure is small enough for identifiers and descriptors.

    Parameters
    ----------
    mol : Molecule
        Structure to judge.

    Returns
    -------
    bool
        False for proteins and nucleic acids (three or more residues) and for
        structures with more than :data:`MAX_HEAVY_ATOMS` heavy atoms.
    """
    heavy = int((mol.numbers != 1).sum())
    names = mol.column("name")
    guides = np.flatnonzero(((names == "CA") & mol.mask("resname", AMINO_ACIDS))
                            | ((names == "P") & mol.mask("resname", NUCLEOTIDES)))
    chains, numbers = mol.column("chain")[guides].tolist(), mol.column("resseq")[guides].tolist()
    residues = set(zip(chains, numbers, strict=True))
    return heavy <= MAX_HEAVY_ATOMS and len(residues) < MIN_POLYMER_RESIDUES


def hydrogens_missing(mol: Molecule) -> bool:
    """Whether some carbon atom has fewer bonds than its four valences.

    Aromatic bonds count 1.5, so ring carbons of benzene and its relatives
    come out complete.

    Parameters
    ----------
    mol : Molecule
        Structure to check.

    Returns
    -------
    bool
        True if any carbon's bond orders add up to less than 3.5.
    """
    valence = np.zeros(len(mol))
    np.add.at(valence, mol.bond_pairs.reshape(-1), np.repeat(mol.bond_orders, 2))
    return bool(np.any(valence[mol.numbers == 6] < 3.5))
