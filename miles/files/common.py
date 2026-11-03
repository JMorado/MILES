"""Types shared by the file readers and writers."""

import os

# JobInfo is documented in miles.core.jobinfo; the readers import it from here
from miles.core.jobinfo import JobInfo as JobInfo
from miles.core.molecule import Molecule

__all__ = ["FileFormatError", "default_multiplicity", "extension"]


class FileFormatError(ValueError):
    """Raised when a file or text cannot be read or written."""


def extension(path: str) -> str:
    """Lower-case file extension without the dot.

    Parameters
    ----------
    path : str
        File name or path.

    Returns
    -------
    str
        Extension such as ``"pdb"``, or ``""`` when there is none.
    """
    return os.path.splitext(path)[1].lower().lstrip(".")


def default_multiplicity(mol: Molecule, charge: int = 0) -> int:
    """Lowest spin multiplicity compatible with the electron count.

    Parameters
    ----------
    mol : Molecule
        Structure whose electrons are counted.
    charge : int, optional
        Total charge.

    Returns
    -------
    int
        1 for an even number of electrons, otherwise 2.
    """
    return 1 if mol.electron_count(charge) % 2 == 0 else 2
