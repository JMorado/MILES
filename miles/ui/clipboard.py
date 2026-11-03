"""Copying molecules to and from the clipboard.

A copied molecule goes on the clipboard twice: in a private JSON format that
keeps bond orders and atom types exactly, and as XYZ text for other programs.
"""

import contextlib
import json

from PySide6.QtCore import QMimeData
from PySide6.QtWidgets import QApplication

from miles.core.molecule import Atom, Molecule
from miles.files.common import FileFormatError
from miles.files.formats import parse_text
from miles.files.openbabel import write_text

MIME_TYPE = "application/x-miles-molecule"


def copy_molecule(mol: Molecule) -> None:
    """Put a molecule on the clipboard."""
    data = {
        "symbols": mol.symbols,
        "coords": mol.coords.tolist(),
        "bonds": [[i, j, order] for (i, j), order in mol.bonds.items()],
        "types": [[a.geom, a.nsites, a.valence, a.inferred] for a in mol.atoms],
    }
    mime = QMimeData()
    mime.setData(MIME_TYPE, json.dumps(data).encode())
    with contextlib.suppress(FileFormatError):
        mime.setText(write_text(mol, "xyz"))
    QApplication.clipboard().setMimeData(mime)


def pasted_molecule() -> Molecule | None:
    """Molecule on the clipboard, in the private format or any text :func:`parse_text` reads.

    Returns
    -------
    Molecule or None
        None when the clipboard holds neither.

    Raises
    ------
    FileFormatError
        If the clipboard text is not a structure.
    """
    mime = QApplication.clipboard().mimeData()
    if mime is None:
        return None
    if mime.hasFormat(MIME_TYPE):
        data = json.loads(bytes(mime.data(MIME_TYPE).data()).decode())
        mol = Molecule()
        for symbol, position, (geom, nsites, valence, inferred) in zip(
                data["symbols"], data["coords"], data["types"], strict=True):
            mol.add_atom(Atom(symbol, geom, nsites, valence, inferred), position)
        for i, j, order in data["bonds"]:
            mol.add_bond(i, j, order)
        return mol
    if mime.hasText():
        return parse_text(mime.text())
    return None
