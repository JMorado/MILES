"""Writing volumetric data as Gaussian cube files.

Cube files are read by most programs that draw orbitals and densities, so a
grid on screen, including an orbital computed from a Molden file, can be
handed on. Reading goes through Open Babel (:mod:`miles.files.openbabel`).

The layout is the one Gaussian's ``cubegen`` writes: two comment lines (the
reader takes the second as the grid's name), the number of atoms and the
grid origin, the number of points and step vector of each axis, one line per
atom (atomic number, nuclear charge, position), then the values with the
last axis running fastest, six per line. Lengths are in bohr; the values are
written as they are stored.
"""


import numpy as np

from miles.core.elements import element
from miles.core.molecule import Molecule
from miles.core.units import BOHR
from miles.core.volume import Volume

PER_LINE = 6


def cube_text(mol: Molecule, volume: Volume, title: str = "") -> str:
    """Text of a cube file holding one grid and the atoms it belongs to.

    Parameters
    ----------
    mol : Molecule
        Atoms written into the header.
    volume : Volume
        The grid; its step vectors need not be orthogonal.
    title : str, optional
        First comment line, e.g. the document's name. The second comment
        line is the grid's name.

    Returns
    -------
    str
        The file's contents, ending with a newline.

    Raises
    ------
    ValueError
        If the values do not form a three-dimensional grid.
    """
    values = np.asarray(volume.values, dtype=np.float64)
    if values.ndim != 3:
        raise ValueError(f"A grid needs three axes, not {values.ndim}")
    origin = np.asarray(volume.origin, dtype=np.float64) / BOHR
    axes = np.asarray(volume.axes, dtype=np.float64) / BOHR
    lines = [_one_line(title) or "Grid", _one_line(volume.name) or "Grid",
             f"{len(mol):5d}" + "".join(f"{x:12.6f}" for x in origin)]
    for n, axis in zip(values.shape, axes, strict=True):
        lines.append(f"{n:5d}" + "".join(f"{x:12.6f}" for x in axis))
    for atom, position in zip(mol.atoms, np.asarray(mol.coords, dtype=np.float64) / BOHR, strict=True):
        z = element(atom.symbol).number
        lines.append(f"{z:5d}{float(z):12.6f}" + "".join(f"{x:12.6f}" for x in position))
    rows = values.reshape(-1, values.shape[2])
    for row in rows:
        for start in range(0, len(row), PER_LINE):
            lines.append("".join(f"{v:13.5E}" for v in row[start:start + PER_LINE]))
    return "\n".join(lines) + "\n"


def write_cube(path: str, mol: Molecule, volume: Volume, title: str = "") -> None:
    """Write one grid and the atoms it belongs to as a cube file.

    Parameters
    ----------
    path : str
        File to write (usually ending in ``.cube``).
    mol : Molecule
        Atoms written into the header.
    volume : Volume
        The grid.
    title : str, optional
        First comment line.

    Raises
    ------
    OSError
        If the file cannot be written.
    """
    with open(path, "w", encoding="ascii", errors="replace") as f:
        f.write(cube_text(mol, volume, title))


def _one_line(text: str) -> str:
    """Text fit for a comment line: one line of plain ASCII."""
    return " ".join(text.split()).encode("ascii", "replace").decode("ascii")
