"""Reading Gaussian input files (``.gjf``, ``.com``), which Open Babel can only write.

An input file has Link 0 lines (``%chk=...``), a route section (``#p
b3lyp/6-31g(d) opt``), a title, the charge and multiplicity, and the
molecule. Each atom of the molecule is given either by Cartesian
coordinates, optionally after a freeze code::

    C   0.000   0.000   0.000
    O  -1   0.000   0.000   1.210

or as a Z-matrix line, with atoms referred to by number or label and values
given as numbers or as variables defined after the molecule::

    C
    O  1  r1
    H  1  r2  2  a1
    H  1  r2  2  a1  3  180.0

    r1  1.21
    r2  1.10
    a1  121.0

Dummy atoms (``X``) and ghost atoms (``Bq``) place other atoms but are not
part of the structure. Only the first job of a file with several
(``--Link1--``) is read.
"""

import math
import os
import re

import numpy as np

from miles.core.elements import normalize_symbol
from miles.core.geometry import place_atom
from miles.core.molecule import Molecule
from miles.core.units import BOHR
from miles.files import openbabel
from miles.files.common import FileFormatError, JobInfo

PLACEHOLDERS = {"X", "Bq"}
"""Dummy and ghost atoms: used to place other atoms, then left out."""


def _is_number(token: str) -> bool:
    try:
        float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return False
    return True


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _element(label: str) -> str:
    """Element of an atom label such as ``C``, ``C12``, ``C(Fragment=1)``, ``C-CA--0.1`` or ``6``.

    Dummy atoms (``X``, ``X5``) give ``"X"`` and ghost atoms (``Bq``,
    ``C-Bq``) give ``"Bq"``.
    """
    base = label.split("(")[0]
    found = re.match(r"[A-Za-z]+", base.split("-")[0])
    letters = found.group().lower() if found else ""
    if letters == "bq" or base.lower().endswith("-bq"):
        return "Bq"
    if letters == "x":
        return "X"
    base = base.split("-")[0]
    symbol = normalize_symbol(base)
    if symbol is None:
        raise FileFormatError(f"Unknown atom {label!r}")
    return symbol


def _sections(text: str) -> tuple[str, str, list[str], list[str]]:
    """Route, title, molecule lines and variable lines of the first job."""
    job = re.split(r"^\s*--link1--\s*$", text, flags=re.I | re.M)[0]
    lines = [line.split("!")[0].rstrip() for line in job.replace("\r", "").splitlines()]
    k = 0
    while k < len(lines) and (not lines[k].strip() or lines[k].lstrip().startswith("%")):
        k += 1
    route: list[str] = []
    while k < len(lines) and lines[k].strip():
        route.append(lines[k].strip())
        k += 1
    if not route or not route[0].startswith("#"):
        raise FileFormatError("This is not a Gaussian input file: "
                              "it has no route section (a line starting with #)")
    blocks: list[list[str]] = []
    current: list[str] = []
    for line in lines[k + 1:]:
        if line.strip():
            current.append(line.strip())
        elif current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    if len(blocks) < 2:
        raise FileFormatError("The input has no molecule (it may take the geometry from a checkpoint file)")
    title, molecule = " ".join(blocks[0]), blocks[1]
    variables = [line for block in blocks[2:] for line in block]
    return " ".join(route), title, molecule, variables


def _variables(lines: list[str]) -> dict[str, float]:
    values = {}
    for line in lines:
        fields = [f for f in re.split(r"[\s=,]+", line) if f]
        if len(fields) >= 2 and _is_number(fields[1]):
            values[fields[0]] = _number(fields[1])
    return values


def parse(text: str) -> tuple[str, int, int, list[str], np.ndarray]:
    """Read the structure of a Gaussian input.

    Parameters
    ----------
    text : str
        Contents of the input file.

    Returns
    -------
    title : str
        The title section.
    charge, multiplicity : int
        From the line after the title.
    symbols : list of str
        Elements, without dummy and ghost atoms.
    coords : numpy.ndarray
        Positions in angstrom, shape ``(n_atoms, 3)``.

    Raises
    ------
    FileFormatError
        If the text is not a Gaussian input with a molecule in it.
    """
    route, title, molecule, variable_lines = _sections(text)
    head = molecule[0].replace(",", " ").split()
    if len(head) < 2 or not all(re.fullmatch(r"[-+]?\d+", f) for f in head[:2]):
        raise FileFormatError("The line after the title must give the charge and multiplicity")
    charge, multiplicity = int(head[0]), int(head[1])
    values = _variables(variable_lines)
    length = BOHR if re.search(r"units\s*[=(]\s*(au|bohr)", route, re.I) else 1.0
    degrees = not re.search(r"units\s*[=(][^)\s]*rad", route, re.I)
    labels: list[str] = []
    elements: list[str] = []
    positions: list[np.ndarray] = []

    def value(token: str) -> float:
        if _is_number(token):
            return _number(token)
        name = token.lstrip("+-")
        if name not in values:
            raise FileFormatError(f"The variable {name!r} is not defined")
        return -values[name] if token.startswith("-") else values[name]

    def known(token: str) -> bool:
        return _is_number(token) or token.lstrip("+-") in values

    def atom(token: str) -> np.ndarray:
        if re.fullmatch(r"\d+", token):
            index = int(token) - 1
        elif token in labels:
            index = labels.index(token)
        else:
            raise FileFormatError(f"Unknown atom {token!r} in the Z-matrix")
        if not 0 <= index < len(positions):
            raise FileFormatError(f"The Z-matrix refers to atom {token}, which is not defined before")
        return positions[index]

    for line in molecule[1:]:
        fields = [f for f in re.split(r"[\s,]+", line) if f]
        label, rest = fields[0], fields[1:]
        if len(rest) >= 4 and rest[0] in ("0", "-1") and all(known(f) for f in rest[1:4]):
            position = np.array([value(f) for f in rest[1:4]]) * length   # freeze code first; ONIOM adds more
        elif len(rest) == 3 and all(known(f) for f in rest):
            position = np.array([value(f) for f in rest]) * length
        elif len(rest) == 0 and not positions:
            position = np.zeros(3)
        elif len(rest) in (2, 4, 6, 7):
            r = value(rest[1]) * length
            if len(rest) == 2:
                position = atom(rest[0]) + np.array([0.0, 0.0, r])
            else:
                theta = value(rest[3]) if degrees else math.degrees(value(rest[3]))
                phi = 0.0
                if len(rest) >= 6:
                    phi = value(rest[5]) if degrees else math.degrees(value(rest[5]))
                position = place_atom(atom(rest[0]), atom(rest[2]), atom(rest[4]) if len(rest) >= 6 else None,
                                      r, theta, phi)
        else:
            raise FileFormatError(f"Cannot read the atom line {line!r}")
        labels.append(label)
        elements.append(_element(label))
        positions.append(np.asarray(position, dtype=float))
    keep = [k for k, symbol in enumerate(elements) if symbol not in PLACEHOLDERS]
    if not keep:
        raise FileFormatError("The input has no atoms")
    return title, charge, multiplicity, [elements[k] for k in keep], np.array([positions[k] for k in keep])


def read_text(text: str, name: str = "") -> tuple[Molecule, JobInfo]:
    """Structure, title, charge and multiplicity of a Gaussian input given as text.

    Bonds are worked out from the distances, as for an XYZ file.

    Parameters
    ----------
    text : str
        Contents of an input file.
    name : str, optional
        Title to use when the input has none.

    Returns
    -------
    mol : Molecule
        The structure.
    info : JobInfo
        Title, charge and multiplicity.
    """
    title, charge, multiplicity, symbols, coords = parse(text)
    rows = "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, (x, y, z) in zip(symbols, coords, strict=True))
    mol = openbabel.read_text(f"{len(symbols)}\n{title}\n{rows}", "xyz")
    return mol, JobInfo(title=title or name, charge=charge, multiplicity=multiplicity)


def read(path: str) -> tuple[Molecule, JobInfo]:
    """Read a Gaussian input file.

    Parameters
    ----------
    path : str
        A ``.gjf`` or ``.com`` file.

    Returns
    -------
    mol : Molecule
        The structure.
    info : JobInfo
        Title, charge and multiplicity.

    Raises
    ------
    FileFormatError
        If the file is not a Gaussian input with a molecule in it.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    return read_text(text, os.path.splitext(os.path.basename(path))[0])
