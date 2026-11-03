"""Details of quantum chemistry outputs that Open Babel misses or reads wrongly.

Open Babel reads the structure of most outputs; these readers take from the
text of the file what it leaves out:

* the geometry of every step of an ORCA or Q-Chem job (optimisation, scan,
  dynamics) with its energy, of a GAMESS or Firefly optimisation, and of a
  Gaussian unrelaxed scan, where Open Babel keeps only some of them;
* the total charge and spin multiplicity, which Open Babel reads only from
  Gaussian outputs;
* the normal modes of ORCA, whose infrared intensities Open Babel takes
  from the wrong column in ORCA 5 and later, and which it cannot read at
  all when ORCA labels them by symmetry;
* Gaussian formatted checkpoint files: Mulliken charges, dipole moment,
  normal modes, charge and multiplicity, and every geometry of the files
  Q-Chem writes for an optimisation.
"""

import re
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from miles.core.elements import element, normalize_symbol
from miles.core.units import BOHR, HARTREE
from miles.core.units import DEBYE_PER_AU as DEBYE
from miles.core.vibrations import VibrationalMode

MAX_STEPS = 20000

GAUSSIAN_FORMATS = {"g03", "g09", "g16", "g92", "g94", "g98", "gal"}
GAMESS_FORMATS = {"gam", "gamess", "gamout"}
FCHK_FORMATS = {"fch", "fchk", "fck"}

Steps = tuple[np.ndarray, list[float] | None]
"""Coordinates of every step, shape ``(n_steps, n_atoms, 3)``, and their energies in kcal/mol (or None)."""

Row = Callable[[list[str]], tuple[int, list[float]]]
"""Reads the fields of an atom line into its atomic number and position."""

NUMBER = r"(-?\d+\.\d+)"
STEP_ENERGIES = {
    "orca": re.compile(r"FINAL SINGLE POINT ENERGY\s+" + NUMBER),
    "qcout": re.compile(r"Total energy in the final basis set =\s*" + NUMBER),
    **dict.fromkeys(GAUSSIAN_FORMATS, re.compile(r"SCF Done:\s+E\([^)]*\)\s*=\s*" + NUMBER)),
}
# A GAMESS step starts with "NSERCH=" (GAMESS-US prefixes "BEGINNING GEOMETRY
# SEARCH POINT"); the geometry printed after "EQUILIBRIUM GEOMETRY LOCATED"
# repeats the last step, and Firefly's "NSERCH= 3 ENERGY=" lines sum a step up.
GAMESS_STEP = re.compile(r"^\s*(?:BEGINNING GEOMETRY SEARCH POINT\s+)?NSERCH=\s*\d+\s*(?:\.\.\.)?\s*$")

CHARGE_PATTERNS = {
    # GAMESS-US prints "SPIN MULTIPLICITY", Firefly only its $CONTRL summary
    # ("MULT = 2  ICHARG = 1"); a bare MULT= may be that of excited states
    **dict.fromkeys(GAMESS_FORMATS, (r"CHARGE OF MOLECULE\s*=\s*(-?\d+)",
                                     r"SPIN MULTIPLICITY\s*=\s*(\d+)|\bMULT\s*=\s*(\d+)\s+ICHARG")),
    "gukout": (r"charge of molecule\s+(-?\d+)", r"state multiplicity\s+(\d+)"),
    "nwo": (r"^\s*Charge\s*:\s*(-?\d+)\s*$", r"Spin multiplicity:\s*(\d+)"),
    "orca": (r"Total Charge\s+Charge\s+\.+\s*(-?\d+)", r"Multiplicity\s+Mult\s+\.+\s*(\d+)"),
    "jout": (r"net molecular charge:\s*(-?\d+)", r"^\s*multiplicity:\s*(\d+)"),
    "mopout": (r"CHARGE ON SYSTEM\s*=\s*(-?\d+)", r"SPIN STATE DEFINED AS A (\w+)"),
}
SPIN_NAMES = {"SINGLET": 1, "DOUBLET": 2, "TRIPLET": 3, "QUARTET": 4, "QUINTET": 5, "SEXTET": 6,
              "SEPTET": 7, "OCTET": 8, "NONET": 9}


def _number_of(token: str) -> int:
    """Atomic number of an element symbol or label (0 when it is not an element)."""
    symbol = normalize_symbol(token)
    return element(symbol).number if symbol else 0


def _blocks(lines: list[str], header: str, skip: int, numbers: list[int], row: Row) -> list[np.ndarray]:
    """Coordinate blocks that follow every line containing ``header``.

    ``skip`` lines come between the header and the first atom; ``row``
    turns the fields of an atom line into its atomic number and position,
    or raises ValueError or IndexError. Blocks whose atoms differ from
    ``numbers`` are left out.
    """
    blocks = []
    count = len(numbers)
    for k, line in enumerate(lines):
        if header not in line:
            continue
        rows = lines[k + 1 + skip:k + 1 + skip + count]
        try:
            parsed = [row(text.split()) for text in rows]
        except (ValueError, IndexError):
            continue
        if len(parsed) == count and [z for z, _ in parsed] == numbers:
            blocks.append(np.array([xyz for _, xyz in parsed]))
    return blocks


def _gaussian_row(fields: list[str]) -> tuple[int, list[float]]:
    return int(fields[1]), [float(v) for v in fields[-3:]]


def _symbol_row(column: int) -> Row:
    def row(fields: list[str]) -> tuple[int, list[float]]:
        return _number_of(fields[column]), [float(v) for v in fields[column + 1:column + 4]]
    return row


def _gamess_row(fields: list[str]) -> tuple[int, list[float]]:
    return round(float(fields[1])), [float(v) for v in fields[2:5]]


def gamess_steps(text: str, numbers: list[int]) -> list[np.ndarray]:
    """Geometries of the steps of a GAMESS or Firefly optimisation or saddle point search.

    Open Babel keeps the energy of every step of a GAMESS output but only the
    last geometry. Every step starts with an ``NSERCH=`` line and prints its
    geometry under "COORDINATES OF ALL ATOMS ARE (ANGS)".

    Parameters
    ----------
    text : str
        The output.
    numbers : list of int
        Atomic numbers of the structure; blocks with other atoms are skipped.

    Returns
    -------
    list of numpy.ndarray
        Coordinates in angstrom of each step, shape ``(n_atoms, 3)``.
    """
    steps = []
    lines = text.splitlines()
    waiting = False
    for k, line in enumerate(lines):
        if GAMESS_STEP.match(line):
            waiting = True
        elif waiting and "COORDINATES OF ALL ATOMS ARE (ANGS)" in line:
            waiting = False
            steps += _blocks(lines[k:k + len(numbers) + 3], "COORDINATES OF ALL ATOMS ARE (ANGS)", 2,
                             numbers, _gamess_row)
    return steps


def job_steps(text: str, fmt: str, numbers: list[int]) -> Steps | None:
    """Geometry of every step of a job, with its energy where each step prints one.

    Parameters
    ----------
    text : str
        The output.
    fmt : str
        Open Babel format code of the output.
    numbers : list of int
        Atomic numbers of the structure; geometries with other atoms are skipped.

    Returns
    -------
    Steps or None
        None when the output has fewer than two steps or is not of ORCA,
        Q-Chem, GAMESS or Gaussian. Gaussian steps are only returned when
        every geometry has its own energy, as in an unrelaxed scan (an
        optimisation prints its last geometry twice).
    """
    lines = text.splitlines()
    if fmt in GAMESS_FORMATS:
        blocks = gamess_steps(text, numbers)
    elif fmt == "orca":
        blocks = _blocks(lines, "CARTESIAN COORDINATES (ANGSTROEM)", 1, numbers, _symbol_row(0))
    elif fmt == "qcout":
        blocks = _blocks(lines, "Standard Nuclear Orientation (Angstroms)", 2, numbers, _symbol_row(1))
    elif fmt in GAUSSIAN_FORMATS:
        header = "Standard orientation:" if "Standard orientation:" in text else "Input orientation:"
        blocks = _blocks(lines, header, 4, numbers, _gaussian_row)
    else:
        return None
    if len(blocks) < 2:
        return None
    pattern = STEP_ENERGIES.get(fmt)
    energies = [float(v) * HARTREE for v in pattern.findall(text)] if pattern is not None else []
    if len(energies) != len(blocks):
        if fmt in GAUSSIAN_FORMATS:
            return None
        energies = []
    return np.array(blocks[:MAX_STEPS]), energies[:MAX_STEPS] or None


def charge_and_multiplicity(text: str, fmt: str) -> tuple[int | None, int | None]:
    """Total charge and spin multiplicity stated in an output.

    Parameters
    ----------
    text : str
        The output.
    fmt : str
        Open Babel format code of the output (GAMESS, GAMESS-UK, NWChem,
        ORCA, Q-Chem, Jaguar and MOPAC are understood).

    Returns
    -------
    tuple of (int or None, int or None)
        The first charge and multiplicity the output states, each None when
        it states none.
    """
    if fmt == "qcout":
        found = re.search(r"^\s*\$molecule\s*\n\s*([-+]?\d+)\s+(\d+)\s*$", text, re.M | re.I)
        return (int(found.group(1)), int(found.group(2))) if found else (None, None)
    patterns = CHARGE_PATTERNS.get(fmt)
    if patterns is None:
        return None, None
    charge_found = re.search(patterns[0], text, re.M)
    spin_found = re.search(patterns[1], text, re.M)
    charge = int(charge_found.group(1)) if charge_found else None
    multiplicity = None
    if spin_found:
        value = next(g for g in spin_found.groups() if g is not None)
        multiplicity = SPIN_NAMES.get(value.upper()) if not value.isdigit() else int(value)
    return charge, multiplicity


# -- ORCA normal modes -----------------------------------------------------------------------------

ORCA_FREQUENCY = re.compile(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1")


def orca_modes(text: str, count: int) -> list[VibrationalMode] | None:
    """Read the normal modes of the last frequency calculation of an ORCA output.

    Parameters
    ----------
    text : str
        The output.
    count : int
        Number of atoms.

    Returns
    -------
    list of VibrationalMode or None
        The modes with non-zero frequencies, with infrared intensities in
        km/mol when the output has an IR spectrum in that unit (ORCA 5 and
        later); None when the output has no complete set of normal modes.
    """
    start = text.rfind("VIBRATIONAL FREQUENCIES")
    table = text.find("NORMAL MODES", start)
    if start < 0 or table < 0:
        return None
    frequencies = {}
    for line in text[start:table].splitlines():
        found = ORCA_FREQUENCY.match(line)
        if found:
            frequencies[int(found.group(1))] = float(found.group(2))
    size = 3 * count
    vectors = np.full((size, size), np.nan)
    columns: list[int] = []
    for line in text[table:].splitlines()[1:]:
        fields = line.split()
        if not fields:
            continue
        if "SPECTRUM" in line or "THERMOCHEMISTRY" in line:
            break
        if all(f.isdigit() for f in fields):
            columns = [int(f) for f in fields]
            continue
        if columns and fields[0].isdigit() and len(fields) == len(columns) + 1:
            try:
                values = [float(f) for f in fields[1:]]
            except ValueError:
                continue   # a row of symmetry labels
            row = int(fields[0])
            if row < size and max(columns) < size:
                vectors[row, columns] = values
    if np.isnan(vectors).any() or len(frequencies) != size:
        return None
    intensities = _orca_intensities(text[table:])
    return [VibrationalMode(frequency, intensities.get(k) if intensities else None,
                            vectors[:, k].reshape(count, 3))
            for k, frequency in sorted(frequencies.items()) if frequency != 0.0]


def _orca_intensities(text: str) -> dict[int, float]:
    """Infrared intensities in km/mol by mode number, from ORCA's IR SPECTRUM table (empty without one)."""
    start = text.find("IR SPECTRUM")
    if start < 0:
        return {}
    lines = text[start:].splitlines()[1:12]
    units = next((line.split() for line in lines if "km/mol" in line), None)
    if units is None:
        return {}
    column = 1 + units.index("km/mol")   # after the "6:" that numbers the mode
    found = {}
    for line in text[start:].splitlines()[1:]:
        fields = line.split()
        if fields and re.fullmatch(r"\d+:", fields[0]) and len(fields) > column:
            try:
                found[int(fields[0][:-1])] = float(fields[column])
            except ValueError:
                break
        elif found:
            break
    return found


# -- formatted checkpoint files --------------------------------------------------------------------

FCHK_RECORD = re.compile(r"^(\S.{0,40}?)\s+([IRCLH])\s+(N=)?\s*(\S+)\s*$")


def fchk_records(text: str) -> dict[str, list[Any]]:
    """Every record of a formatted checkpoint file, by name.

    Parameters
    ----------
    text : str
        Contents of the file.

    Returns
    -------
    dict
        Name to a list with one value per occurrence, in file order (Q-Chem
        repeats records for every step of a job): an int or float for a
        number, a numpy.ndarray for an integer or real array; text records
        are left out.
    """
    records: dict[str, list[Any]] = {}
    lines = text.splitlines()
    k = 2   # title and job type
    while k < len(lines):
        found = FCHK_RECORD.match(lines[k])
        k += 1
        if not found:
            continue
        name, kind, is_array, value = found.group(1).strip(), found.group(2), found.group(3), found.group(4)
        if not is_array:
            if kind == "I":
                records.setdefault(name, []).append(int(value))
            elif kind == "R":
                records.setdefault(name, []).append(float(value.replace("D", "E")))
            continue
        size = int(value)
        data: list[str] = []
        while k < len(lines) and lines[k].startswith(" "):
            if kind in "IR":
                data.extend(lines[k].split())
            k += 1
        if kind in "IR" and len(data) >= size:
            dtype = int if kind == "I" else float
            records.setdefault(name, []).append(np.array([dtype(v.replace("D", "E")) for v in data[:size]]))
    return records


@dataclass
class Checkpoint:
    """What a formatted checkpoint file adds to the structure Open Babel reads.

    Attributes
    ----------
    frames : numpy.ndarray or None
        Every geometry in angstrom when the file holds several (Q-Chem writes
        one per optimisation step).
    energies : list of float or None
        Energy of each frame in kcal/mol.
    charge, multiplicity : int or None
        Total charge and spin multiplicity.
    charges : numpy.ndarray or None
        Mulliken charges.
    dipole : numpy.ndarray or None
        Dipole moment in debye.
    vibrations : list of VibrationalMode
        Normal modes of a frequency calculation.
    """

    frames: np.ndarray | None = None
    energies: list[float] | None = None
    charge: int | None = None
    multiplicity: int | None = None
    charges: np.ndarray | None = None
    dipole: np.ndarray | None = None
    vibrations: list[VibrationalMode] = field(default_factory=list)


def checkpoint(text: str, count: int) -> Checkpoint:
    """Read what Open Babel leaves out of a formatted checkpoint file.

    Parameters
    ----------
    text : str
        Contents of the file.
    count : int
        Number of atoms.

    Returns
    -------
    Checkpoint
        Missing parts stay None (or empty).
    """
    records = fchk_records(text)
    found = Checkpoint()
    geometries = [g.reshape(-1, 3) * BOHR for g in records.get("Current cartesian coordinates", [])
                  if g.size == 3 * count]
    if len(geometries) > 1:
        found.frames = np.array(geometries[:MAX_STEPS])
        energies = records.get("Total Energy", [])
        if len(energies) == len(geometries):
            found.energies = [float(e) * HARTREE for e in energies[:MAX_STEPS]]
    if records.get("Charge"):
        found.charge = int(records["Charge"][0])
    if records.get("Multiplicity"):
        found.multiplicity = int(records["Multiplicity"][0])
    charges = records.get("Mulliken Charges", [])
    if charges and charges[-1].size == count and np.any(charges[-1]):
        found.charges = charges[-1].astype(float)
    dipole = records.get("Dipole Moment", [])
    if dipole and dipole[-1].size == 3:
        found.dipole = dipole[-1].astype(float) * DEBYE
    found.vibrations = _checkpoint_modes(records, count)
    return found


def _checkpoint_modes(records: dict[str, list[Any]], count: int) -> list[VibrationalMode]:
    """Read the normal modes from the Vib-E2 and Vib-Modes records.

    Vib-E2 holds the frequencies, reduced masses, force constants and
    infrared intensities of the modes, in blocks of one value per mode.
    """
    if not records.get("Vib-E2") or not records.get("Vib-Modes"):
        return []
    values, vectors = records["Vib-E2"][-1], records["Vib-Modes"][-1]
    modes = vectors.size // (3 * count)
    if modes == 0 or vectors.size != modes * 3 * count or values.size < 4 * modes:
        return []
    shapes = vectors.reshape(modes, count, 3)
    return [VibrationalMode(float(values[k]), float(values[3 * modes + k]), shapes[k].astype(float))
            for k in range(modes)]
