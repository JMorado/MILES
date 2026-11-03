"""Chemistry file formats through Open Babel (about 150 readable formats)."""

import gzip
import math
import os
import re
from collections import defaultdict, deque
from collections.abc import Iterable
from dataclasses import replace

import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel

from miles.core.chem import RESIDUE_FORMATS, from_openbabel, openbabel_coordinates, to_openbabel
from miles.core.elements import element, normalize_symbol
from miles.core.molecule import RESIDUE_KEY, Molecule
from miles.core.orbitals import OrbitalSet
from miles.core.units import ELECTRONVOLT, HARTREE
from miles.core.vibrations import VibrationalMode
from miles.core.volume import Volume
from miles.files import outputs
from miles.files.common import FileFormatError, JobInfo, extension
from miles.files.molden import molden_atoms, read_molden_orbitals
from miles.files.outputs import FCHK_FORMATS, GAMESS_FORMATS, GAUSSIAN_FORMATS

MAX_FRAMES = 20000

# Quantum chemistry outputs share the .log/.out extensions; the program is
# recognised from a phrase near the top of the file. The order matters:
# GAMESS-UK outputs also say "GAMESS", and other outputs may mention MOPAC.
QM_SIGNATURES = (
    ("Entering Gaussian System", "g09"),
    ("Gaussian, Inc.", "g09"),
    ("O   R   C   A", "orca"),
    ("Welcome to Q-Chem", "qcout"),
    ("Northwest Computational Chemistry Package", "nwo"),
    ("GAMESS-UK", "gukout"),
    ("GAMESS", "gamout"),
    ("Jaguar", "jout"),
    ("PROGRAM SYSTEM MOLPRO", "mpo"),
    ("MOPAC", "mopout"),
)
# Outputs of programs Open Babel cannot read (it finds no structure, or a
# wrong one), recognised the same way to say so.
UNREADABLE_SIGNATURES = (
    ("Psi4: An Open-Source Ab Initio", "Psi4"),
    ("TURBOMOLE", "Turbomole"),
    ("Dalton - An Electronic Structure Program", "DALTON"),
    ("MOLCAS", "Molcas"),
    ("x T B", "xtb"),
    ("Serenity", "Serenity"),
    ("CP2K", "CP2K"),
)
OUTPUT_EXTENSIONS = ("log", "out", "output")
QM_FORMATS = {code for _, code in QM_SIGNATURES} | GAUSSIAN_FORMATS | GAMESS_FORMATS | {"adfout"}
MOLDEN_FORMATS = {"mold", "molden", "molf"}
CUBE_FORMATS = {"cub", "cube"}
TEXT_FORMATS = QM_FORMATS | FCHK_FORMATS
"""Formats whose text is also read by :mod:`miles.files.outputs`."""
# Open Babel's GAMESS-UK reader aborts the whole program on the normal modes
# of a frequency calculation, so it is only given the text before them.
CUT_BEFORE = {"gukout": "cartesians to normal mode coordinates"}

CHARGE_FORMATS = {"sdf", "sd", "mol", "mdl", "cml", "smi", "smiles", "can", "ism", "inchi"}
"""Formats that state formal charges, so a charge of 0 there means neutral rather than not given.

PDB files state them too, in columns 79-80, but leave the field blank when
they do not; theirs are read from the records (see :func:`add_pdb_columns`).
"""
BOND_FORMATS = {"sdf", "sd", "mol", "mdl", "mol2", "ml2", "cml", "mrv"}
"""Formats whose bonds come from the file, so structures with other bonds are other molecules."""
# Readers that keep the energy of every optimisation step (OBConformerData),
# with the unit each one uses.
STEP_ENERGY_UNITS = (dict.fromkeys(("g03", "g09", "g16", "g92", "g94", "g98", "gal", "nwo"), 1.0)
                     | dict.fromkeys(("gam", "gamess", "gamout"), HARTREE))
# Energies on the comment lines of multi-structure XYZ files, for the
# programs whose format is fixed; anything else is left alone.
NUMBER = r"(-?\d+\.\d+(?:[eE][-+]?\d+)?)"
XYZ_ENERGIES = (
    (re.compile(r"\benergy:\s*" + NUMBER), HARTREE),                        # xtb
    (re.compile(r"ORCA-job.*\bE\s+" + NUMBER), HARTREE),                    # ORCA trajectories
    (re.compile(r"^(?=.*\bProperties=).*?\benergy=" + NUMBER), ELECTRONVOLT),   # extended XYZ (ASE)
    (re.compile(r"^\s*(-?\d+\.\d{4,})\s*$"), HARTREE),                     # CREST: the energy alone
)


def read_extensions() -> list[str]:
    """Extensions Open Babel can read.

    Returns
    -------
    list of str
        Format codes, which double as file extensions (``"pdb"``, ``"sdf"``, ...).
    """
    return sorted(pybel.informats)


def write_extensions() -> list[str]:
    """Extensions Open Babel can write.

    Returns
    -------
    list of str
        Format codes usable as file extensions.
    """
    return sorted(pybel.outformats)


def input_format(path: str) -> str | None:
    """Open Babel format code for a file, looking inside QM output files.

    Parameters
    ----------
    path : str
        File to inspect.

    Returns
    -------
    str or None
        Format code, or ``None`` if Open Babel has no reader for it.
    """
    ext = _inner_extension(path)
    if ext in OUTPUT_EXTENSIONS:
        head = _text(path, 200_000)
        for phrase, code in QM_SIGNATURES:
            if phrase in head:
                return code
        if any(phrase in head for phrase, _ in UNREADABLE_SIGNATURES):
            return None
    if ext in pybel.informats:
        return ext
    found = ob.OBConversion.FormatFromExt(path)
    if found is None or found.Flags() & ob.NOTREADABLE:
        return None   # e.g. Gaussian input, which Open Babel can only write
    return str(found.GetID())


def _inner_extension(path: str) -> str:
    """Extension of a file, looking through ``.gz`` (Open Babel reads gzipped files)."""
    ext = extension(path)
    return extension(path[:-3]) if ext == "gz" else ext


def _text(path: str, limit: int = -1) -> str:
    """Contents of a text file (gzipped or not), up to ``limit`` characters."""
    opener = gzip.open if extension(path) == "gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        return str(fh.read(limit))


def unreadable_program(path: str) -> str | None:
    """Name of the program that wrote an output Open Babel cannot read, if it is one of those known.

    Parameters
    ----------
    path : str
        A ``.log``, ``.out`` or ``.output`` file.

    Returns
    -------
    str or None
        E.g. ``"Psi4"``, or None.
    """
    if _inner_extension(path) not in OUTPUT_EXTENSIONS:
        return None
    head = _text(path, 200_000)
    if any(phrase in head for phrase, _ in QM_SIGNATURES):
        return None
    return next((name for phrase, name in UNREADABLE_SIGNATURES if phrase in head), None)


def _same_molecule_key(obmol: "ob.OBMol", bonds: bool) -> tuple[list[int], list[tuple[int, int, int]]]:
    """Return the elements in order, and the bonds when the file gives them: what makes one molecule."""
    elements = [a.GetAtomicNum() for a in ob.OBMolAtomIter(obmol)]
    if not bonds:
        return elements, []
    ends = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondOrder()) for b in ob.OBMolBondIter(obmol)]
    return elements, sorted((min(i, j), max(i, j), order) for i, j, order in ends)


def _frames_from(conv: "ob.OBConversion", first: "ob.OBMol", single: bool = False,
                 bonds: bool = False) -> tuple[np.ndarray | None, list[str], int]:
    """Coordinates of further structures in the file (or conformers) of the same molecule.

    A later structure is a frame of the first when it has the same elements
    in the same order and, for formats whose bonds come from the file
    (``bonds``), the same bonds. Other structures, such as the next
    molecule of an SDF library, are counted but not read.

    Returns
    -------
    frames : numpy.ndarray or None
        Coordinates of every frame, the first structure included; None for one structure.
    titles : list of str
        Title of each frame when the file holds several structures (empty for conformers).
    others : int
        Number of structures of other molecules after the first.
    """
    frames: list[np.ndarray] = []
    titles: list[str] = []
    others = 0
    if first.NumConformers() > 1:
        for k in range(first.NumConformers()):
            first.SetConformer(k)
            frames.append(openbabel_coordinates(first))
        first.SetConformer(first.NumConformers() - 1)
    else:
        frames.append(openbabel_coordinates(first))
        titles.append(first.GetTitle())
        molecule = _same_molecule_key(first, bonds)
        nxt = ob.OBMol()
        reads = MAX_FRAMES if single else 0
        # some readers (GAMESS and MOPAC input) keep reporting empty structures
        # at the end of the file, so an empty one ends the file; the reads are
        # capped too
        while reads < MAX_FRAMES and conv.Read(nxt) and nxt.NumAtoms() > 0:
            reads += 1
            if _same_molecule_key(nxt, bonds) == molecule:
                frames.append(openbabel_coordinates(nxt))
                titles.append(nxt.GetTitle())
            else:
                others += 1
            nxt = ob.OBMol()
    return (np.array(frames), titles, others) if len(frames) > 1 else (None, [], others)


def xyz_comment_energies(comments: list[str]) -> list[float] | None:
    """Energies written on the comment lines of an XYZ file with several structures.

    Recognises xtb (``energy: -5.07``), ORCA trajectories
    (``Coordinates from ORCA-job ... E -76.3``), extended XYZ
    (``energy=-12.3`` with ``Properties=``, in eV) and CREST conformer files
    (the energy alone). Every line must follow the same convention.

    Parameters
    ----------
    comments : list of str
        Second line of each structure.

    Returns
    -------
    list of float or None
        Energies in kcal/mol, or None when the lines hold no energy in a
        known format.
    """
    for pattern, unit in XYZ_ENERGIES:
        found = [pattern.search(text) for text in comments]
        if comments and all(found):
            return [float(m.group(1)) * unit for m in found if m is not None]
    return None


def frame_energies(obmol: "ob.OBMol", fmt: str, titles: list[str], count: int) -> list[float] | None:
    """Energy of each frame, when the file gives one for every frame.

    Parameters
    ----------
    obmol : openbabel.OBMol
        First structure as read by Open Babel (QM outputs keep the energies
        of the optimisation steps with it).
    fmt : str
        Open Babel format code.
    titles : list of str
        Title (comment line) of each structure of a multi-structure file.
    count : int
        Number of frames.

    Returns
    -------
    list of float or None
        Energies in kcal/mol, or None.
    """
    unit = STEP_ENERGY_UNITS.get(fmt)
    if unit is not None:
        data = [ob.toConformerData(d) for d in obmol.GetData() if d.GetDataType() == ob.ConformerData]
        values = [float(e) for e in data[0].GetEnergies()] if data else []
        if len(values) == count:
            return [e * unit for e in values]
    if fmt == "xyz" and len(titles) == count:
        return xyz_comment_energies(titles)
    return None


def _cannot_read(path: str, fmt: str | None, conv: "ob.OBConversion | None" = None) -> FileFormatError:
    """Make the error for a file nothing could be read from, saying as much as is known about why."""
    name = os.path.basename(path)
    program = unreadable_program(path)
    if program is not None:
        return FileFormatError(f"{name} is output of {program}, which MILES cannot read. Save the "
                               f"structure from {program} in another format, such as XYZ or Molden, "
                               "and open that instead.")
    if fmt is None:
        return FileFormatError(f"No reader for .{_inner_extension(path)} files")
    found = conv.GetInFormat() if conv is not None else None
    kind = str(found.Description()).splitlines()[0].strip() if found is not None else ""
    return FileFormatError(f"No structure found in {name}" + (f" (read as {kind})" if kind else ""))


PDB_RECORD_FORMATS = {"pdb", "ent", "pdbqt"}
"""Formats with PDB ATOM and HETATM records, whose occupancy, B-factor and alternate locations are read."""

_PDB_CHARGE = re.compile(r"^(\d)([+-])$|^([+-])(\d)$")


def pdb_records(lines: Iterable[str], charges: bool = True) -> list[dict]:
    """Read what the ATOM and HETATM records of a PDB file say beyond what Open Babel keeps.

    Parameters
    ----------
    lines : iterable of str
        Lines of the file; only the first model is read.
    charges : bool, optional
        Read the formal charge from columns 79-80 (not for PDBQT, which
        uses those columns for atom types).

    Returns
    -------
    list of dict
        One per record in file order: ``key`` (chain, residue number,
        insertion code, residue name, atom name), ``altloc``, ``occupancy``,
        ``b_factor`` (NaN when missing) and ``charge`` (None when the field
        is blank).
    """
    records = []
    for line in lines:
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        key = (line[21].strip(), resseq, line[26].strip(), line[17:20].strip(), line[12:16].strip())
        values = []
        for field in (line[54:60], line[60:66]):
            try:
                values.append(float(field))
            except ValueError:
                values.append(math.nan)
        charge = None
        found = _PDB_CHARGE.match(line[78:80].strip()) if charges else None
        if found is not None:
            digit = found.group(1) or found.group(4)
            sign = found.group(2) or found.group(3)
            charge = int(digit) * (1 if sign == "+" else -1)
        records.append({"key": key, "altloc": line[16].strip(), "occupancy": values[0],
                        "b_factor": values[1], "charge": charge})
    return records


def add_pdb_columns(source: str | Iterable[str], mol: Molecule, charges: bool = True) -> None:
    """Add what the ATOM and HETATM records give beyond Open Babel: occupancy, B-factor, altloc, charge.

    Open Babel's Python bindings give none of these back, so they are read
    from the records' fixed columns (first model only). Open Babel makes one
    atom of each record in file order; the records are matched to the atoms
    by chain, residue number, insertion code, residue name and atom name,
    taking repeated keys (alternate locations) in order. The occupancy and
    B-factor become the properties ``"occupancy"`` and ``"b_factor"``
    (NaN for atoms without a record; not added when no record has them),
    the alternate location goes to the ``altloc`` column, and a charge
    written in columns 79-80 becomes the atom's stated formal charge.

    Parameters
    ----------
    source : str or iterable of str
        The PDB file the molecule was read from, or its lines.
    mol : Molecule
        The molecule read from it, with residue data.
    charges : bool, optional
        Read formal charges (not for PDBQT files).
    """
    if isinstance(source, str):
        with open(source, encoding="utf-8", errors="replace") as fh:
            records = pdb_records(fh, charges)
    else:
        records = pdb_records(source, charges)
    if not records:
        return
    queues: dict[tuple, deque[dict]] = defaultdict(deque)
    for record in records:
        queues[record["key"]].append(record)
    columns = [mol.column(c).tolist() for c in (*RESIDUE_KEY, "name")]
    matched: list[dict | None] = []
    for chain, resseq, icode, resname, name in zip(*columns, strict=True):
        queue = queues.get((chain, resseq, icode, resname, name))
        matched.append(queue.popleft() if queue else None)
    rows = np.array([(m["occupancy"], m["b_factor"]) if m else (math.nan, math.nan) for m in matched],
                    dtype=float).reshape(-1, 2)
    if not np.isnan(rows).all():
        mol.set_property("occupancy", rows[:, 0])
        mol.set_property("b_factor", rows[:, 1])
    altloc = [m["altloc"] if m else "" for m in matched]
    if any(altloc):
        mol.set_column("altloc", altloc)
    stated = [(i, m["charge"]) for i, m in enumerate(matched) if m and m["charge"] is not None]
    if stated:
        atoms, values = zip(*stated, strict=True)
        mol.set_column("formal_charge", np.array(values, dtype=np.int8), list(atoms))


def read(path: str) -> tuple[Molecule, JobInfo]:
    """Read the first structure of a file, and any further ones as frames.

    Grids (cube files), normal modes (frequency calculations) and molecular
    orbitals (Molden files) are read too. For quantum chemistry outputs and
    formatted checkpoint files, :mod:`miles.files.outputs` adds what
    Open Babel leaves out: steps of a job, charge and multiplicity, and so on.

    Parameters
    ----------
    path : str
        File in any format Open Babel can read.

    Returns
    -------
    mol : Molecule
        The structure (the last geometry for QM outputs).
    info : JobInfo
        Title, charge, multiplicity, trajectory frames with their energies
        when the file has them, volumes, normal modes and orbitals.

    Raises
    ------
    FileFormatError
        If the format is unknown or the file holds no structure.
    """
    fmt = input_format(path)
    if fmt is None:
        raise _cannot_read(path, fmt)
    conv = ob.OBConversion()
    conv.SetInFormat(fmt)
    obmol = ob.OBMol()
    text = _text(path) if fmt in TEXT_FORMATS else None
    if text is not None and fmt in CUT_BEFORE:
        cut = text.find(CUT_BEFORE[fmt])
        ok = conv.ReadString(obmol, text if cut < 0 else text[:text.rfind("\n", 0, cut) + 1])
    elif fmt in MOLDEN_FORMATS and (fixed := molden_with_symbols(_text(path))) is not None:
        ok = conv.ReadString(obmol, fixed)
    else:
        ok = conv.ReadFile(obmol, path)
    if not ok or obmol.NumAtoms() == 0:
        raise _cannot_read(path, fmt, conv)
    frames, titles, others = _frames_from(conv, obmol, single=fmt in FCHK_FORMATS, bonds=fmt in BOND_FORMATS)
    mol = from_openbabel(obmol, residues=fmt in RESIDUE_FORMATS, charges=fmt in CHARGE_FORMATS)
    energies = None
    steps = outputs.job_steps(text, fmt, [element(s).number for s in mol.symbols]) if text else None
    if steps is not None:
        coords, step_energies = steps
        if frames is None or len(coords) > len(frames):
            frames = coords
        if step_energies is not None and len(step_energies) == len(frames):
            energies = step_energies
    if frames is not None:
        mol.coords = frames[-1 if fmt in QM_FORMATS else 0].copy()
    info = JobInfo(title=obmol.GetTitle().strip() or os.path.splitext(os.path.basename(path))[0],
                   charge=obmol.GetTotalCharge(), frames=frames, volumes=volumes(obmol),
                   vibrations=vibrations(obmol))
    if frames is not None:
        info.energies = energies if energies is not None else frame_energies(obmol, fmt, titles, len(frames))
    info.charges, info.dipole = file_charges(obmol)
    spin = obmol.GetTotalSpinMultiplicity()
    # only Open Babel's Gaussian reader takes the multiplicity from the file
    info.multiplicity = spin if spin > 0 and fmt in GAUSSIAN_FORMATS else None
    if text is not None:
        _add_output_details(text, fmt, mol, info)
    if fmt in PDB_RECORD_FORMATS:
        add_pdb_columns(path, mol, charges=fmt != "pdbqt")
    if others:
        info.notes.append(f"The file holds {others} more {'structure' if others == 1 else 'structures'} of "
                          "other molecules; only the first is shown.")
    if fmt in MOLDEN_FORMATS:
        try:
            info.orbitals = molden_orbitals(path, mol)
        except ValueError as exc:
            info.notes.append(f"The molecular orbitals were left out ({exc}).")
        if info.vibrations and all(m.intensity is None for m in info.vibrations):
            intensities = _molden_intensities(path)
            if len(intensities) == len(info.vibrations):
                info.vibrations = [replace(mode, intensity=intensity)
                                   for mode, intensity in zip(info.vibrations, intensities, strict=True)]
    if fmt in CUBE_FORMATS and len(info.volumes) == 1:
        info.volumes[0].name = _cube_comment(path) or info.volumes[0].name
    return mol, info


def _add_output_details(text: str, fmt: str, mol: Molecule, info: JobInfo) -> None:
    """Fill in what Open Babel leaves out of a quantum chemistry output or checkpoint file."""
    charge, multiplicity = outputs.charge_and_multiplicity(text, fmt)
    if fmt == "orca":
        modes = outputs.orca_modes(text, len(mol))
        if modes:
            info.vibrations = modes
    if fmt in FCHK_FORMATS:
        found = outputs.checkpoint(text, len(mol))
        charge, multiplicity = found.charge, found.multiplicity
        if found.frames is not None:
            info.frames, info.energies = found.frames, found.energies
            mol.coords = found.frames[-1].copy()
        if found.charges is not None:
            info.charges = found.charges
        info.dipole = found.dipole   # Open Babel gives the dipole in atomic units, not debye
        if not info.vibrations:
            info.vibrations = found.vibrations
    if charge is not None:
        info.charge = charge
    if multiplicity is not None:
        info.multiplicity = multiplicity


def molden_with_symbols(text: str) -> str | None:
    """Write the atom labels of the [Atoms] section of a Molden file as element symbols.

    Psi4 writes labels in capitals (``CU``), which Open Babel reads as
    dummy atoms.

    Parameters
    ----------
    text : str
        Contents of a Molden file.

    Returns
    -------
    str or None
        The text with labels such as ``CU`` or ``Cu1`` replaced by ``Cu``, or
        None when every label already is a symbol.
    """
    lines = text.splitlines(keepends=True)
    inside, changed = False, False
    for k, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("["):
            inside = stripped.lower().startswith("[atoms]")
            continue
        fields = line.split()
        if not inside or len(fields) < 6:
            continue
        symbol = normalize_symbol(fields[0])
        if symbol is not None and symbol != fields[0]:
            start = line.index(fields[0])
            lines[k] = line[:start] + symbol + line[start + len(fields[0]):]
            changed = True
    return "".join(lines) if changed else None


def _molden_intensities(path: str) -> list[float]:
    """Infrared intensities from the ``[INT]`` section of a Molden file.

    Open Babel reads the section only when it comes before the normal
    modes, but many programs write it after them.

    Parameters
    ----------
    path : str
        The Molden file.

    Returns
    -------
    list of float
        Intensities in file order (the first number on each line); empty when
        the file has no ``[INT]`` section.
    """
    values: list[float] = []
    inside = False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            text = line.strip()
            if text.startswith("["):
                if inside:
                    break
                inside = text.upper().startswith("[INT]")
                continue
            if inside and text:
                try:
                    values.append(float(text.split()[0]))
                except ValueError:
                    break
    return values


def _cube_comment(path: str) -> str:
    """Second line of a cube file, which says what the grid holds (the first is the job title)."""
    with open(path, encoding="utf-8", errors="replace") as f:
        f.readline()
        return f.readline().strip()


def _vector(v: "ob.vector3") -> np.ndarray:
    return np.array([v.GetX(), v.GetY(), v.GetZ()])


def volumes(obmol: "ob.OBMol") -> list[Volume]:
    """Read the grids attached to an Open Babel molecule (cube files hold one or more).

    Parameters
    ----------
    obmol : openbabel.OBMol
        Molecule as read by Open Babel.

    Returns
    -------
    list of Volume
        Every grid, in file order.
    """
    found = []
    grids = [ob.toGridData(d) for d in obmol.GetData() if d.GetDataType() == ob.GridData]
    for k, grid in enumerate(grids, start=1):
        axes = np.array([_vector(grid.GetXAxis()), _vector(grid.GetYAxis()), _vector(grid.GetZAxis())])
        origin = _vector(grid.GetOriginVector())
        # the step counts are not reachable from Python; the far corner gives them
        steps = np.rint(np.linalg.solve(axes.T, _vector(grid.GetMaxVector()) - origin)).astype(int) + 1
        values = np.asarray(grid.GetValues(), dtype=np.float32)
        if values.size != int(np.prod(steps)):
            continue
        name = grid.GetAttribute().strip() or f"Grid {k}"
        found.append(Volume(name, origin, axes, values.reshape(tuple(int(n) for n in steps))))
    return found


def file_charges(obmol: "ob.OBMol") -> tuple[np.ndarray | None, np.ndarray | None]:
    """Read the partial charges and dipole moment stored in a file, if any.

    Must be called before anything asks Open Babel for partial charges,
    which would otherwise compute Gasteiger charges and report those.

    Parameters
    ----------
    obmol : openbabel.OBMol
        Molecule as read by Open Babel.

    Returns
    -------
    charges : numpy.ndarray or None
        Charge of each atom (Mulliken charges of quantum chemistry outputs,
        MOL2 charges), or None when the file has none.
    dipole : numpy.ndarray or None
        Dipole moment in debye, or None when the file does not report one.
    """
    charges = None
    if obmol.HasPartialChargesPerceived():
        values = np.array([a.GetPartialCharge() for a in ob.OBMolAtomIter(obmol)], dtype=float)
        if values.any() and np.all(np.isfinite(values)):
            charges = values
    dipole = None
    data = obmol.GetData("Dipole Moment")
    if data is not None:
        dipole = _vector(ob.toVectorData(data).GetData())
    return charges, dipole


def vibrations(obmol: "ob.OBMol") -> list[VibrationalMode]:
    """Read the normal modes attached to an Open Babel molecule (frequency calculations).

    Parameters
    ----------
    obmol : openbabel.OBMol
        Molecule as read by Open Babel.

    Returns
    -------
    list of VibrationalMode
        Modes in file order; empty if the file has none.
    """
    data = [ob.toVibrationData(d) for d in obmol.GetData() if d.GetDataType() == ob.VibrationData]
    if not data:
        return []
    vib, n = data[0], obmol.NumAtoms()
    frequencies = list(vib.GetFrequencies())
    intensities = list(vib.GetIntensities())
    modes = []
    for k, (frequency, displacements) in enumerate(zip(frequencies, vib.GetLx(), strict=False)):
        vectors = np.array([_vector(v) for v in list(displacements)[:n]])
        if vectors.shape != (n, 3):
            continue
        intensity = float(intensities[k]) if k < len(intensities) else None
        modes.append(VibrationalMode(float(frequency), intensity, vectors))
    return _last_copy(modes)


def _last_copy(modes: list[VibrationalMode]) -> list[VibrationalMode]:
    """Keep the last of several copies of the same modes, which Open Babel joins into one list.

    Q-Chem prints the harmonic analysis again in an anharmonic calculation.
    Copies repeat the displacements too, unlike degenerate modes.
    """
    count = len(modes)
    for size in range(1, count // 2 + 1):
        if count % size == 0 and all(abs(modes[k].frequency - modes[k % size].frequency) < 1e-6
                                     and np.allclose(modes[k].displacements, modes[k % size].displacements)
                                     for k in range(count)):
            return modes[-size:]
    return modes


def molden_orbitals(path: str, mol: Molecule) -> OrbitalSet | None:
    """Read the molecular orbitals of a Molden file whose structure has been read.

    The orbitals are centred on the atoms of the [Atoms] section. If those
    are the atoms of ``mol`` in the same order, ``mol`` takes their positions
    so the orbitals line up with the structure; otherwise the orbitals are
    left out.

    Parameters
    ----------
    path : str
        Molden file.
    mol : Molecule
        Structure read from the same file; its coordinates may be updated.

    Returns
    -------
    OrbitalSet or None
        The orbitals, or None if the file has none.

    Raises
    ------
    ValueError
        If the orbitals cannot be read or do not belong to the structure;
        the message says why.
    """
    text = _text(path)
    found = read_molden_orbitals(text)
    listed = molden_atoms(text)
    if found is None or listed is None:
        return None
    labels, coords = listed
    if [normalize_symbol(label) for label in labels] != mol.symbols:
        raise ValueError("the atoms they are centred on are not those of the structure")
    if not np.allclose(coords, mol.coords, atol=1e-4):
        mol.coords = coords.copy()
    return found


def read_text(text: str, fmt: str) -> Molecule:
    """Read a structure from a string.

    Parameters
    ----------
    text : str
        File contents.
    fmt : str
        Open Babel format code.

    Returns
    -------
    Molecule
        The first structure in the text.

    Raises
    ------
    FileFormatError
        If nothing could be read.
    """
    conv = ob.OBConversion()
    if not conv.SetInFormat(fmt):
        raise FileFormatError(f"Unknown format {fmt}")
    obmol = ob.OBMol()
    if not conv.ReadString(obmol, text) or obmol.NumAtoms() == 0:
        raise FileFormatError("No structure found")
    mol = from_openbabel(obmol, residues=fmt in RESIDUE_FORMATS, charges=fmt in CHARGE_FORMATS)
    if fmt in PDB_RECORD_FORMATS:
        add_pdb_columns(text.splitlines(), mol, charges=fmt != "pdbqt")
    return mol


def write_text(mol: Molecule, fmt: str, title: str = "") -> str:
    """Write a structure to a string.

    Parameters
    ----------
    mol : Molecule
        Structure to write.
    fmt : str
        Open Babel format code, e.g. ``"xyz"``.
    title : str, optional
        Title line, for formats that have one.

    Returns
    -------
    str
        The file contents.

    Raises
    ------
    FileFormatError
        If Open Babel cannot write that format.
    """
    conv = ob.OBConversion()
    if not conv.SetOutFormat(fmt):
        raise FileFormatError(f"Cannot write {fmt}")
    return conv.WriteString(to_openbabel(mol, title))


def write(path: str, mol: Molecule, info: JobInfo) -> None:
    """Write a structure in the format given by the file extension.

    Parameters
    ----------
    path : str
        Destination; its extension selects the format.
    mol : Molecule
        Structure to write.
    info : JobInfo
        Supplies the title.

    Raises
    ------
    FileFormatError
        If Open Babel cannot write that format.
    """
    fmt = extension(path)
    conv = ob.OBConversion()
    if not conv.SetOutFormat(fmt):
        raise FileFormatError(f"Cannot write .{fmt} files")
    if not conv.WriteFile(to_openbabel(mol, info.title), path):
        raise FileFormatError(f"Could not write {path}")
    conv.CloseOutFile()
