"""Sessions: a structure with its styles, measurements, frames, view, selection, grids, modes and orbitals.

A session file (``.miles``) is gzip-compressed JSON. Coordinates, frames,
grids, mode displacements and orbital coefficients are stored as
base64-encoded binary arrays: float64, so a session gives back exactly what
was saved, except grids, which are float32 in memory too. Every column of
the atom table is saved by name, so columns added to
:data:`~miles.core.molecule.ATOM_COLUMNS` later are saved without changes
here. Fields added after the first version are optional, so older sessions
still load.

Saving writes a temporary file next to the session and puts it in place
only once it is complete, so a failed save leaves the previous file as it
was.
"""

import base64
import contextlib
import gzip
import json
import os
import tempfile
from typing import Any

import numpy as np

from miles.core.molecule import ATOM_COLUMNS, Molecule
from miles.core.orbitals import Orbital, OrbitalSet, Shell, Spin
from miles.core.vibrations import VibrationalMode
from miles.core.volume import Volume
from miles.files.common import FileFormatError, JobInfo

EXTENSION = "miles"
FORMAT = "miles-session"
VERSION = 2
"""Version 2 stores coordinates as float64 and every atom column; version 1 files still load."""

SKIPPED_COLUMNS = ("number", "uid", "coords")
"""Atom columns not saved by name: elements go as symbols, coordinates as an array, identifiers are new."""

DISPLAY_FIELDS = ("rep", "hidden")
"""Per-atom display fields; sessions written before the display moved out of the molecule
keep them with the atoms."""


def _json_values(values: np.ndarray) -> list[Any]:
    """Values as JSON-ready lists, with None for NaN."""
    if values.dtype.kind == "f":
        return [None if v != v else v for v in values.tolist()]
    return values.tolist()


def _json_default(value: Any) -> Any:
    """Turn NumPy scalars and arrays, e.g. in plugin state, into plain JSON values."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} cannot be saved in a session")


def _columns_to_dict(mol: Molecule) -> dict[str, list[Any]]:
    """Every column of the atom table besides elements, coordinates and identifiers."""
    return {name: _json_values(mol.column(name)) for name in ATOM_COLUMNS if name not in SKIPPED_COLUMNS}


def _columns_from_dict(fields: dict[str, list[Any]], n: int) -> dict[str, list[Any]]:
    """Columns this version knows, with saved nulls turned back into each column's missing value."""
    columns = {}
    for name, values in fields.items():
        spec = ATOM_COLUMNS.get(name)
        if spec is None or name in SKIPPED_COLUMNS or len(values) != n:
            continue
        if spec.dtype.kind != "O":   # NaN was saved as null, and version 1 saved unset site counts as null
            values = [spec.default if v is None else v for v in values]
        columns[name] = values
    return columns


def _properties_to_dict(mol: Molecule) -> dict[str, Any]:
    """Property columns read from files or added by plugins (not those worked out, such as charges)."""
    saved = {}
    for name in mol.properties:
        if name in mol.property_sources:
            continue
        spec = mol.property_spec(name)
        default = spec.default
        if isinstance(default, float) and default != default:
            default = None
        saved[name] = {"dtype": spec.dtype.str if spec.dtype != object else "object", "default": default,
                       "values": _json_values(mol.column(name))}
    return saved


def _properties_from_dict(mol: Molecule, data: dict[str, Any]) -> None:
    for name, entry in data.items():
        dtype = np.dtype(object if entry["dtype"] == "object" else entry["dtype"])
        values = entry["values"]
        default = entry.get("default")
        if dtype.kind == "f":
            values = [np.nan if v is None else v for v in values]
            default = np.nan if default is None else default
        if len(values) == len(mol):
            mol.add_property(name, dtype, default, values)


def _pack(array: Any, dtype: Any = np.float64) -> dict[str, Any]:
    a = np.ascontiguousarray(np.asarray(array), dtype=dtype)
    data = base64.b64encode(a.tobytes()).decode("ascii")
    return {"shape": list(a.shape), "dtype": a.dtype.str, "data": data}


def _unpack(blob: dict[str, Any]) -> np.ndarray:
    raw = base64.b64decode(blob["data"])
    dtype = np.dtype(blob.get("dtype", "<f4"))   # version 1 stored float32
    return np.frombuffer(raw, dtype=dtype).reshape(blob["shape"]).astype(float)


def _orbitals_to_dict(orbitals: OrbitalSet) -> dict[str, Any]:
    return {
        "orbitals": [[o.energy, o.occupation, str(o.spin), o.symmetry] for o in orbitals.orbitals],
        "shells": [{"center": s.center.tolist(), "angmom": s.angmom, "exponents": s.exponents.tolist(),
                    "coefficients": s.coefficients.tolist(), "spherical": s.spherical}
                   for s in orbitals.shells],
        "coefficients": _pack(orbitals.coefficients),
        "atoms": np.asarray(orbitals.atoms, dtype=float).tolist(),
    }


def _orbitals_from_dict(data: dict[str, Any]) -> OrbitalSet:
    orbitals = [Orbital(float(energy), float(occupation), Spin(spin), symmetry)
                for energy, occupation, spin, symmetry in data["orbitals"]]
    shells = [Shell(np.array(s["center"], dtype=float), int(s["angmom"]),
                    np.array(s["exponents"], dtype=float), np.array(s["coefficients"], dtype=float),
                    bool(s["spherical"])) for s in data["shells"]]
    return OrbitalSet(orbitals, shells, _unpack(data["coefficients"]),
                      np.array(data["atoms"], dtype=float).reshape(-1, 3))


def to_dict(mol: Molecule, info: JobInfo) -> dict[str, Any]:
    """Describe a structure and its session state as JSON-ready data.

    Parameters
    ----------
    mol : Molecule
        Structure, including pinned measurements.
    info : JobInfo
        Title, charge, frames, results, display (per-atom styles) and view
        settings, as :meth:`miles.core.document.Document.session_info` gives them.

    Returns
    -------
    dict
        JSON-compatible dictionary.
    """
    frames = info.frames
    if frames is not None and not isinstance(frames, np.ndarray):
        frames = np.array([frames[k] for k in range(len(frames))])
    return {
        "format": FORMAT,
        "version": VERSION,
        "title": info.title,
        "charge": info.charge,
        "multiplicity": info.multiplicity,
        "symbols": mol.symbols,
        "coords": _pack(mol.coords),
        "atoms": _columns_to_dict(mol),
        "display": info.display,
        "properties": _properties_to_dict(mol),
        "bonds": [[i, j, order] for (i, j), order in mol.bonds.items()],
        "measurements": [list(m) for m in mol.measurements],
        "frames": None if frames is None else _pack(frames),
        "frame_labels": info.frame_labels,
        "energies": None if info.energies is None else [float(e) for e in info.energies],
        "view": info.view,
        "volumes": [{"name": v.name, "origin": np.asarray(v.origin, dtype=float).tolist(),
                     "axes": np.asarray(v.axes, dtype=float).tolist(), "values": _pack(v.values, np.float32),
                     "both_signs": v.both_signs}
                    for v in info.volumes],
        "vibrations": [{"frequency": m.frequency, "intensity": m.intensity,
                        "displacements": _pack(m.displacements)} for m in info.vibrations],
        "orbitals": None if info.orbitals is None else _orbitals_to_dict(info.orbitals),
        "charges": None if info.charges is None else np.asarray(info.charges, dtype=float).tolist(),
        "dipole": None if info.dipole is None else np.asarray(info.dipole, dtype=float).tolist(),
    }


def from_dict(data: dict[str, Any]) -> tuple[Molecule, JobInfo]:
    """Inverse of :func:`to_dict`.

    Parameters
    ----------
    data : dict
        Session dictionary.

    Returns
    -------
    mol : Molecule
        The structure.
    info : JobInfo
        Title, charge, frames, results, display and view settings.

    Raises
    ------
    FileFormatError
        If the dictionary is not a session, or a part of it is damaged.
    """
    if not isinstance(data, dict) or data.get("format") != FORMAT:
        raise FileFormatError("Not a session file")
    try:
        return _from_dict(data)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as exc:
        raise FileFormatError(f"The session file is damaged: {type(exc).__name__}: {exc}") from exc


def _from_dict(data: dict[str, Any]) -> tuple[Molecule, JobInfo]:
    symbols = list(data["symbols"])
    fields = data["atoms"]
    mol = Molecule.from_columns(_unpack(data["coords"]).reshape(len(symbols), 3), symbol=symbols,
                                **_columns_from_dict(fields, len(symbols)))
    display = data.get("display")
    if display is None and any(f in fields for f in DISPLAY_FIELDS):
        n = len(symbols)
        display = {"rep": fields.get("rep", [""] * n), "hidden": fields.get("hidden", [False] * n)}
    mol.set_bonds({(int(i), int(j)): float(order) for i, j, order in data["bonds"]})
    mol.measurements = [tuple(m) for m in data.get("measurements", [])]
    _properties_from_dict(mol, data.get("properties") or {})
    info = JobInfo(title=data.get("title", ""), charge=data.get("charge", 0),
                   multiplicity=data.get("multiplicity"), view=dict(data.get("view") or {}), display=display)
    if data.get("frames"):
        info.frames = _unpack(data["frames"])
        info.frame_labels = data.get("frame_labels")
        info.energies = data.get("energies")
    info.volumes = [Volume(v["name"], np.array(v["origin"], dtype=float), np.array(v["axes"], dtype=float),
                           _unpack(v["values"]).astype(np.float32), v.get("both_signs"))
                    for v in data.get("volumes", [])]
    info.vibrations = [VibrationalMode(float(m["frequency"]), m.get("intensity"), _unpack(m["displacements"]))
                       for m in data.get("vibrations", [])]
    if data.get("orbitals"):
        info.orbitals = _orbitals_from_dict(data["orbitals"])
    if data.get("charges") is not None:
        info.charges = np.array(data["charges"], dtype=float)
    if data.get("dipole") is not None:
        info.dipole = np.array(data["dipole"], dtype=float)
    return mol, info


def write(path: str, mol: Molecule, info: JobInfo) -> None:
    """Write a session file, replacing an older one only once the new one is complete.

    Parameters
    ----------
    path : str
        Destination (``.miles``).
    mol : Molecule
        Structure to save.
    info : JobInfo
        Title, frames, results and view settings to save with it.

    Raises
    ------
    FileFormatError
        If something in the session cannot be written (the old file is kept).
    OSError
        If the file cannot be written.
    """
    try:
        text = json.dumps(to_dict(mol, info), separators=(",", ":"), default=_json_default)
    except (TypeError, ValueError) as exc:
        raise FileFormatError(f"Could not save the session: {exc}") from exc
    target = os.path.abspath(path)
    mode = (os.stat(target).st_mode & 0o777) if os.path.exists(target) else 0o644
    handle, temporary = tempfile.mkstemp(prefix=".", suffix=".miles-part", dir=os.path.dirname(target))
    try:
        with os.fdopen(handle, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw) as fh:
            fh.write(text.encode("utf-8"))
        os.chmod(temporary, mode)
        os.replace(temporary, target)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(temporary)
        raise


def read(path: str) -> tuple[Molecule, JobInfo]:
    """Read a session file.

    Parameters
    ----------
    path : str
        Session file.

    Returns
    -------
    mol : Molecule
        The structure.
    info : JobInfo
        Title, frames, results and view settings.

    Raises
    ------
    FileFormatError
        If the file is not a session or is damaged.
    """
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError, EOFError) as exc:
        raise FileFormatError(f"Not a session file: {exc}") from exc
    return from_dict(data)
