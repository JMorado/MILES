r"""Molecular orbitals and basis sets from Molden files.

A Molden file lists the atoms, the basis set ([GTO]) and the orbitals
([MO]) of a quantum chemistry calculation. Programs bend the format in a
few known ways, which the reader undoes:

* The format gives contraction coefficients for normalised primitives,
  but some programs (ORCA, xtb) write them for unnormalised ones.
* Cartesian d, f and g functions are normalised differently by CFOUR,
  Turbomole and Psi4 up to 1.3.2 (see :func:`_component_scales`).
* ORCA writes the spherical functions with ``|m|`` of 3 or 4 (in f, g and
  h shells) with the opposite sign; this is undone for files made by
  ``orca_2mkl``.

The reading is chosen by the orbitals themselves: the first one that makes
every orbital normalised is taken. The basis functions and orbitals are
those of :mod:`miles.core.orbitals`.
"""

import math
import re

import numpy as np

from miles.core.orbitals import (
    CARTESIAN,
    SHELL_TYPES,
    Orbital,
    OrbitalSet,
    Shell,
    Spin,
    cartesian_factor,
    contracted_shell,
    double_factorial,
    overlap_matrix,
    primitive_norm,
    self_overlap,
)
from miles.core.units import BOHR

NORM_TOLERANCE = 1e-4
"""Largest error in the norm of an orbital for a reading of a Molden file to be taken as right."""

NORM_FAILURE = 1e-2
"""Orbitals whose norms are off by more than this in every known reading are not read."""

CHECKED_ORBITALS = 300
"""Number of orbitals whose norms are checked (spread over the list) when a file has more."""

_SECTION = re.compile(r"^\s*\[([^\]]+)\](.*)$")

RawShell = tuple[np.ndarray, int, np.ndarray, np.ndarray]
"""Centre, angular momentum, exponents and coefficients as written in the file."""


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _sections(text: str) -> dict[str, tuple[str, list[str]]]:
    """Rest of the header line and body lines of each [Section], keyed by lower-case name."""
    found: dict[str, tuple[str, list[str]]] = {}
    body: list[str] | None = None
    for line in text.splitlines():
        match = _SECTION.match(line)
        if match:
            name = match.group(1).strip().lower()
            body = [] if name not in found else None   # keep the first of repeated sections
            if body is not None:
                found[name] = (match.group(2).strip(), body)
        elif body is not None:
            body.append(line)
    return found


def _atom_rows(header: str, lines: list[str]) -> tuple[list[str], np.ndarray]:
    """Labels and positions in bohr from the [Atoms] section."""
    rows = [line.split() for line in lines if len(line.split()) >= 6]
    scale = 1.0 / BOHR if "angs" in header.lower() else 1.0
    return [row[0] for row in rows], np.array([[_number(v) for v in row[3:6]] for row in rows]) * scale


def molden_atoms(text: str) -> tuple[list[str], np.ndarray] | None:
    """Atoms listed in the [Atoms] section of a Molden file.

    Parameters
    ----------
    text : str
        Contents of the file.

    Returns
    -------
    tuple of (list of str, numpy.ndarray) or None
        Atom labels as written (usually element symbols) and positions in
        angstrom, or None when there is no [Atoms] section.
    """
    sections = _sections(text)
    if "atoms" not in sections:
        return None
    labels, coords = _atom_rows(*sections["atoms"])
    return labels, coords.reshape(-1, 3) * BOHR


def _basis_rows(lines: list[str], atoms: np.ndarray) -> list[RawShell]:
    """Shells of a [GTO] section, as written."""
    shells: list[RawShell] = []
    tokens = [line.split() for line in lines]
    k, center = 0, np.zeros(3)
    while k < len(tokens):
        row = tokens[k]
        k += 1
        if not row:
            continue
        if row[0].isdigit():
            center = atoms[int(row[0]) - 1]
            continue
        kind, count = row[0].lower(), int(row[1])
        scale = _number(row[2]) if len(row) > 2 else 1.0
        data = np.array([[_number(v) for v in tokens[k + n][:3]] for n in range(count)])
        k += count
        exponents = data[:, 0] * scale**2
        if kind == "sp":
            shells.append((center, 0, exponents, data[:, 1]))
            shells.append((center, 1, exponents, data[:, 2]))
        elif kind in SHELL_TYPES:
            shells.append((center, SHELL_TYPES[kind], exponents, data[:, 1]))
        else:
            raise ValueError(f"Shells of type {row[0]!r} are not supported")
    return shells


def _primitives_normalized(shells: list[RawShell]) -> bool:
    """Whether the contraction coefficients refer to normalised primitives.

    Programs normalise their contracted functions in their own convention, so
    the reading that brings the contractions closer to unit norm is taken.
    Shells of a single primitive say nothing and are skipped.
    """
    balance = 0.0
    for _, angmom, exponents, coefficients in shells:
        if len(exponents) < 2:
            continue
        normalized = self_overlap(exponents, coefficients * primitive_norm(exponents, angmom), angmom)
        bare = self_overlap(exponents, coefficients, angmom)
        if normalized > 0.0 and bare > 0.0:
            balance += abs(math.log(normalized)) - abs(math.log(bare))
    return balance <= 0.0


def _orbital_rows(lines: list[str]) -> tuple[list[Orbital], list[dict[int, float]]]:
    """Orbitals and their sparse coefficient vectors (0-based basis index to value) from [MO]."""
    orbitals: list[Orbital] = []
    vectors: list[dict[int, float]] = []
    keys: dict[str, str] = {}
    current: dict[int, float] | None = None
    for line in lines:
        text = line.strip()
        if not text:
            continue
        if "=" in text:
            if current is not None:   # a keyword after coefficients starts the next orbital
                current, keys = None, {}
            key, _, value = text.partition("=")
            keys[key.strip().lower()] = value.strip()
            continue
        if current is None:
            spin = keys.get("spin", "").lower()
            orbitals.append(Orbital(_number(keys.get("ene", "0")), _number(keys.get("occup", "0")),
                                    Spin(spin) if spin in ("alpha", "beta") else Spin.NONE,
                                    keys.get("sym", "")))
            current = {}
            vectors.append(current)
        index, value = text.split()[:2]
        current[int(index) - 1] = _number(value)
    return orbitals, vectors


def read_molden_orbitals(text: str) -> OrbitalSet | None:
    """Molecular orbitals and basis set from the text of a Molden file.

    Parameters
    ----------
    text : str
        Contents of the file.

    Returns
    -------
    OrbitalSet or None
        None when the file has no [Atoms], [GTO] or [MO] section, or no orbitals.

    Raises
    ------
    ValueError
        If the sections cannot be read or do not fit together.
    """
    sections = _sections(text)
    if not {"atoms", "gto", "mo"} <= sections.keys():
        return None
    try:
        _, atoms = _atom_rows(*sections["atoms"])
        raw = _basis_rows(sections["gto"][1], atoms.reshape(-1, 3))
        orbitals, vectors = _orbital_rows(sections["mo"][1])
    except (IndexError, KeyError) as exc:
        raise ValueError(f"Malformed Molden file: {exc}") from exc
    if not orbitals or not raw:
        return None
    spherical = {2: bool(sections.keys() & {"5d", "5d7f", "5d10f"}),
                 3: bool(sections.keys() & {"5d", "5d7f", "7f"}),
                 4: "9g" in sections,
                 5: True}   # Molden has no flag for h shells; programs write them spherical
    guess = _primitives_normalized(raw)
    shells = [contracted_shell(center, angmom, exponents, weights, spherical.get(angmom, False), guess)
              for center, angmom, exponents, weights in raw]
    n_basis = sum(s.size for s in shells)
    coefficients = np.zeros((len(orbitals), n_basis))
    for k, vector in enumerate(vectors):
        for i, value in vector.items():
            if not 0 <= i < n_basis:
                raise ValueError(f"Orbital {k + 1} refers to basis function {i + 1}, but there are {n_basis}")
            coefficients[k, i] = value
    title = sections.get("title", ("", []))
    if "orca_2mkl" in " ".join([title[0], *title[1]]).lower():
        _undo_orca_signs(shells, coefficients)
    shells, scales = _choose_reading(raw, shells, coefficients, guess)
    return OrbitalSet(orbitals, shells, coefficients * scales, atoms.reshape(-1, 3))


CARTESIAN_CONVENTIONS = ("each", "product", "x**l", "shell")
"""Ways programs normalise cartesian d, f and g functions in Molden files (see :func:`_component_scales`)."""


def _component_scales(shell: Shell, convention: str) -> np.ndarray:
    """Factors that turn coefficients of a program's cartesian functions into coefficients of ours.

    ``"each"``: every function normalised, as the Molden format specifies.
    ``"product"``: the function ``x**a y**b z**c`` is
    ``sqrt((2a-1)!! (2b-1)!! (2c-1)!!)`` times a normalised one (CFOUR).
    ``"x**l"``: all functions of the shell share the normalisation of
    ``x**l`` (Psi4 up to 1.3.2). ``"shell"``: every function is
    ``sqrt((2l-1)!!)`` times a normalised one (Turbomole).
    """
    if convention == "each" or shell.spherical or shell.angmom < 2:
        return np.ones(shell.size)
    if convention == "x**l":
        return np.array([1.0 / cartesian_factor(*powers) for powers in CARTESIAN[shell.angmom]])
    if convention == "shell":
        return np.full(shell.size, math.sqrt(double_factorial(2 * shell.angmom - 1)))
    return np.array([math.sqrt(double_factorial(2 * a - 1) * double_factorial(2 * b - 1)
                               * double_factorial(2 * c - 1)) for a, b, c in CARTESIAN[shell.angmom]])


def _norm_error(coefficients: np.ndarray, overlap: np.ndarray) -> float:
    """Largest ``|<psi|psi> - 1|`` over a spread of the orbitals (rows of zeros are skipped)."""
    rows = np.flatnonzero(np.abs(coefficients).max(axis=1) > 0.0)
    if rows.size > CHECKED_ORBITALS:
        rows = rows[np.linspace(0, rows.size - 1, CHECKED_ORBITALS).astype(int)]
    if rows.size == 0:
        return 0.0
    chosen = coefficients[rows]
    norms = np.einsum("ij,ij->i", chosen @ overlap, chosen)
    return float(np.abs(norms - 1.0).max())


def _choose_reading(raw: list[RawShell], shells: list[Shell], coefficients: np.ndarray,
                    guess: bool) -> tuple[list[Shell], np.ndarray]:
    """Work out how the program wrote its basis functions, from the norms of the orbitals.

    Programs differ from the Molden format in two ways: the primitives may
    be unnormalised (ORCA, xtb), and cartesian functions may be normalised
    in another way (:func:`_component_scales`). The readings are tried in
    turn, starting with the Molden format itself, and the first one that
    makes every orbital normalised is taken; otherwise the closest one.

    Returns
    -------
    shells : list of Shell
        The basis set to use.
    scales : numpy.ndarray
        Factor for each basis function's coefficient, shape ``(n_basis,)``.

    Raises
    ------
    ValueError
        If no reading brings the orbitals close to normalised.
    """
    best: tuple[float, list[Shell], np.ndarray] | None = None
    for normalized in (guess, not guess):
        if normalized != guess:
            try:
                shells = [contracted_shell(center, angmom, exponents, weights, shell.spherical, normalized)
                          for (center, angmom, exponents, weights), shell in zip(raw, shells, strict=True)]
            except ValueError:
                continue
        overlap = overlap_matrix(shells)
        tried: list[np.ndarray] = []
        for convention in CARTESIAN_CONVENTIONS:
            scales = np.concatenate([_component_scales(shell, convention) for shell in shells])
            if any(np.array_equal(scales, other) for other in tried):
                continue   # no cartesian d, f or g shells: the conventions agree
            tried.append(scales)
            error = _norm_error(coefficients * scales, overlap)
            if error <= NORM_TOLERANCE:
                return shells, scales
            if best is None or error < best[0]:
                best = (error, shells, scales)
    assert best is not None
    error, shells, scales = best
    if error > NORM_FAILURE:
        raise ValueError(f"they are not normalised in any way of writing them that is known; "
                         f"the norm of an orbital is off by {error:.2g}")
    return shells, scales


def _undo_orca_signs(shells: list[Shell], coefficients: np.ndarray) -> None:
    """Flip the coefficients of spherical functions with ``|m|`` of 3 or 4, whose sign ORCA reverses.

    In Molden order (m = 0, +1, -1, +2, -2, ...) these are positions 5 to 8
    of f, g and h shells; h functions with ``|m| = 5`` keep their sign.
    """
    offset = 0
    for shell in shells:
        if shell.spherical and shell.angmom >= 3:
            coefficients[:, offset + 5:offset + min(shell.size, 9)] *= -1.0
        offset += shell.size
