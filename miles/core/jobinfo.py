"""What a reader returns with a structure besides its atoms: title, charge, frames, grids, modes, orbitals.

A :class:`~miles.core.document.Document` takes copies of what it needs from
a ``JobInfo`` and leaves it as it was: the trajectory becomes a
:class:`~miles.core.frames.Trajectory`, the grids, modes, orbitals,
charges and dipole become :class:`~miles.core.results.Results`.
"""

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from miles.core.orbitals import OrbitalSet
from miles.core.vibrations import VibrationalMode
from miles.core.volume import Volume


@dataclass
class JobInfo:
    """Information that travels with a structure besides its atoms.

    Attributes
    ----------
    title : str
        Name or title line of the structure.
    charge : int
        Total charge.
    multiplicity : int or None
        Spin multiplicity; ``None`` means "work it out from the electron count".
    frames : array_like or None
        Trajectory coordinates with shape ``(n_frames, n_atoms, 3)``.
    frame_labels : list of str or None
        Optional text for each frame, e.g. conformer energies.
    energies : list of float or None
        Energy of each frame in kcal/mol, e.g. the steps of an optimisation.
    view : dict
        Saved view settings (camera, styles), filled by sessions.
    volumes : list of Volume
        Grids of values read with the structure, e.g. orbitals from a cube file.
    vibrations : list of VibrationalMode
        Normal modes from a frequency calculation, for the geometry shown.
    orbitals : OrbitalSet or None
        Molecular orbitals with their basis set, e.g. from a Molden file.
    charges : numpy.ndarray or None
        Partial charge of each atom from the file (Mulliken charges of a
        quantum chemistry output, MOL2 charges), in elementary charges.
    dipole : numpy.ndarray or None
        Dipole moment in debye reported by the calculation, for the geometry read.
    notes : list of str
        What the reader could not take from the file although it was there,
        e.g. orbitals in a convention it does not know, as sentences to show
        the user.
    display : dict or None
        Per-atom representation and focus-mode hiding saved by a session
        (see :meth:`miles.core.display.AtomDisplay.to_lists`).
    """

    title: str = ""
    charge: int = 0
    multiplicity: int | None = None
    frames: Any = None
    frame_labels: list[str] | None = None
    energies: list[float] | None = None
    view: dict[str, Any] = field(default_factory=dict)
    volumes: list[Volume] = field(default_factory=list)
    vibrations: list[VibrationalMode] = field(default_factory=list)
    orbitals: OrbitalSet | None = None
    charges: np.ndarray | None = None
    dipole: np.ndarray | None = None
    notes: list[str] = field(default_factory=list)
    display: dict[str, list[Any]] | None = None
