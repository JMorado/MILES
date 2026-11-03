"""Results of a calculation read with a structure: grids, normal modes, orbitals, charges and dipole.

They belong to the structure they were computed for. :class:`Results` is
stamped with the topology token of that molecule
(:attr:`miles.core.molecule.Molecule.topology`); once an edit changes the
atoms or bonds, the token differs and the results are stale: the document
keeps them, so undo brings them back, but nothing uses them until the
structure is the one they were computed for again.
"""

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from miles.core.jobinfo import JobInfo
from miles.core.molecule import Molecule
from miles.core.orbitals import OrbitalSet
from miles.core.vibrations import VibrationalMode
from miles.core.volume import Volume


@dataclass(frozen=True, eq=False)
class Results:
    """What a calculation reported about one structure.

    The lists are treated as fixed once the results are made: the document
    replaces the whole object to add or remove a grid, so undo steps can
    keep the old one.

    Attributes
    ----------
    topology : int or None
        Topology token of the molecule the results belong to; None for no results.
    volumes : list of Volume
        Grids of values, e.g. an orbital or the density from a cube file.
    vibrations : list of VibrationalMode
        Normal modes of a frequency calculation.
    orbitals : OrbitalSet or None
        Molecular orbitals with their basis set.
    charges : numpy.ndarray or None
        Partial charge of each atom reported by the file.
    dipole : numpy.ndarray or None
        Dipole moment in debye.
    """

    topology: int | None = None
    volumes: list[Volume] = field(default_factory=list)
    vibrations: list[VibrationalMode] = field(default_factory=list)
    orbitals: OrbitalSet | None = None
    charges: np.ndarray | None = None
    dipole: np.ndarray | None = None

    @classmethod
    def from_info(cls, info: JobInfo, mol: Molecule) -> "Results":
        """Take the results a reader returned for a molecule (copied, so the reader's lists stay theirs)."""
        charges = None
        if info.charges is not None and len(info.charges) == len(mol):
            charges = np.array(info.charges, dtype=float)
        return cls(mol.topology, list(info.volumes), list(info.vibrations), info.orbitals, charges,
                   None if info.dipole is None else np.array(info.dipole, dtype=float))

    def is_empty(self) -> bool:
        """Whether there is nothing in them."""
        return not (self.volumes or self.vibrations or self.orbitals is not None
                    or self.charges is not None or self.dipole is not None)

    def belong_to(self, mol: Molecule) -> bool:
        """Whether they were computed for this molecule's atoms and bonds."""
        return self.topology is not None and self.topology == mol.topology

    def names(self) -> list[str]:
        """Return what they hold, in words: ``["normal modes", "orbitals", ...]``."""
        found = [("grids", bool(self.volumes)), ("normal modes", bool(self.vibrations)),
                 ("orbitals", self.orbitals is not None), ("charges", self.charges is not None),
                 ("dipole moment", self.dipole is not None)]
        return [name for name, present in found if present]

    def with_volumes(self, volumes: list[Volume], topology: int) -> "Results":
        """Return results with other grids, belonging to the molecule with ``topology``."""
        return replace(self, topology=topology, volumes=list(volumes))

    def as_info_fields(self) -> dict[str, Any]:
        """Return the fields of :class:`JobInfo` these results fill in."""
        return {"volumes": self.volumes, "vibrations": self.vibrations, "orbitals": self.orbitals,
                "charges": self.charges, "dipole": self.dipole}
