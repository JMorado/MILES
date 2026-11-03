"""Vibrational normal modes as frequency calculations report them.

Frequency calculations report, for every normal mode, its wavenumber, its
infrared intensity and the Cartesian displacement of every atom. Imaginary
frequencies, which mark a transition state or an unfinished optimisation,
are stored as negative wavenumbers.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class VibrationalMode:
    """One normal mode.

    Attributes
    ----------
    frequency : float
        Wavenumber in cm-1; negative for an imaginary frequency.
    intensity : float or None
        Infrared intensity in km/mol, or None if the calculation gave none.
    displacements : numpy.ndarray
        Displacement of every atom, shape ``(n_atoms, 3)``; only the
        direction and relative size matter.
    """

    frequency: float
    intensity: float | None
    displacements: np.ndarray

    @property
    def imaginary(self) -> bool:
        """Whether the frequency is imaginary."""
        return self.frequency < 0
