"""Animating normal modes, drawing them as arrows, and their infrared spectrum."""

from collections.abc import Sequence

import numpy as np

from miles import VibrationalMode


def mode_frames(coords: np.ndarray, mode: VibrationalMode, amplitude: float = 0.4,
                count: int = 24) -> np.ndarray:
    """Coordinates for one period of a vibration, for animation.

    Parameters
    ----------
    coords : numpy.ndarray
        Equilibrium coordinates, shape ``(n_atoms, 3)``.
    mode : VibrationalMode
        The mode.
    amplitude : float, optional
        Largest displacement of any atom, in angstrom.
    count : int, optional
        Number of frames; frame ``count`` would equal frame 0, so the frames
        loop smoothly.

    Returns
    -------
    numpy.ndarray
        Coordinates, shape ``(count, n_atoms, 3)``: the equilibrium geometry
        displaced by ``sin(2 pi k / count)`` times the scaled mode.
    """
    coords = np.asarray(coords, dtype=float)
    displacements = np.asarray(mode.displacements, dtype=float)
    largest = float(np.linalg.norm(displacements, axis=1).max()) if len(displacements) else 0.0
    scaled = displacements * (amplitude / largest) if largest > 0 else np.zeros_like(coords)
    phases = np.sin(2.0 * np.pi * np.arange(count) / count)
    return coords[None, :, :] + phases[:, None, None] * scaled[None, :, :]


def ir_spectrum(modes: Sequence[VibrationalMode], wavenumbers: np.ndarray, width: float = 20.0) -> np.ndarray:
    """Infrared absorbance broadened with Lorentzians.

    Parameters
    ----------
    modes : sequence of VibrationalMode
        Normal modes; imaginary modes and modes without intensity are skipped.
    wavenumbers : numpy.ndarray
        Wavenumbers (cm-1) at which to evaluate the spectrum.
    width : float, optional
        Full width at half maximum of each band, cm-1.

    Returns
    -------
    numpy.ndarray
        Absorbance at each wavenumber, scaled so the highest point is 1;
        zeros when no mode has an intensity.
    """
    x = np.asarray(wavenumbers, dtype=float)
    spectrum = np.zeros_like(x)
    half = width / 2.0
    for mode in modes:
        if mode.imaginary or not mode.intensity:
            continue
        spectrum += mode.intensity * half ** 2 / ((x - mode.frequency) ** 2 + half ** 2)
    peak = float(spectrum.max()) if spectrum.size else 0.0
    return spectrum / peak if peak > 0 else spectrum


def mode_arrows(mode: VibrationalMode, length: float = 1.1, threshold: float = 0.05) -> np.ndarray:
    """Displacement vectors of a mode scaled for drawing as arrows on the atoms.

    Parameters
    ----------
    mode : VibrationalMode
        The mode.
    length : float, optional
        Length in angstrom of the arrow of the atom that moves most.
    threshold : float, optional
        Atoms moving less than this fraction of the largest displacement get
        no arrow (a zero vector), so the picture shows what moves.

    Returns
    -------
    numpy.ndarray
        Vectors in angstrom, shape ``(n_atoms, 3)``; all zero for a mode
        without displacements.
    """
    displacements = np.asarray(mode.displacements, dtype=float).reshape(-1, 3)
    sizes = np.linalg.norm(displacements, axis=1)
    largest = float(sizes.max()) if len(sizes) else 0.0
    if largest == 0.0:
        return np.zeros_like(displacements)
    arrows = displacements * (length / largest)
    arrows[sizes < threshold * largest] = 0.0
    return arrows
