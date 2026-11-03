"""Quantities along a trajectory: distances, angles, dihedrals and RMSD in every frame.

The functions work on a whole stack of frames at once, shape
``(n_frames, n_atoms, 3)``, so a long trajectory is measured in one NumPy
pass. Angles are in degrees, lengths in angstrom.
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike

MIN_RMSD_ATOMS = 3
"""RMSD uses the heavy atoms when there are at least this many, otherwise all atoms."""


def _unit(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(norm > 0, norm, 1.0)


def measurement_series(frames: ArrayLike, atoms: Sequence[int]) -> np.ndarray:
    """Distance, angle or dihedral of 2, 3 or 4 atoms in every frame.

    Parameters
    ----------
    frames : array_like
        Coordinates, shape ``(n_frames, n_atoms, 3)``.
    atoms : sequence of int
        Two atoms for a distance, three for an angle (vertex in the middle)
        or four for a dihedral.

    Returns
    -------
    numpy.ndarray
        One value per frame: angstrom for distances, degrees for angles, and
        degrees in (-180, 180] for dihedrals (see :func:`unwrap_degrees`).

    Raises
    ------
    ValueError
        If ``atoms`` does not have 2, 3 or 4 entries.
    """
    x = np.asarray(frames, dtype=float)[:, list(atoms)]
    if len(atoms) == 2:
        return np.linalg.norm(x[:, 0] - x[:, 1], axis=-1)
    if len(atoms) == 3:
        u, v = _unit(x[:, 0] - x[:, 1]), _unit(x[:, 2] - x[:, 1])
        return np.degrees(np.arccos(np.clip((u * v).sum(axis=-1), -1.0, 1.0)))
    if len(atoms) == 4:
        b1 = _unit(x[:, 2] - x[:, 1])
        p, s = x[:, 0] - x[:, 1], x[:, 3] - x[:, 2]
        v = p - (p * b1).sum(axis=-1, keepdims=True) * b1
        w = s - (s * b1).sum(axis=-1, keepdims=True) * b1
        return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(axis=-1), (v * w).sum(axis=-1)))
    raise ValueError(f"A measurement needs 2, 3 or 4 atoms, not {len(atoms)}")


def unwrap_degrees(values: ArrayLike) -> np.ndarray:
    """Remove the jumps of 360 degrees from a series of dihedrals.

    A dihedral that turns past 180 degrees jumps to -180; unwrapped, it
    goes on to 181, 182 and so on, so a plot draws one continuous line.

    Parameters
    ----------
    values : array_like
        Angles in degrees, in order.

    Returns
    -------
    numpy.ndarray
        The same angles with multiples of 360 added so that consecutive
        values differ by less than 180, shifted as a whole so their mean
        lies within [-180, 180].
    """
    unwrapped = np.degrees(np.unwrap(np.radians(np.asarray(values, dtype=float))))
    if not unwrapped.size:
        return unwrapped
    return unwrapped - 360.0 * np.round(unwrapped.mean() / 360.0)


def rmsd_atoms(symbols: Sequence[str]) -> list[int]:
    """Atoms the RMSD of a trajectory is computed over.

    Parameters
    ----------
    symbols : sequence of str
        Element symbols.

    Returns
    -------
    list of int
        The heavy atoms, or every atom when there are fewer than
        :data:`MIN_RMSD_ATOMS` heavy atoms.
    """
    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    return heavy if len(heavy) >= MIN_RMSD_ATOMS else list(range(len(symbols)))


def rmsd_series(frames: ArrayLike, atoms: Sequence[int] | None = None, reference: int = 0) -> np.ndarray:
    """Root-mean-square deviation of every frame from a reference frame.

    Each frame is first moved onto the reference by the rotation and
    translation that fit it best (Kabsch), so only changes of shape count,
    not the molecule drifting or turning.

    Parameters
    ----------
    frames : array_like
        Coordinates, shape ``(n_frames, n_atoms, 3)``.
    atoms : sequence of int, optional
        Atoms compared; all atoms by default (see :func:`rmsd_atoms`).
    reference : int, optional
        Frame the others are compared with.

    Returns
    -------
    numpy.ndarray
        RMSD in angstrom, one value per frame (0 for the reference).
    """
    x = np.asarray(frames, dtype=float)
    if atoms is not None:
        x = x[:, list(atoms)]
    n = x.shape[1]
    if n == 0:
        return np.zeros(len(x))
    moving = x - x.mean(axis=1, keepdims=True)
    fixed = moving[reference]
    covariance = np.einsum("fni,nj->fij", moving, fixed)
    u, s, vt = np.linalg.svd(covariance)
    # a reflection would fit better than any rotation: use the rotation, flipping the smallest axis
    sign = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", u, vt)))
    sign[sign == 0] = 1.0
    overlap = s[:, 0] + s[:, 1] + sign * s[:, 2]
    squared = ((moving ** 2).sum(axis=(1, 2)) + (fixed ** 2).sum() - 2.0 * overlap) / n
    return np.sqrt(np.clip(squared, 0.0, None))


def align_frames(frames: ArrayLike, atoms: Sequence[int] | None = None, reference: int = 0) -> np.ndarray:
    """Move every frame onto a reference frame by its best rotation and translation.

    This removes the overall drift and tumbling of a trajectory, so only
    changes of shape remain. Each frame is fitted with the Kabsch algorithm
    on the chosen atoms and then moved as a whole, so its internal geometry
    is unchanged.

    Parameters
    ----------
    frames : array_like
        Coordinates, shape ``(n_frames, n_atoms, 3)``.
    atoms : sequence of int, optional
        Atoms the fit uses (at least three, not all on a line); all atoms by
        default (see :func:`rmsd_atoms`).
    reference : int, optional
        Frame the others are moved onto; it stays where it is.

    Returns
    -------
    numpy.ndarray
        The moved frames, same shape as ``frames``.
    """
    x = np.asarray(frames, dtype=float)
    fit = x if atoms is None else x[:, list(atoms)]
    if not fit.shape[1]:
        return x.copy()
    centers = fit.mean(axis=1)
    moving = fit - centers[:, None]
    covariance = np.einsum("fni,nj->fij", moving, moving[reference])
    u, _, vt = np.linalg.svd(covariance)
    v = np.transpose(vt, (0, 2, 1))
    ut = np.transpose(u, (0, 2, 1))
    sign = np.sign(np.linalg.det(v @ ut))
    sign[sign == 0] = 1.0
    v[:, :, 2] *= sign[:, None]   # a proper rotation, never a reflection
    rotations = v @ ut
    return np.einsum("fij,fnj->fni", rotations, x - centers[:, None]) + centers[reference]
