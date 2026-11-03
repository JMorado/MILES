"""Trajectories as the document keeps them: frames that are never changed in place.

A :class:`Trajectory` holds one read-only coordinate array per frame,
with optional text and energy for each frame. Changing a frame makes a new
trajectory that shares every other frame with the old one, so an undo step
that keeps the old trajectory costs one frame, not the whole trajectory.

The document shows one frame at a time in its molecule; see
:mod:`miles.core.document` for how edits, playback and undo use them.
"""

from collections.abc import Sequence
from dataclasses import dataclass, replace
from functools import cached_property

import numpy as np
from numpy.typing import ArrayLike


@dataclass(frozen=True, eq=False)
class Trajectory:
    """Frames of a structure: coordinates, labels and energies.

    Attributes
    ----------
    frames : tuple of numpy.ndarray
        Coordinates of each frame, shape ``(n_atoms, 3)``, read-only.
    labels : tuple of str or None
        Text for each frame, e.g. conformer energies or scan values.
    energies : tuple of float or None
        Energy of each frame in kcal/mol; None unless every frame has one.
    loop : bool
        The motion is periodic (a vibration), so playback may blend the last
        frame into the first.
    """

    frames: tuple[np.ndarray, ...]
    labels: tuple[str, ...] | None = None
    energies: tuple[float, ...] | None = None
    loop: bool = False

    @classmethod
    def from_array(cls, frames: ArrayLike, labels: Sequence[str] | None = None,
                   energies: Sequence[float] | None = None, loop: bool = False) -> "Trajectory | None":
        """Build a trajectory from coordinates of shape ``(n_frames, n_atoms, 3)``.

        The coordinates are copied, so the caller's array stays theirs.
        Energies are kept only when there is one for every frame.

        Returns
        -------
        Trajectory or None
            None for an empty array: no frames means no trajectory.

        Raises
        ------
        ValueError
            If the array does not have the shape of a trajectory.
        """
        array = np.array(frames, dtype=np.float64)
        if array.ndim != 3 or array.shape[2] != 3:
            raise ValueError(f"Frames must have the shape (n_frames, n_atoms, 3), not {array.shape}")
        if not len(array):
            return None
        array.flags.writeable = False
        kept = None
        if energies is not None and len(energies) == len(array):
            kept = tuple(float(e) for e in energies)
        return cls(tuple(array[k] for k in range(len(array))),
                   None if labels is None else tuple(str(label) for label in labels), kept, bool(loop))

    def __len__(self) -> int:
        return len(self.frames)

    @property
    def atom_count(self) -> int:
        """Number of atoms in each frame."""
        return int(self.frames[0].shape[0])

    @cached_property
    def array(self) -> np.ndarray:
        """All frames as one read-only array, shape ``(n_frames, n_atoms, 3)``, made once when asked for."""
        stacked = np.stack(self.frames)
        stacked.flags.writeable = False
        return stacked

    def coords(self, k: int) -> np.ndarray:
        """Coordinates of frame ``k`` (counted from 0, wrapping around)."""
        return self.frames[int(k) % len(self.frames)]

    def with_frame(self, k: int, coords: ArrayLike) -> "Trajectory":
        """Return a trajectory whose frame ``k`` holds other coordinates; the other frames are shared."""
        new = np.array(coords, dtype=np.float64)
        if new.shape != self.frames[0].shape:
            raise ValueError(f"Frame coordinates must have the shape {self.frames[0].shape}, not {new.shape}")
        new.flags.writeable = False
        k = int(k) % len(self.frames)
        return replace(self, frames=(*self.frames[:k], new, *self.frames[k + 1:]))

    def blend(self, position: float) -> tuple[int, float, np.ndarray]:
        """Coordinates at a point between two frames, for smooth playback.

        ``position`` counts frames from 0; its whole part is the frame
        before. Between the last frame and the first, coordinates are blended
        only for periodic motion; otherwise the last frame is held.

        Returns
        -------
        tuple of (int, float, numpy.ndarray)
            The frame before, how far past it (0 to 1) and the coordinates.
        """
        count = len(self.frames)
        position %= count
        k = int(position)
        t = position - k
        following = k + 1 if k + 1 < count else (0 if self.loop else k)
        return k, t, self.frames[k] * (1.0 - t) + self.frames[following] * t
