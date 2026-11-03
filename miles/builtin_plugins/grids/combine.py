"""Combining two grids point by point: their sum or difference (Combine... in the grid panel)."""

from enum import StrEnum

import numpy as np

from miles import Volume


class Combination(StrEnum):
    """How :func:`combine_volumes` joins two grids."""

    DIFFERENCE = "difference"   # A - B
    SUM = "sum"                 # A + B


def same_points(a: Volume, b: Volume, tolerance: float = 1e-4) -> bool:
    """Whether two grids have the same points (shape, origin and step vectors).

    Parameters
    ----------
    a, b : Volume
        The grids.
    tolerance : float, optional
        Largest difference in angstrom that still counts as the same.

    Returns
    -------
    bool
        True when the values can be combined point by point.
    """
    return (a.values.shape == b.values.shape
            and bool(np.allclose(a.origin, b.origin, atol=tolerance))
            and bool(np.allclose(a.axes, b.axes, atol=tolerance)))


def combine_volumes(a: Volume, b: Volume, combination: Combination = Combination.DIFFERENCE,
                    name: str = "") -> Volume:
    """Sum or difference of two grids, on the points of the first.

    Parameters
    ----------
    a, b : Volume
        The grids. When ``b`` has other points it is interpolated at the
        points of ``a`` (:meth:`Volume.sample`); where ``a`` reaches beyond
        ``b``, the values at the edge of ``b`` are used.
    combination : Combination, optional
        ``a - b`` or ``a + b``.
    name : str, optional
        Name of the result (default: e.g. ``"A - B"`` with the grids' names).

    Returns
    -------
    Volume
        The result. A difference is marked as taking both signs
        (:attr:`Volume.both_signs`), so it is drawn with a positive and a
        negative lobe even where one sign is small.
    """
    if same_points(a, b):
        other = np.asarray(b.values, dtype=np.float64)
    else:
        nx, ny, nz = a.values.shape
        j, k = np.meshgrid(np.arange(ny), np.arange(nz), indexing="ij")
        plane = np.stack([j.ravel(), k.ravel()], axis=1) @ np.asarray(a.axes, dtype=float)[1:3]
        other = np.empty(a.values.shape)
        for i in range(nx):   # one plane at a time keeps memory small for large grids
            points = np.asarray(a.origin, dtype=float) + i * np.asarray(a.axes, dtype=float)[0] + plane
            other[i] = b.sample(points).reshape(ny, nz)
    sign = -1.0 if combination == Combination.DIFFERENCE else 1.0
    values = np.asarray(a.values, dtype=np.float64) + sign * other
    symbol = "-" if combination == Combination.DIFFERENCE else "+"
    both = True if combination == Combination.DIFFERENCE else None
    return Volume(name or f"{a.name} {symbol} {b.name}", np.array(a.origin, dtype=float),
                  np.array(a.axes, dtype=float), values.astype(np.float32), both)
