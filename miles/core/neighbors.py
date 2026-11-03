"""Fast neighbour searches on a uniform grid.

Points are binned into cubic cells as wide as the cutoff, so every pair
within the cutoff lies in the same or adjacent cells. Candidates come from
the 27 surrounding cells and are filtered by distance; the searches are
vectorised over all points for each of the 27 cell offsets.

Every search here includes the cutoff itself: two points exactly
``cutoff`` apart are within it, as "within 5 angstrom" reads.
"""

import itertools

import numpy as np
from numpy.typing import ArrayLike

OFFSETS = np.array(sorted(itertools.product((-1, 0, 1), repeat=3), key=lambda o: sum(map(abs, o))),
                   dtype=np.int64)
"""The 27 cell offsets around a cell, the cell itself first."""

PAIR_LIMIT = 2_000_000
"""Candidate pairs looked at in one go by :func:`near_mask`, which bounds its memory."""

_NO_PAIRS = (np.zeros(0, np.int64), np.zeros(0, np.int64))


def _cells(points: np.ndarray, origin: np.ndarray, size: float) -> np.ndarray:
    return np.floor((points - origin) / size).astype(np.int64) + 1


def _key(cells: np.ndarray, dims: np.ndarray) -> np.ndarray:
    return (cells[:, 0] * dims[1] + cells[:, 1]) * dims[2] + cells[:, 2]


def _expand(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Row of each candidate and its rank among the candidates of that row."""
    rows = np.repeat(np.arange(len(counts)), counts)
    rank = np.arange(int(counts.sum())) - np.repeat(np.cumsum(counts) - counts, counts)
    return rows, rank


def cross_pairs(a: ArrayLike, b: ArrayLike, cutoff: float) -> tuple[np.ndarray, np.ndarray]:
    """Find all pairs of points from two sets within a cutoff of each other.

    Parameters
    ----------
    a, b : array_like
        Point sets, shape ``(n, 3)`` and ``(m, 3)``.
    cutoff : float
        Distance limit (inclusive).

    Returns
    -------
    i, j : numpy.ndarray
        Indices into ``a`` and ``b`` of each close pair.
    """
    pa = np.asarray(a, dtype=float).reshape(-1, 3)
    pb = np.asarray(b, dtype=float).reshape(-1, 3)
    if not len(pa) or not len(pb) or cutoff <= 0:
        return _NO_PAIRS
    origin = np.minimum(pa.min(axis=0), pb.min(axis=0))
    cells_a, cells_b = _cells(pa, origin, cutoff), _cells(pb, origin, cutoff)
    dims = np.maximum(cells_a.max(axis=0), cells_b.max(axis=0)) + 2
    order = np.argsort(_key(cells_b, dims), kind="stable")
    sorted_keys = _key(cells_b, dims)[order]
    found_i, found_j = [], []
    for offset in OFFSETS:
        keys = _key(cells_a + offset, dims)
        lo = np.searchsorted(sorted_keys, keys, "left")
        counts = np.searchsorted(sorted_keys, keys, "right") - lo
        if not counts.any():
            continue
        i, rank = _expand(counts)
        j = order[lo[i] + rank]
        d = pa[i] - pb[j]
        close = np.einsum("ij,ij->i", d, d) <= cutoff * cutoff
        found_i.append(i[close])
        found_j.append(j[close])
    if not found_i:
        return _NO_PAIRS
    return np.concatenate(found_i), np.concatenate(found_j)


def pairs_within(points: ArrayLike, cutoff: float) -> np.ndarray:
    """Find all pairs ``i < j`` of points within a cutoff of each other.

    Parameters
    ----------
    points : array_like
        Coordinates, shape ``(n, 3)``.
    cutoff : float
        Distance limit (inclusive).

    Returns
    -------
    numpy.ndarray
        Index pairs, shape ``(k, 2)``.
    """
    i, j = cross_pairs(points, points, cutoff)
    keep = i < j
    return np.stack([i[keep], j[keep]], axis=1)


def _batches(counts: np.ndarray) -> list[np.ndarray]:
    """Split the rows with candidates into batches of at most about :data:`PAIR_LIMIT` candidate pairs."""
    rows = np.flatnonzero(counts)
    if not len(rows):
        return []
    totals = np.cumsum(counts[rows])
    cuts = np.searchsorted(totals, np.arange(PAIR_LIMIT, int(totals[-1]), PAIR_LIMIT), "right")
    return [batch for batch in np.split(rows, cuts) if len(batch)]


def near_mask(points: ArrayLike, centers: ArrayLike, cutoff: float) -> np.ndarray:
    """Mark the points within a cutoff of any centre.

    The smaller set is walked cell by cell and the larger one looked up in
    sorted cell order. Candidate pairs are looked at in batches, so memory
    stays bounded even when both sets are large; when the points are
    walked, those found are dropped from the remaining cells.

    Parameters
    ----------
    points, centers : array_like
        Coordinates, shape ``(n, 3)`` and ``(m, 3)``.
    cutoff : float
        Distance limit (inclusive).

    Returns
    -------
    numpy.ndarray
        Boolean mask over ``points``.
    """
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    centers = np.asarray(centers, dtype=float).reshape(-1, 3)
    found = np.zeros(len(points), dtype=bool)
    if not len(points) or not len(centers) or cutoff <= 0:
        return found
    origin = np.minimum(points.min(axis=0), centers.min(axis=0))
    cell_points, cell_centers = _cells(points, origin, cutoff), _cells(centers, origin, cutoff)
    dims = np.maximum(cell_points.max(axis=0), cell_centers.max(axis=0)) + 2
    limit = cutoff * cutoff
    if len(centers) < len(points):
        # few centres: walk them and look the points up
        point_keys = _key(cell_points, dims)
        order = np.argsort(point_keys, kind="stable")
        sorted_keys, sorted_points = point_keys[order], points[order]
        for offset in OFFSETS:
            keys = _key(cell_centers + offset, dims)
            low = np.searchsorted(sorted_keys, keys, "left")
            counts = np.searchsorted(sorted_keys, keys, "right") - low
            for batch in _batches(counts):
                c = counts[batch]
                row = np.repeat(batch, c)
                target = low[row] + np.arange(int(c.sum())) - np.repeat(np.cumsum(c) - c, c)
                d = centers[row] - sorted_points[target]
                found[order[target[np.einsum("ij,ij->i", d, d) <= limit]]] = True
        return found
    order = np.argsort(_key(cell_centers, dims), kind="stable")
    sorted_keys = _key(cell_centers, dims)[order]
    sorted_centers = centers[order]
    pending = np.arange(len(points))
    for offset in OFFSETS:
        if not len(pending):
            break
        keys = _key(cell_points[pending] + offset, dims)
        low = np.searchsorted(sorted_keys, keys, "left")
        counts = np.searchsorted(sorted_keys, keys, "right") - low
        for batch in _batches(counts):
            c = counts[batch]
            row = np.repeat(batch, c)
            rank = np.arange(int(c.sum())) - np.repeat(np.cumsum(c) - c, c)
            d = points[pending[row]] - sorted_centers[low[row] + rank]
            close = row[np.einsum("ij,ij->i", d, d) <= limit]
            found[pending[close]] = True
        pending = pending[~found[pending]]
    return found


class PointGrid:
    """Points binned once for many "which points are near here" questions; more can be added.

    Used where one search per atom would otherwise look at every atom, such
    as when hydrogens are added one atom at a time.

    Parameters
    ----------
    points : array_like
        Starting points, shape ``(n, 3)``; they get the indices ``0 .. n - 1``.
    reach : float
        Largest distance asked about (the cell size).
    """

    def __init__(self, points: ArrayLike, reach: float) -> None:
        self.points = np.asarray(points, dtype=float).reshape(-1, 3)
        self.reach = float(reach)
        self._base = len(self.points)
        self._origin = self.points.min(axis=0) if self._base else np.zeros(3)
        keys = self._keys(self._cells(self.points)) if self._base else np.zeros(0, np.int64)
        self._order = np.argsort(keys, kind="stable")
        self._sorted = keys[self._order]
        self._added = np.zeros((16, 3))   # room for added points, filled up to ``_count``
        self._count = 0
        self._added_cells: dict[int, list[int]] = {}

    def _cells(self, points: np.ndarray) -> np.ndarray:
        return np.floor((points - self._origin) / self.reach).astype(np.int64)

    @staticmethod
    def _keys(cells: np.ndarray) -> np.ndarray:
        shifted = cells + (1 << 20)   # room for cells on either side of the origin
        return (shifted[:, 0] << 42) | (shifted[:, 1] << 21) | shifted[:, 2]

    def __len__(self) -> int:
        return self._base + self._count

    def add(self, point: ArrayLike) -> int:
        """Add a point; returns its index."""
        p = np.asarray(point, dtype=float).reshape(3)
        if self._count == len(self._added):
            self._added = np.concatenate([self._added, np.zeros_like(self._added)])
        self._added[self._count] = p
        self._count += 1
        key = int(self._keys(self._cells(p[None]))[0])
        self._added_cells.setdefault(key, []).append(self._count - 1)
        return self._base + self._count - 1

    def position(self, indices: ArrayLike) -> np.ndarray:
        """Return the coordinates of points by index, shape ``(k, 3)``."""
        idx = np.asarray(indices, dtype=np.intp).reshape(-1)
        out = np.empty((len(idx), 3))
        base = idx < self._base
        out[base] = self.points[idx[base]]
        if (~base).any():
            out[~base] = self._added[idx[~base] - self._base]
        return out

    def near(self, center: ArrayLike, distance: float | None = None) -> np.ndarray:
        """Return the points within ``distance`` (at most the grid's reach) of a point, ascending by index."""
        c = np.asarray(center, dtype=float).reshape(3)
        limit = self.reach if distance is None else min(float(distance), self.reach)
        keys = self._keys(self._cells(c[None]) + OFFSETS)
        found = []
        if self._base:
            lo = np.searchsorted(self._sorted, keys, "left")
            hi = np.searchsorted(self._sorted, keys, "right")
            ranges = [self._order[a:b] for a, b in zip(lo.tolist(), hi.tolist(), strict=True) if b > a]
            if ranges:
                found.append(np.concatenate(ranges))
        if self._added_cells:
            extra = [k for key in keys.tolist() for k in self._added_cells.get(key, ())]
            if extra:
                found.append(np.asarray(extra, dtype=np.intp) + self._base)
        if not found:
            return np.zeros(0, dtype=np.intp)
        idx = np.concatenate(found)
        d = self.position(idx) - c
        return np.sort(idx[np.einsum("ij,ij->i", d, d) <= limit * limit])
