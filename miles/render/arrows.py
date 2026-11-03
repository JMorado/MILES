"""Arrows in 3D, such as the dipole moment or the displacements of a vibration.

An arrow is a cylinder (the shaft) ending in a cone (the head). Arrows are
given as rows of :data:`COLUMNS` numbers: the tail, the tip, the shaft
radius and the colour. :func:`arrow_mesh` turns them into one closed,
shaded triangle mesh for the OpenGL renderer; :func:`arrow_parts` gives the
same arrows as a shaft and a head for the QPainter renderer.
"""

import numpy as np
from numpy.typing import ArrayLike

from miles.render.cartoon import Mesh

COLUMNS = 10
"""Numbers per arrow: tail (3), tip (3), shaft radius, colour (3)."""

SIDES = 18
"""Facets around the shaft and the head."""

HEAD_WIDTH = 2.3
"""Radius of the head as a multiple of the shaft radius."""

HEAD_LENGTH = 5.5
"""Length of the head as a multiple of the shaft radius (at most half the arrow)."""


def arrow_rows(tails: ArrayLike, tips: ArrayLike, radius: float, color: ArrayLike) -> np.ndarray:
    """Rows for arrows of one radius and colour.

    Parameters
    ----------
    tails, tips : array_like
        Start and end of each arrow in angstrom, shape ``(n, 3)``.
    radius : float
        Shaft radius in angstrom.
    color : array_like
        RGB colour, components from 0 to 1.

    Returns
    -------
    numpy.ndarray
        Shape ``(n, COLUMNS)``, float32.
    """
    tails = np.asarray(tails, dtype=float).reshape(-1, 3)
    tips = np.asarray(tips, dtype=float).reshape(-1, 3)
    rows = np.zeros((len(tails), COLUMNS), np.float32)
    rows[:, 0:3], rows[:, 3:6], rows[:, 6] = tails, tips, radius
    rows[:, 7:10] = np.asarray(color, dtype=float)[:3]
    return rows


def arrow_parts(arrows: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split arrows into a shaft and a head.

    Parameters
    ----------
    arrows : numpy.ndarray
        Arrow rows, shape ``(n, COLUMNS)``.

    Returns
    -------
    tuple of numpy.ndarray
        The tail, the base of the head and the tip (each ``(n, 3)``), and the
        head radius (``(n,)``). Arrows too short for their head get a head
        half their length.
    """
    tails, tips = arrows[:, 0:3].astype(float), arrows[:, 3:6].astype(float)
    radius = arrows[:, 6].astype(float)
    axis = tips - tails
    length = np.linalg.norm(axis, axis=1)
    head = np.minimum(HEAD_LENGTH * radius, 0.5 * length)
    unit = axis / np.maximum(length, 1e-12)[:, None]
    bases = tips - unit * head[:, None]
    return tails, bases, tips, HEAD_WIDTH * radius


def _frames(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors perpendicular to each axis and to each other."""
    helper = np.where(np.abs(axis[:, [0]]) < 0.9, [[1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0]])
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u, axis=1)[:, None]
    return u, np.cross(axis, u)


def arrow_mesh(arrows: np.ndarray) -> Mesh | None:
    """Build one closed triangle mesh of arrows, drawn like the other meshes.

    Each arrow has a capped shaft and a cone with a flat base; facets have
    their own normals along the sides, so the shading is smooth around the
    arrow and sharp at its rims.

    Parameters
    ----------
    arrows : numpy.ndarray
        Arrow rows, shape ``(n, COLUMNS)``.

    Returns
    -------
    Mesh or None
        The mesh, or None when there are no arrows of any length.
    """
    arrows = np.asarray(arrows, dtype=np.float32).reshape(-1, COLUMNS)
    keep = np.linalg.norm(arrows[:, 3:6] - arrows[:, 0:3], axis=1) > 1e-6
    arrows = arrows[keep]
    if not len(arrows):
        return None
    tails, bases, tips, head_radius = arrow_parts(arrows)
    radius = arrows[:, 6].astype(float)
    axis = tips - tails
    axis /= np.linalg.norm(axis, axis=1)[:, None]
    u, v = _frames(axis)
    angles = np.linspace(0.0, 2.0 * np.pi, SIDES, endpoint=False)
    ring = np.cos(angles)[None, :, None] * u[:, None, :] + np.sin(angles)[None, :, None] * v[:, None, :]
    n = len(arrows)
    slope = head_radius / np.maximum(np.linalg.norm(tips - bases, axis=1), 1e-12)
    cone_normal = ring + slope[:, None, None] * axis[:, None, :]
    cone_normal /= np.linalg.norm(cone_normal, axis=2)[:, :, None]
    back = np.repeat(-axis[:, None, :], SIDES, axis=1)

    # vertex groups per arrow, each SIDES long (in this order):
    groups = [
        (tails[:, None] + radius[:, None, None] * ring, back),             # 0 tail cap rim
        (tails[:, None] + radius[:, None, None] * ring, ring),             # 1 shaft at the tail
        (bases[:, None] + radius[:, None, None] * ring, ring),             # 2 shaft at the head
        (bases[:, None] + radius[:, None, None] * ring, back),             # 3 inner rim of the head's base
        (bases[:, None] + head_radius[:, None, None] * ring, back),        # 4 outer rim of the head's base
        (bases[:, None] + head_radius[:, None, None] * ring, cone_normal),  # 5 cone at its base
        (np.repeat(tips[:, None], SIDES, axis=1), cone_normal),            # 6 tip, once per facet
    ]
    positions = np.concatenate([g[0] for g in groups] + [tails[:, None]], axis=1)   # + tail cap centre
    normals = np.concatenate([g[1] for g in groups] + [-axis[:, None]], axis=1)
    per = positions.shape[1]

    k = np.arange(SIDES)
    nxt = (k + 1) % SIDES
    tail_center = len(groups) * SIDES

    def band(a: int, b: int) -> np.ndarray:
        """Join two rings with triangles facing away from the axis (or backwards, for the base)."""
        return np.concatenate([np.stack([a * SIDES + k, a * SIDES + nxt, b * SIDES + nxt], axis=1),
                               np.stack([a * SIDES + k, b * SIDES + nxt, b * SIDES + k], axis=1)])

    cap = np.stack([np.full(SIDES, tail_center), nxt, k], axis=1)              # facing backwards
    cone = np.stack([5 * SIDES + k, 5 * SIDES + nxt, 6 * SIDES + k], axis=1)   # one triangle per facet
    triangles = np.concatenate([cap, band(1, 2), band(3, 4), cone])
    indices = (triangles[None, :, :] + (np.arange(n) * per)[:, None, None]).reshape(-1, 3)
    colors = np.repeat(arrows[:, 7:10], per, axis=0).astype(np.float32)
    return Mesh(vertices=positions.reshape(-1, 3).astype(np.float32),
                normals=normals.reshape(-1, 3).astype(np.float32),
                colors=colors, base_colors=colors, indices=indices.astype(np.uint32))
