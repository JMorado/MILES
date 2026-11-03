"""Triangle meshes for the QPainter renderer: simplification and per-triangle shading.

The QPainter renderer draws surfaces and isosurfaces one triangle at a time,
sorted by depth together with the atoms and bonds. That costs a few
microseconds per triangle, so meshes with more triangles than
:data:`MAX_TRIANGLES` are first simplified by vertex clustering: vertices
that fall into the same cell of a regular grid are merged, and triangles
that collapse are dropped. The result is cached per mesh. Shading follows
the OpenGL renderer's light (:func:`shade`), with one colour per triangle.
"""

from dataclasses import dataclass

import numpy as np

from miles.render.camera import Camera
from miles.render.cartoon import Mesh

MAX_TRIANGLES = 24000
"""Meshes are simplified to at most this many triangles (about half face the viewer)."""

LIGHT = np.array([-0.38, 0.52, 0.76]) / np.linalg.norm([-0.38, 0.52, 0.76])
HALF = (LIGHT + np.array([0.0, 0.0, 1.0])) / np.linalg.norm(LIGHT + np.array([0.0, 0.0, 1.0]))

_cache: dict[tuple[int, int, int], tuple[Mesh, Mesh]] = {}
_CACHE_SIZE = 8


@dataclass
class ScreenTriangles:
    """Triangles of a mesh that face the viewer, ready to draw.

    Attributes
    ----------
    x, y : numpy.ndarray
        Pixel coordinates of the corners, shape ``(n, 3)`` each.
    depth : numpy.ndarray
        View-space depth of each triangle's centre (larger is nearer).
    rgb : numpy.ndarray
        Shaded colour of each triangle, shape ``(n, 3)``, from 0 to 1.
    alpha : float
        Opacity of the mesh.
    """

    x: np.ndarray
    y: np.ndarray
    depth: np.ndarray
    rgb: np.ndarray
    alpha: float


def cluster_vertices(mesh: Mesh, cell: float) -> Mesh:
    """Simplify a mesh by merging the vertices that share a grid cell.

    Parameters
    ----------
    mesh : Mesh
        Triangle mesh.
    cell : float
        Edge of the grid cells in angstrom.

    Returns
    -------
    Mesh
        The simplified mesh: one vertex per occupied cell at the mean
        position, with the mean colour and the normalised sum of the normals;
        triangles whose corners merge are dropped.
    """
    vertices = mesh.vertices.astype(np.float64)
    keys = np.floor((vertices - vertices.min(axis=0)) / cell).astype(np.int64)
    _, cluster, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    cluster = cluster.reshape(-1)
    n = len(counts)

    def mean(values: np.ndarray) -> np.ndarray:
        out = np.zeros((n, values.shape[1]))
        np.add.at(out, cluster, values)
        return out / counts[:, None]

    positions = mean(vertices)
    normals = mean(mesh.normals.astype(np.float64))
    normals /= np.maximum(np.linalg.norm(normals, axis=1), 1e-12)[:, None]
    colors = mean(mesh.colors.astype(np.float64))
    triangles = cluster[mesh.indices.reshape(-1, 3)]
    keep = ((triangles[:, 0] != triangles[:, 1]) & (triangles[:, 1] != triangles[:, 2])
            & (triangles[:, 0] != triangles[:, 2]))
    triangles = triangles[keep]
    _, first = np.unique(np.sort(triangles, axis=1), axis=0, return_index=True)
    triangles = triangles[np.sort(first)]
    return Mesh(vertices=positions.astype(np.float32), normals=normals.astype(np.float32),
                colors=colors.astype(np.float32), base_colors=colors.astype(np.float32),
                indices=triangles.astype(np.uint32).reshape(-1))


def simplified(mesh: Mesh, max_triangles: int = MAX_TRIANGLES) -> Mesh:
    """Return ``mesh``, or a simplified copy with at most ``max_triangles`` triangles (cached).

    The grid cell starts at the mesh's mean edge length and grows until the
    clustered mesh is small enough.
    """
    count = mesh.indices.size // 3
    if count <= max_triangles:
        return mesh
    key = (id(mesh), id(mesh.colors), count)
    hit = _cache.get(key)
    if hit is not None and hit[0] is mesh:
        return hit[1]
    corners = mesh.vertices[mesh.indices.reshape(-1, 3)].astype(np.float64)
    edge = float(np.linalg.norm(corners[:, 1] - corners[:, 0], axis=1).mean())
    cell = edge * np.sqrt(count / max_triangles)
    result = cluster_vertices(mesh, cell)
    for _ in range(12):
        if result.indices.size // 3 <= max_triangles:
            break
        cell *= 1.25
        result = cluster_vertices(mesh, cell)
    if len(_cache) >= _CACHE_SIZE:
        _cache.pop(next(iter(_cache)))
    _cache[key] = (mesh, result)
    return result


def shade(colors: np.ndarray, normals: np.ndarray) -> np.ndarray:
    """Colours lit as the OpenGL renderer lights them, for view-space unit normals.

    Parameters
    ----------
    colors : numpy.ndarray
        Base colours, shape ``(n, 3)``.
    normals : numpy.ndarray
        Unit normals in view space (z towards the viewer), shape ``(n, 3)``.

    Returns
    -------
    numpy.ndarray
        Lit colours, clipped to [0, 1].
    """
    wrap = np.clip(normals @ LIGHT * 0.5 + 0.5, 0.0, 1.0)
    lit = colors * (0.46 + 0.58 * wrap * wrap)[:, None]
    lit = lit + 0.14 * (np.maximum(normals @ HALF, 0.0) ** 38)[:, None]
    facing = np.clip(normals[:, 2] / 0.55, 0.0, 1.0)
    facing = facing * facing * (3 - 2 * facing)   # smoothstep(0, 0.55, n.z)
    lit = lit * (0.84 + 0.16 * facing)[:, None]
    return np.clip(lit, 0.0, 1.0)


def screen_triangles(mesh: Mesh, camera: Camera, width: float, height: float, alpha: float = 1.0,
                     max_triangles: int = MAX_TRIANGLES) -> ScreenTriangles:
    """Project, cull and shade the triangles of a mesh.

    Parameters
    ----------
    mesh : Mesh
        Closed triangle mesh (a surface or an isosurface).
    camera : Camera
        View.
    width, height : float
        Size of the picture in pixels.
    alpha : float, optional
        Opacity of the mesh.
    max_triangles : int, optional
        Simplify larger meshes first (see :func:`simplified`).

    Returns
    -------
    ScreenTriangles
        The triangles facing the viewer and at least partly on the picture.
        Back faces are left out: the meshes are closed, so they are hidden
        behind front faces, as the see-through drawing of the OpenGL
        renderer also leaves them out.
    """
    mesh = simplified(mesh, max_triangles)
    triangles = mesh.indices.reshape(-1, 3)
    view = camera.to_view(mesh.vertices)
    normals = mesh.normals.astype(np.float64) @ camera.rotation.T
    face = normals[triangles].sum(axis=1)
    face /= np.maximum(np.linalg.norm(face, axis=1), 1e-12)[:, None]
    sx, sy = camera.to_screen(view, width, height)
    x, y = sx[triangles], sy[triangles]
    shown = ((face[:, 2] > 0.0) & (x.max(axis=1) >= 0) & (x.min(axis=1) <= width)
             & (y.max(axis=1) >= 0) & (y.min(axis=1) <= height))
    colors = mesh.colors.astype(np.float64)[triangles[shown]].mean(axis=1)
    return ScreenTriangles(x[shown], y[shown], view[triangles[shown], 2].mean(axis=1),
                           shade(colors, face[shown]), alpha)
