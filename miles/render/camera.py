"""Orthographic camera, smooth camera moves and framing helpers."""

import math
from dataclasses import dataclass

import numpy as np

from miles.core.geometry import rotation_matrix
from miles.core.molecule import Molecule, Representation
from miles.render.scene import atom_radius

MIN_SCALE = 2.0
MAX_SCALE = 500.0


class Camera:
    """Orthographic camera.

    The scene is rotated about ``center``, scaled to ``scale`` pixels per
    angstrom and shifted by ``pan`` pixels. View space has x to the right,
    y up and z towards the viewer.
    """

    def __init__(self) -> None:
        self.rotation = np.eye(3)
        self.center = np.zeros(3)
        self.scale = 60.0
        self.pan = np.zeros(2)

    def copy(self) -> "Camera":
        """Independent copy."""
        c = Camera()
        c.rotation = self.rotation.copy()
        c.center = self.center.copy()
        c.scale = self.scale
        c.pan = self.pan.copy()
        return c

    def to_view(self, coords: np.ndarray) -> np.ndarray:
        """World coordinates to view coordinates."""
        return (np.asarray(coords, dtype=float) - self.center) @ self.rotation.T

    def to_screen(self, view: np.ndarray, width: float, height: float) -> tuple[np.ndarray, np.ndarray]:
        """View coordinates to pixel coordinates (y down)."""
        sx = width / 2 + self.pan[0] + view[..., 0] * self.scale
        sy = height / 2 + self.pan[1] - view[..., 1] * self.scale
        return sx, sy

    def project(self, point: np.ndarray, width: float, height: float) -> np.ndarray:
        """Pixel position of one world point."""
        sx, sy = self.to_screen(self.to_view(point), width, height)
        return np.array([float(sx), float(sy)])

    def unproject(self, x: float, y: float, width: float, height: float, depth: float = 0.0) -> np.ndarray:
        """World point under a pixel, at a view-space depth."""
        v = np.array([(x - width / 2 - self.pan[0]) / self.scale,
                      -(y - height / 2 - self.pan[1]) / self.scale, depth])
        return self.center + self.rotation.T @ v

    def screen_vector_to_world(self, dx: float, dy: float) -> np.ndarray:
        """World displacement matching a mouse movement in pixels."""
        return self.rotation.T @ np.array([dx / self.scale, -dy / self.scale, 0.0])

    def turn(self, ax: float, ay: float, az: float = 0.0) -> None:
        """Rotate the scene about the screen axes (radians)."""
        R = rotation_matrix([1, 0, 0], ax) @ rotation_matrix([0, 1, 0], ay) @ rotation_matrix([0, 0, 1], az)
        u, _, vt = np.linalg.svd(R @ self.rotation)
        self.rotation = u @ vt

    def zoom_at(self, x: float, y: float, width: float, height: float, factor: float) -> None:
        """Zoom by ``factor`` keeping the point under the cursor fixed."""
        new_scale = max(MIN_SCALE, min(MAX_SCALE, self.scale * factor))
        vx = (x - width / 2 - self.pan[0]) / self.scale
        vy = (y - height / 2 - self.pan[1]) / self.scale
        self.scale = new_scale
        self.pan[0] = x - width / 2 - vx * new_scale
        self.pan[1] = y - height / 2 - vy * new_scale

    def fit(self, coords: np.ndarray, width: float, height: float, fill: float = 0.8,
            pad: float = 1.2, max_scale: float = 200.0) -> None:
        """Centre on a point set and zoom so it fills a fraction of the view."""
        self.pan[:] = 0
        if len(coords) == 0:
            return
        self.center = (coords.min(axis=0) + coords.max(axis=0)) / 2
        radius = float(np.linalg.norm(coords - self.center, axis=1).max()) + pad
        self.scale = max(MIN_SCALE, min(max_scale, fill * min(width, height) / (2 * radius)))

    def view_matrix(self) -> np.ndarray:
        """4x4 world-to-view matrix."""
        m = np.eye(4)
        m[:3, :3] = self.rotation
        m[:3, 3] = -self.rotation @ self.center
        return m

    def projection_matrix(self, width: float, height: float, depth: float) -> np.ndarray:
        """4x4 orthographic projection; view z in ``[-depth, depth]`` maps to the depth range."""
        s = self.scale
        return np.array([
            [2 * s / width, 0, 0, 2 * self.pan[0] / width],
            [0, 2 * s / height, 0, -2 * self.pan[1] / height],
            [0, 0, -1.0 / depth, 0],
            [0, 0, 0, 1],
        ])


@dataclass(frozen=True, slots=True)
class DepthRange:
    """View-space depth of a scene.

    Attributes
    ----------
    near, far : float
        View z of the farthest and nearest points (used for depth fading).
    extent : float
        Half-depth of the projection volume.
    """

    near: float
    far: float
    extent: float


def depth_range(camera: Camera, center: np.ndarray, radius: float) -> DepthRange:
    """Depth range of a bounding sphere seen by a camera."""
    zc = float(camera.to_view(center[None, :])[0, 2])
    return DepthRange(zc - radius, zc + radius, abs(zc) + radius + 10.0)


# -- smooth moves --------------------------------------------------------------------------

def _quaternion(R: np.ndarray) -> np.ndarray:
    w = math.sqrt(max(0.0, 1.0 + R[0, 0] + R[1, 1] + R[2, 2])) / 2
    x = math.copysign(math.sqrt(max(0.0, 1.0 + R[0, 0] - R[1, 1] - R[2, 2])) / 2, R[2, 1] - R[1, 2])
    y = math.copysign(math.sqrt(max(0.0, 1.0 - R[0, 0] + R[1, 1] - R[2, 2])) / 2, R[0, 2] - R[2, 0])
    z = math.copysign(math.sqrt(max(0.0, 1.0 - R[0, 0] - R[1, 1] + R[2, 2])) / 2, R[1, 0] - R[0, 1])
    q = np.array([w, x, y, z])
    return q / np.linalg.norm(q)


def _matrix(q: np.ndarray) -> np.ndarray:
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def _slerp(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    d = float(a @ b)
    if d < 0:
        b, d = -b, -d
    if d > 0.9995:
        q = a + t * (b - a)
        return q / np.linalg.norm(q)
    th = math.acos(d)
    return (math.sin((1 - t) * th) * a + math.sin(t * th) * b) / math.sin(th)


def interpolate(start: Camera, end: Camera, t: float) -> Camera:
    """Camera part of the way between two others.

    Rotations are interpolated along the shortest arc, zoom geometrically.

    Parameters
    ----------
    start, end : Camera
        The two views.
    t : float
        Fraction of the way, from 0 to 1.

    Returns
    -------
    Camera
        The intermediate view.
    """
    cam = Camera()
    cam.rotation = _matrix(_slerp(_quaternion(start.rotation), _quaternion(end.rotation), t))
    cam.center = start.center + (end.center - start.center) * t
    cam.scale = start.scale * (end.scale / start.scale) ** t
    cam.pan = start.pan + (end.pan - start.pan) * t
    return cam


# -- framing --------------------------------------------------------------------------------

def nice_rotation(coords: np.ndarray, tilt: bool = True) -> np.ndarray:
    """View rotation that shows a fragment face-on, optionally tilted a little."""
    if len(coords) < 3:
        R = np.eye(3)
    else:
        x = coords - coords.mean(axis=0)
        _, vecs = np.linalg.eigh(x.T @ x)
        ex, ey = vecs[:, 2], vecs[:, 1]
        R = np.array([ex, ey, np.cross(ex, ey)])
    if tilt:
        R = rotation_matrix([1, 0, 0], math.radians(-20)) @ rotation_matrix([0, 1, 0], math.radians(24)) @ R
    return R


def thumbnail_camera(mol: Molecule, width: float, height: float, mode: Representation = Representation.BALL,
                     rotation: np.ndarray | None = None, fill: float = 0.9) -> Camera:
    """Camera that frames a whole molecule in a small picture.

    Parameters
    ----------
    mol : Molecule
        Structure to frame.
    width, height : float
        Picture size in pixels.
    mode : Representation, optional
        Representation, which sets the margin for the atom radii.
    rotation : numpy.ndarray, optional
        View rotation (default: :func:`nice_rotation`).
    fill : float, optional
        Fraction of the picture the molecule may cover.

    Returns
    -------
    Camera
        The camera.
    """
    cam = Camera()
    cam.rotation = nice_rotation(mol.coords) if rotation is None else rotation
    if len(mol):
        view = cam.to_view(mol.coords)
        lo, hi = view.min(axis=0), view.max(axis=0)
        cam.center = cam.rotation.T @ ((lo + hi) / 2)
        pad = 2 * max(atom_radius(symbol, mode) for symbol in set(mol.symbols)) + 0.25
        cam.scale = min(fill * width / max(hi[0] - lo[0] + pad, 1.0),
                        fill * height / max(hi[1] - lo[1] + pad, 1.0), 70.0)
    return cam
