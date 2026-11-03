"""QPainter drawing of scenes: the fallback renderer and thumbnails.

Spheres become shaded discs, each bond half a shaded quad and each triangle
of a surface or isosurface a flat-shaded polygon, all sorted back to front.
The look follows the OpenGL renderer without its exact intersections.
"""

import math
from collections.abc import Callable, Sequence
from typing import NamedTuple

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import (
    QBrush,
    QColor,
    QLinearGradient,
    QPainter,
    QPainterPath,
    QPen,
    QPixmap,
    QPolygonF,
    QRadialGradient,
    QTransform,
)

from miles.core.molecule import Molecule, Representation
from miles.render.arrows import arrow_parts
from miles.render.camera import Camera, thumbnail_camera
from miles.render.cartoon import Mesh
from miles.render.scene import CylinderKind, Highlight, Scene, Style, build_scene, trace_cylinders
from miles.render.theme import LIGHT, Theme, rgb
from miles.render.triangles import ScreenTriangles, screen_triangles

FADE = 0.28
OUTLINE_PX = 1.15
HALO_PX = 3.2
WHITE, BLACK = QColor(255, 255, 255), QColor(0, 0, 0)


class _Piece(NamedTuple):
    """Something to draw at a depth: a sphere row, half a cylinder row or an arrow head (view space)."""

    depth: float
    is_sphere: bool
    row: int
    half: int
    start: np.ndarray
    end: np.ndarray
    head: int = -1       # arrow row for an arrow head (start: centre of its base, end: tip)


Shader = Callable[[np.ndarray, float], QColor]


def _color(values: np.ndarray) -> QColor:
    return QColor.fromRgbF(float(values[0]), float(values[1]), float(values[2]))


def _mix(a: QColor, b: QColor, t: float) -> QColor:
    t = max(0.0, min(1.0, t))
    return QColor.fromRgbF(a.redF() + (b.redF() - a.redF()) * t,
                           a.greenF() + (b.greenF() - a.greenF()) * t,
                           a.blueF() + (b.blueF() - a.blueF()) * t)


def _with_trace(scene: Scene) -> Scene:
    """Scene with the cartoon replaced by a thick backbone trace."""
    if scene.mesh is None or scene.mesh.empty:
        return scene
    trace = trace_cylinders(scene.mesh)
    return Scene(spheres=scene.spheres, sphere_atoms=scene.sphere_atoms,
                 cylinders=np.concatenate([scene.cylinders, trace]).astype(np.float32),
                 cylinder_trim=np.concatenate([scene.cylinder_trim, np.zeros((len(trace), 2), np.float32)]))


def _with_arrow_shafts(scene: Scene) -> Scene:
    """Scene with the shafts of its arrows added as cylinders; the heads are drawn as cones."""
    if not len(scene.arrows):
        return scene
    tails, bases, _, _ = arrow_parts(scene.arrows)
    rows = np.zeros((len(tails), 14), np.float32)
    rows[:, 0:3], rows[:, 3:6], rows[:, 6] = tails, bases, scene.arrows[:, 6]
    rows[:, 7:10] = rows[:, 10:13] = scene.arrows[:, 7:10]
    rows[:, 13] = CylinderKind.SOLID
    return Scene(spheres=scene.spheres, sphere_atoms=scene.sphere_atoms,
                 cylinders=np.concatenate([scene.cylinders, rows]).astype(np.float32),
                 cylinder_trim=np.concatenate([scene.cylinder_trim, np.zeros((len(rows), 2), np.float32)]),
                 mesh=scene.mesh, arrows=scene.arrows)


def _pieces(scene: Scene, camera: Camera, trim_bonds: bool) -> list[_Piece]:
    """Spheres, bond halves and arrow heads sorted back to front."""
    S, C = scene.spheres, scene.cylinders
    centers = camera.to_view(S[:, 0:3])
    pieces = [_Piece(float(v[2]), True, k, 0, v, v) for k, v in enumerate(centers)]
    if len(C):
        A, B = camera.to_view(C[:, 0:3]), camera.to_view(C[:, 3:6])
        for k in range(len(C)):
            for half, (a, b) in enumerate(((A[k], B[k]), (B[k], A[k]))):
                length = float(np.linalg.norm(b - a))
                if length < 1e-6:
                    continue
                trim = 0.0
                if trim_bonds:   # hide the part of the bond inside the ball
                    r = float(scene.cylinder_trim[k, half])
                    trim = math.sqrt(max(r * r - C[k, 6] ** 2, 0.0)) * 0.97
                if trim >= length / 2:
                    continue
                start = a + (b - a) / length * trim
                mid = (a + b) / 2
                pieces.append(_Piece(float(start[2] + mid[2]) / 2, False, k, half, start, mid))
    if len(scene.arrows):
        _, bases, tips, _ = arrow_parts(scene.arrows)
        for k, (base, tip) in enumerate(zip(camera.to_view(bases), camera.to_view(tips), strict=True)):
            pieces.append(_Piece(float(base[2] + tip[2]) / 2, False, -1, 0, base, tip, k))
    pieces.sort(key=lambda piece: (piece.depth, piece.is_sphere))
    return pieces


def paint_scene(p: QPainter, width: float, height: float, scene: Scene, camera: Camera, theme: Theme,
                style: Style, background: bool = True,
                meshes: Sequence[tuple[Mesh | None, float]] = ()) -> None:
    """Draw a scene with QPainter.

    Parameters
    ----------
    p : QPainter
        Active painter.
    width, height : float
        Size of the area in pixels.
    scene : Scene
        What to draw; a cartoon is drawn as a thick backbone trace.
    camera : Camera
        View.
    theme : Theme
        Colours.
    style : Style
        Display options (depth fading, ball-and-stick bond trimming).
    background : bool, optional
        Fill the area with the theme background first.
    meshes : sequence of (Mesh or None, float), optional
        Closed surfaces to draw with their opacity: the molecular surface and
        the isosurface of volumetric data. Large ones are simplified first
        (:mod:`miles.render.triangles`).
    """
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    if background:
        p.fillRect(QRectF(0, 0, width, height), theme.background)
    layers = [screen_triangles(mesh, camera, width, height, alpha) for mesh, alpha in meshes
              if mesh is not None and not mesh.empty]
    layers = [layer for layer in layers if len(layer.depth)]
    if scene.empty and not layers:
        return
    scene = _with_arrow_shafts(_with_trace(scene))
    pieces = _pieces(scene, camera, style.mode == Representation.BALL)
    if not pieces and not layers:
        return
    S, C = scene.spheres, scene.cylinders
    sx, sy = camera.to_screen(camera.to_view(S[:, 0:3]), width, height)
    z = np.concatenate([camera.to_view(S[:, 0:3])[:, 2], camera.to_view(C[:, 0:3])[:, 2],
                        [piece.depth for piece in pieces], *(layer.depth for layer in layers)])
    z_front = float(z.max())
    span = max(z_front - float(z.min()), 8.0)
    fade = FADE if style.depth_fade else 0.0
    highlight_colors = theme.highlight_colors()
    outline = QPen(theme.outline, OUTLINE_PX)
    outline.setJoinStyle(Qt.PenJoinStyle.RoundJoin)

    def shaded(values: np.ndarray, z: float) -> QColor:
        c = _color(values)
        return _mix(c, theme.background, fade * (z_front - z) / span) if fade else c

    head_radius = arrow_parts(scene.arrows)[3] if len(scene.arrows) else np.zeros(0)
    tri_x, tri_y, tri_colors, tri_depth = _triangle_lists(layers, theme, fade, z_front, span)
    keys = np.concatenate([[piece.depth for piece in pieces], tri_depth])
    smooth = True
    for item in np.argsort(keys, kind="stable").tolist():   # stable: pieces keep their own order on ties
        if item >= len(pieces):   # a triangle of a surface
            if smooth:
                p.setRenderHint(QPainter.RenderHint.Antialiasing, False)   # no seams between triangles
                p.setPen(Qt.PenStyle.NoPen)
                smooth = False
            j = item - len(pieces)
            x, y = tri_x[j], tri_y[j]
            p.setBrush(QColor(*tri_colors[j]))
            p.drawPolygon(QPolygonF([QPointF(x[0], y[0]), QPointF(x[1], y[1]), QPointF(x[2], y[2])]))
            continue
        if not smooth:
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            smooth = True
        piece = pieces[item]
        k = piece.row
        if piece.head >= 0:
            arrow = scene.arrows[piece.head]
            _paint_arrow_head(p, piece, float(head_radius[piece.head]), shaded(arrow[7:10], piece.depth),
                              camera, width, height, outline)
        elif piece.is_sphere:
            center = QPointF(sx[k], sy[k])
            r = -S[k, 3] if S[k, 3] < 0 else S[k, 3] * camera.scale
            c = shaded(S[k, 4:7], piece.depth)
            highlight = int(S[k, 7])
            if highlight != Highlight.NONE:
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(highlight_colors[highlight])
                p.drawEllipse(center, r + HALO_PX, r + HALO_PX)
            g = QRadialGradient(QPointF(sx[k] - 0.3 * r, sy[k] - 0.35 * r), 1.25 * r)
            g.setColorAt(0.0, _mix(c, WHITE, 0.22))
            g.setColorAt(0.6, c)
            g.setColorAt(1.0, _mix(c, BLACK, 0.22))
            p.setBrush(QBrush(g))
            p.setPen(Qt.PenStyle.NoPen if S[k, 3] < 0 else outline)
            p.drawEllipse(center, r, r)
        else:
            _paint_bond_half(p, C[k], piece, camera, width, height, shaded, outline)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    p.setPen(Qt.PenStyle.NoPen)


def _triangle_lists(layers: list[ScreenTriangles], theme: Theme, fade: float, z_front: float,
                    span: float) -> tuple[list[list[float]], list[list[float]], list[list[int]], np.ndarray]:
    """Corners, RGBA colours (0 to 255, depth-faded) and depths of all surface triangles."""
    if not layers:
        return [], [], [], np.zeros(0)
    x = np.concatenate([layer.x for layer in layers])
    y = np.concatenate([layer.y for layer in layers])
    depth = np.concatenate([layer.depth for layer in layers])
    color = np.concatenate([layer.rgb for layer in layers])
    alpha = np.concatenate([np.full(len(layer.depth), layer.alpha) for layer in layers])
    if fade:
        t = np.clip(fade * (z_front - depth) / span, 0.0, 1.0)[:, None]
        color = color + (np.array(rgb(theme.background)) - color) * t
    rgba = np.rint(np.column_stack([color, alpha]) * 255).astype(int)
    return x.tolist(), y.tolist(), rgba.tolist(), depth


def _paint_arrow_head(p: QPainter, piece: _Piece, radius: float, color: QColor, camera: Camera,
                      width: float, height: float, outline: QPen) -> None:
    """Draw a cone as seen from the camera: its outline triangle and its base, in depth order."""
    (bx, tx), (by, ty) = camera.to_screen(np.array([piece.start, piece.end]), width, height)
    dx, dy = tx - bx, ty - by
    span = math.hypot(dx, dy)
    r = radius * camera.scale
    length = float(np.linalg.norm(piece.end - piece.start)) * camera.scale
    ux, uy = (dx / span, dy / span) if span > 1e-6 else (1.0, 0.0)
    squash = math.sqrt(max(0.0, 1.0 - min(1.0, span / max(length, 1e-6)) ** 2))   # base seen edge-on: 0
    base = QPainterPath()
    base.addEllipse(QPointF(0, 0), r * squash, r)
    angle = math.degrees(math.atan2(uy, ux))
    ellipse = QTransform().translate(bx, by).rotate(angle).map(base)
    left, right = QPointF(bx - uy * r, by + ux * r), QPointF(bx + uy * r, by - ux * r)
    cone = QPolygonF([left, QPointF(tx, ty), right])
    g = QLinearGradient(left, right)
    g.setColorAt(0.0, _mix(color, WHITE, 0.18))
    g.setColorAt(1.0, _mix(color, BLACK, 0.22))
    facing = piece.start[2] > piece.end[2]   # the base is nearer the viewer than the tip
    p.setPen(outline)
    if facing:
        p.setBrush(QBrush(g))
        p.drawPolygon(cone)
        p.setBrush(_mix(color, BLACK, 0.12))
        p.drawPath(ellipse)
    else:
        p.setBrush(_mix(color, BLACK, 0.12))
        p.drawPath(ellipse)
        p.setBrush(QBrush(g))
        p.drawPolygon(cone)


def _paint_bond_half(p: QPainter, row: np.ndarray, piece: _Piece, camera: Camera, width: float, height: float,
                     shaded: Shader, outline: QPen) -> None:
    x1, y1 = camera.to_screen(piece.start, width, height)
    x2, y2 = camera.to_screen(piece.end, width, height)
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length < 0.3:
        return
    r = -float(row[6]) if row[6] < 0 else float(row[6]) * camera.scale
    nx, ny = -dy / length * r, dx / length * r
    c = shaded(row[7:10] if piece.half == 0 else row[10:13], piece.depth)
    kind = int(row[13])
    dashes = [(0.0, 0.35), (0.6, 0.95)] if kind == CylinderKind.DASHED else [(0.0, 1.0)]
    for t0, t1 in dashes:
        ax, ay = x1 + dx * t0, y1 + dy * t0
        bx, by = x1 + dx * t1, y1 + dy * t1
        left_a, left_b = QPointF(ax + nx, ay + ny), QPointF(bx + nx, by + ny)
        right_a, right_b = QPointF(ax - nx, ay - ny), QPointF(bx - nx, by - ny)
        g = QLinearGradient(left_a, right_a)
        g.setColorAt(0.0, _mix(c, WHITE, 0.12))
        g.setColorAt(1.0, _mix(c, BLACK, 0.18))
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QBrush(g))
        p.drawPolygon(QPolygonF([left_a, left_b, right_b, right_a]))
        if row[6] > 0 and kind == CylinderKind.SOLID:
            p.setPen(outline)
            p.drawLine(left_a, left_b)
            p.drawLine(right_a, right_b)


def render_pixmap(mol: Molecule, width: int, height: int, theme: Theme = LIGHT, style: Style | None = None,
                  rotation: np.ndarray | None = None, highlights: dict[int, Highlight] | None = None,
                  ratio: float = 2.0, fill: float = 0.9, background: bool = False) -> QPixmap:
    """Picture of a whole molecule, e.g. a fragment thumbnail.

    Parameters
    ----------
    mol : Molecule
        Structure to draw.
    width, height : int
        Size in logical pixels.
    theme : Theme, optional
        Colours.
    style : Style, optional
        Display options.
    rotation : numpy.ndarray, optional
        View rotation (default: face-on).
    highlights : dict of int to Highlight, optional
        Highlighted atoms.
    ratio : float, optional
        Device pixel ratio of the pixmap.
    fill : float, optional
        Fraction of the picture the molecule may cover.
    background : bool, optional
        Fill with the theme background instead of transparency.

    Returns
    -------
    QPixmap
        The picture.
    """
    style = style or Style()
    pm = QPixmap(int(width * ratio), int(height * ratio))
    pm.setDevicePixelRatio(ratio)
    pm.fill(Qt.GlobalColor.transparent)
    camera = thumbnail_camera(mol, width, height, style.mode, rotation, fill)
    p = QPainter(pm)
    scene = build_scene(mol, style, theme, highlights)
    paint_scene(p, width, height, scene, camera, theme, style, background)
    p.end()
    return pm
