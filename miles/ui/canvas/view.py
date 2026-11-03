"""Moving the camera of a canvas: fitting, centring, zooming, and smooth moves between views."""

import time
from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
from PySide6.QtCore import QPointF, QTimer

from miles.render.camera import Camera, interpolate

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

ANIMATION_SECONDS = 0.24

ANIMATION_STEP_MS = 16

FIT_FILL = 0.72

FIT_MAX_SCALE = 95.0

AXIS_VIEWS = {
    "x": np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float),
    "y": np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float),
    "z": np.eye(3),
}
"""View rotations that look down each axis."""


class ViewController:
    """The camera of one canvas and the smooth moves between views.

    It reads the canvas's size, the insets its panels cover and the
    document's coordinates, and asks the canvas to redraw; it changes
    nothing else.

    Parameters
    ----------
    canvas : Canvas
        The canvas whose view this is.

    Attributes
    ----------
    camera : Camera
        The view shown now.
    needs_fit : bool
        The structure is to be fitted as soon as the canvas has a size.
    """

    def __init__(self, canvas: "Canvas") -> None:
        self.canvas = canvas
        self.camera = Camera()
        self.needs_fit = len(canvas.doc.mol) > 0
        self._animation: tuple[Camera, Camera, float] | None = None
        self._timer = QTimer(canvas)
        self._timer.setInterval(ANIMATION_STEP_MS)
        self._timer.timeout.connect(self._step)

    @property
    def animating(self) -> bool:
        """Whether the camera is moving to a new view."""
        return self._animation is not None

    def stop(self) -> None:
        """Stay where the camera is now (a click ends a smooth move)."""
        self._animation = None
        self._timer.stop()

    def visible_coords(self) -> np.ndarray:
        """Coordinates of the atoms not hidden by focus mode (all atoms if every one is hidden)."""
        doc = self.canvas.doc
        shown = ~doc.display.hidden
        return doc.mol.coords[shown] if shown.any() else doc.mol.coords

    def usable_insets(self) -> tuple[int, int]:
        """Insets to keep clear when fitting, dropped when they would leave too little room."""
        canvas = self.canvas
        top, bottom = canvas.view_insets
        if canvas.height() - top - bottom < max(160, canvas.height() // 3):
            return 0, 0
        return top, bottom

    def fit(self, animate: bool = False) -> None:
        """Zoom to show the whole structure (the visible part in focus mode)."""
        canvas = self.canvas
        if canvas.width() < 10:
            self.needs_fit = True
            return
        target = self.camera.copy()
        top, bottom = self.usable_insets()
        target.fit(self.visible_coords(), canvas.width(), canvas.height() - top - bottom, fill=FIT_FILL,
                   max_scale=FIT_MAX_SCALE)
        target.pan[1] = (top - bottom) / 2   # centred in the space the panels leave free
        self.needs_fit = False
        self.move_to(target, animate)

    def fit_if_needed(self) -> None:
        """Fit the structure if that waited for the canvas to get a size."""
        if self.needs_fit and len(self.canvas.doc.mol):
            self.fit()

    def center(self) -> None:
        """Centre the structure without changing zoom or orientation."""
        target = self.camera.copy()
        if len(self.canvas.doc.mol):
            c = self.visible_coords()
            target.center = (c.min(axis=0) + c.max(axis=0)) / 2
        top, bottom = self.usable_insets()
        target.pan[:] = (0.0, (top - bottom) / 2)
        self.move_to(target, animate=True)

    def reset(self) -> None:
        """Return to the default orientation and fit."""
        target = self.camera.copy()
        target.rotation = np.eye(3)
        target.fit(self.visible_coords(), self.canvas.width(), self.canvas.height(), fill=FIT_FILL,
                   max_scale=FIT_MAX_SCALE)
        self.move_to(target, animate=True)

    def look_along(self, axis: str) -> None:
        """Look down the x, y or z axis."""
        target = self.camera.copy()
        target.rotation = AXIS_VIEWS[axis]
        self.move_to(target, animate=True)

    def center_on(self, atoms: Iterable[int]) -> None:
        """Rotate about the centre of some atoms from now on."""
        target = self.camera.copy()
        target.center = self.canvas.doc.mol.centroid(atoms)
        target.pan[:] = 0
        self.move_to(target, animate=True)

    def ensure_visible(self, margin: float = 40.0) -> None:
        """Zoom out (never in) if atoms left the view."""
        canvas, mol = self.canvas, self.canvas.doc.mol
        if not len(mol) or canvas.width() < 50:
            return
        w, h = canvas.width(), canvas.height()
        sx, sy = self.camera.to_screen(self.camera.to_view(mol.coords), w, h)
        if sx.min() >= margin and sx.max() <= w - margin and sy.min() >= margin and sy.max() <= h - margin:
            return
        target = self.camera.copy()
        target.fit(mol.coords, w, h, fill=0.75)
        target.scale = min(target.scale, self.camera.scale)
        self.move_to(target, animate=True)

    def zoom_to(self, atoms: Iterable[int]) -> None:
        """Zoom smoothly so that some atoms fill the view."""
        chosen = list(atoms)
        if not chosen:
            return
        target = self.camera.copy()
        target.fit(self.canvas.doc.mol.coords[chosen], self.canvas.width(), self.canvas.height(),
                   fill=FIT_FILL, max_scale=FIT_MAX_SCALE)
        self.move_to(target, animate=True)

    def zoom_at(self, pos: QPointF, factor: float) -> None:
        """Zoom by a factor, keeping the point under the cursor in place."""
        self.camera.zoom_at(pos.x(), pos.y(), self.canvas.width(), self.canvas.height(), factor)
        self.canvas.redraw()

    def move_to(self, target: Camera, animate: bool) -> None:
        """Show another view, moving there smoothly or at once."""
        if animate:
            self._animation = (self.camera.copy(), target, time.monotonic())
            self._timer.start()
        else:
            self.stop()
            self.camera = target
            self.canvas.redraw()

    def _step(self) -> None:
        if self._animation is None:
            self._timer.stop()
            return
        start, end, t0 = self._animation
        t = min(1.0, (time.monotonic() - t0) / ANIMATION_SECONDS)
        self.camera = interpolate(start, end, 1 - (1 - t) ** 3)
        if t >= 1.0:
            self.camera = end
            self.stop()
        self.canvas.redraw()
