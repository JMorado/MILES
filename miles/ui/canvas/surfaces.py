"""What draws the molecule of a canvas: QPainter inside the canvas, or one OpenGL view per window.

A canvas always paints its own overlay (labels, measurements, the status
line). The molecule itself is drawn by its surface:

* :class:`PainterSurface` paints it with QPainter inside the canvas, for
  machines without OpenGL;
* :class:`GLView` is a single OpenGL widget that the stage puts under its
  canvases. It draws the scene of the canvas shown now, so the tabs share
  one OpenGL context, one set of compiled shaders and one framebuffer,
  however many there are; the canvases above it are transparent.
"""

import importlib
import os
from typing import TYPE_CHECKING, Protocol, override

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QImage, QOpenGLContext, QPainter, QPixmap, QSurfaceFormat
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QWidget

from miles.render.arrows import arrow_mesh
from miles.render.cartoon import Mesh
from miles.render.gl import GLRenderer, MeshSlot, Viewport
from miles.render.painter import paint_scene
from miles.render.scene import build_scene
from miles.render.theme import Background, export_theme

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

SOFTWARE_ENV = "MILES_SOFTWARE"
"""Set this environment variable to draw with QPainter even when OpenGL works."""


class Surface(Protocol):
    """What draws the molecule of the canvases of a window."""

    draws_in_canvas: bool
    """True when the canvas paints the molecule itself (:meth:`paint`); False when it is drawn underneath."""

    def paint(self, p: QPainter, canvas: "Canvas") -> None:
        """Paint the molecule of a canvas inside it (only when :attr:`draws_in_canvas`)."""

    def redraw(self, canvas: "Canvas") -> None:
        """Draw a canvas's molecule again soon, if this surface shows it."""

    def render_export(self, canvas: "Canvas", scale: float, background: Background,
                      ambient_occlusion: bool) -> QImage:
        """Picture of a canvas's view for publication, with its overlay."""


class PainterSurface:
    """Draws each canvas's molecule with QPainter, inside the canvas."""

    draws_in_canvas = True

    def paint(self, p: QPainter, canvas: "Canvas") -> None:
        """Paint the molecule of a canvas."""
        frame = canvas.drawing.frame()
        paint_scene(p, canvas.width(), canvas.height(), frame.scene(canvas.scenes.base), canvas.camera,
                    canvas.theme, canvas.display, meshes=frame.layers())

    def redraw(self, canvas: "Canvas") -> None:
        """Nothing to do: the canvas repaints the molecule with its overlay."""

    def render_export(self, canvas: "Canvas", scale: float, background: Background,
                      ambient_occlusion: bool) -> QImage:
        """Picture of the current view (QPainter has no ambient occlusion)."""
        return painted_export(canvas, scale, background)


def painted_export(canvas: "Canvas", scale: float, background: Background) -> QImage:
    """Picture of a canvas's view drawn with QPainter, with its overlay."""
    theme = export_theme(canvas.theme, background)
    w, h = canvas.width(), canvas.height()
    pm = QPixmap(int(w * scale), int(h * scale))
    pm.setDevicePixelRatio(scale)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    items = canvas.drawing.scene_items(theme, export=True)
    scene = build_scene(canvas.mol, canvas.display, theme, colors=items.atom_colors,
                        display=canvas.doc.display)
    scene.arrows = items.arrows
    layers = [(item.mesh, item.alpha) for item in items.meshes]
    paint_scene(p, w, h, scene, canvas.camera, theme, canvas.display,
                background != Background.TRANSPARENT, meshes=layers)
    p.end()
    image = pm.toImage()
    canvas.drawing.paint_export_overlay(image, theme)
    return image


class GLView(QOpenGLWidget):
    """The one OpenGL widget of a window: draws the molecule of the canvas shown now.

    The stage keeps it under its canvases, at their size. Switching tabs
    sends the new canvas's arrays to the GPU once; the renderer falls back
    to QPainter if the shaders fail.
    """

    draws_in_canvas = False

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.renderer = GLRenderer()
        self.canvas: Canvas | None = None
        self._forget_uploads()

    def _forget_uploads(self) -> None:
        """Send everything again at the next paint (another canvas, or after an export)."""
        self._uploaded = False
        self._shown_arrows: Mesh | None = None
        self._shown_layers: list[int] = []
        self._shown_alphas: list[float] = []
        self._gpu_flags: np.ndarray | None = None   # highlights as last sent to the GPU

    def show_canvas(self, canvas: "Canvas | None") -> None:
        """Draw this canvas's molecule from now on."""
        if canvas is self.canvas:
            return
        self.canvas = canvas
        self._forget_uploads()
        self.update()

    def paint(self, p: QPainter, canvas: "Canvas") -> None:
        """Not used: the molecule is drawn underneath the canvas."""

    def redraw(self, canvas: "Canvas") -> None:
        """Draw again soon if ``canvas`` is the one shown."""
        if canvas is self.canvas:
            self.update()

    @override
    def initializeGL(self) -> None:
        if not self.renderer.initialize():
            print("OpenGL renderer unavailable, drawing with QPainter:", self.renderer.error)
        context = self.context()
        if context is not None:   # when the view moves to another window, Qt makes a new context
            context.aboutToBeDestroyed.connect(self._release_gl)

    @override
    def deleteLater(self) -> None:
        # Qt destroys the context in the base class destructor, after this object's slots are
        # disconnected, so the GPU objects are deleted here, while the context still exists
        self._release_gl()
        super().deleteLater()

    def _release_gl(self) -> None:
        """Delete the renderer's GPU objects, with the view's context current."""
        if not self.renderer.ok:
            return
        self.makeCurrent()
        try:
            self.renderer.release()
        finally:
            self.doneCurrent()
        self._forget_uploads()

    @override
    def paintGL(self) -> None:
        canvas = self.canvas
        if canvas is None:
            return
        canvas.view.fit_if_needed()
        frame = canvas.drawing.frame()
        p = QPainter(self)
        if self.renderer.ok:
            p.beginNativePainting()
            if frame.spheres_changed or not self._uploaded:
                self.renderer.upload_spheres(frame.spheres)
                self._gpu_flags = None
            self._upload_flags(frame.flags)
            if frame.cylinders_changed or not self._uploaded:
                self.renderer.upload_cylinders(frame.cylinders)
            if frame.mesh_changed or not self._uploaded:
                self.renderer.upload_mesh(frame.mesh, frame.mesh_colors)
            if frame.arrows_changed or frame.arrow_mesh is not self._shown_arrows or not self._uploaded:
                self.renderer.upload_mesh(frame.arrow_mesh, slot=MeshSlot.ARROWS)
                self._shown_arrows = frame.arrow_mesh
            layers = [id(item.mesh) for item in frame.meshes]
            alphas = [item.alpha for item in frame.meshes]
            if layers != self._shown_layers or not self._uploaded:
                self.renderer.upload_layers([(item.mesh, item.alpha) for item in frame.meshes])
                self._shown_layers = layers
            elif alphas != self._shown_alphas:   # the same meshes, more or less see-through
                self.renderer.set_layer_alphas(alphas)
            self._shown_alphas = alphas
            self._uploaded = True
            viewport = Viewport(self.width(), self.height(), self.devicePixelRatioF())
            self.renderer.draw(canvas.camera, viewport, canvas.theme, canvas.display.depth_fade,
                               canvas.scenes.depth(canvas.camera))
            p.endNativePainting()
        else:
            scene = frame.scene(canvas.scenes.base)
            paint_scene(p, self.width(), self.height(), scene, canvas.camera, canvas.theme, canvas.display,
                        meshes=frame.layers())
        p.end()

    def _upload_flags(self, flags: np.ndarray) -> None:
        """Send the highlights that differ from those on the GPU (all after new spheres)."""
        shown = self._gpu_flags
        if shown is None or len(shown) != len(flags):
            self.renderer.upload_flags(flags)
        else:
            rows = np.flatnonzero(flags != shown)
            if not len(rows):
                return
            self.renderer.upload_flags(flags, rows)
        self._gpu_flags = flags

    def render_export(self, canvas: "Canvas", scale: float, background: Background,
                      ambient_occlusion: bool) -> QImage:
        """Picture of a canvas's view drawn by the GPU, with its overlay."""
        if not self.renderer.ok:
            return painted_export(canvas, scale, background)
        from OpenGL import GL

        theme = export_theme(canvas.theme, background)
        transparent = background == Background.TRANSPARENT
        self.makeCurrent()
        try:
            items = canvas.drawing.scene_items(theme, export=True)
            scene = build_scene(canvas.mol, canvas.display, theme, colors=items.atom_colors,
                                display=canvas.doc.display)
            self.renderer.upload_spheres(scene.spheres)
            self.renderer.upload_cylinders(scene.cylinders)
            self.renderer.upload_mesh(scene.mesh)
            self.renderer.upload_mesh(arrow_mesh(items.arrows), slot=MeshSlot.ARROWS)
            self.renderer.upload_layers([(item.mesh, item.alpha) for item in items.meshes])
            fade = canvas.display.depth_fade and not transparent
            depth = canvas.scenes.depth(canvas.camera)
            pixels = self.renderer.render_image(canvas.camera, canvas.width(), canvas.height(), scale, theme,
                                                fade, depth, transparent, ambient_occlusion)
        finally:
            self._forget_uploads()   # the shown canvas's arrays go up again at the next paint
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.defaultFramebufferObject())
            self.doneCurrent()
        h, w = pixels.shape[:2]
        image = QImage(pixels.data, w, h, 4 * w, QImage.Format.Format_RGBA8888_Premultiplied).copy()
        canvas.drawing.paint_export_overlay(image, theme)
        self.update()
        return image


def configure_gl() -> None:
    """Ask for an OpenGL 4.1 core context with multisampling (call before creating the application)."""
    fmt = QSurfaceFormat()
    fmt.setVersion(4, 1)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    fmt.setDepthBufferSize(24)
    fmt.setSamples(4)
    QSurfaceFormat.setDefaultFormat(fmt)


_gl_available: bool | None = None


def gl_available() -> bool:
    """Whether an OpenGL 3.3+ context can be created (checked once)."""
    global _gl_available
    if _gl_available is None:
        _gl_available = False
        if not os.environ.get(SOFTWARE_ENV):
            try:
                importlib.import_module("OpenGL.GL")
                context = QOpenGLContext()
                _gl_available = bool(context.create()) and context.format().majorVersion() >= 3
            except Exception:  # missing PyOpenGL or no usable driver
                _gl_available = False
    return _gl_available


def make_surface(parent: QWidget) -> "PainterSurface | GLView":
    """Return the best surface for this machine: one OpenGL view under ``parent``, or QPainter."""
    if gl_available():
        return GLView(parent)
    return PainterSurface()
