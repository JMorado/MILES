"""The canvas: one document's view, put together from parts that each do one job."""

import copy
from collections.abc import Iterable
from typing import Any, cast, override

import numpy as np
from PySide6.QtCore import QEvent, QPointF, Qt, QTimer, Signal
from PySide6.QtGui import (
    QEnterEvent,
    QHelpEvent,
    QImage,
    QKeyEvent,
    QMouseEvent,
    QNativeGestureEvent,
    QPainter,
    QPaintEvent,
    QResizeEvent,
    QWheelEvent,
)
from PySide6.QtWidgets import QWidget

from miles.core import stereo
from miles.core.molecule import Molecule
from miles.core.symmetry import find_point_group
from miles.render.camera import Camera
from miles.render.picking import Hit, Picker
from miles.render.scene import Style, visible_atoms
from miles.render.theme import LIGHT, Background, Theme
from miles.ui.builder import Builder, MouseMode, Tool
from miles.ui.canvas.commands import Commands
from miles.ui.canvas.drawing import Drawing
from miles.ui.canvas.extensions import Extensions
from miles.ui.canvas.input import InputController
from miles.ui.canvas.overlay import paint_overlay
from miles.ui.canvas.scene_cache import SceneCache
from miles.ui.canvas.surfaces import PainterSurface, Surface
from miles.ui.canvas.tools import TOOLS, CanvasTool, Gesture
from miles.ui.canvas.view import ViewController
from miles.ui.document import Change, Document
from miles.ui.tasks import TaskRunner
from miles.ui.text import SEPARATOR, subscript_digits

INFO_DELAY_MS = 120

POINT_GROUP_LIMIT = 500     # larger structures leave the point group out of the status line

MODE_CURSORS = {
    MouseMode.ROTATE: Qt.CursorShape.ArrowCursor,
    MouseMode.TRANSLATE: Qt.CursorShape.SizeAllCursor,
    MouseMode.SCALE: Qt.CursorShape.SizeVerCursor,
}


class Canvas(QWidget):
    """View of one document: what is shown of it, the selection, input and the status line.

    The canvas owns parts that each do one job:

    * ``view`` (:class:`~miles.ui.canvas.view.ViewController`): the
      camera and smooth moves (fit, centre, zoom);
    * ``scenes`` (:class:`~miles.ui.canvas.scene_cache.SceneCache`):
      the arrays to draw, kept until their inputs change;
    * ``drawing`` (:class:`~miles.ui.canvas.drawing.Drawing`): the
      frame and what plugins draw;
    * ``input`` (:class:`~miles.ui.canvas.input.InputController`):
      mouse, trackpad and keyboard;
    * ``commands`` (:class:`~miles.ui.canvas.commands.Commands`): what
      the menus do to the document and the view;
    * ``surface`` (:mod:`~miles.ui.canvas.surfaces`): what draws the
      molecule, QPainter inside the canvas or the window's OpenGL view
      underneath it. The canvas always paints its overlay itself.

    Parameters
    ----------
    doc : miles.ui.document.Document
        The document shown.
    builder : Builder
        The window's tool, fragment and mouse mode, shared by all canvases.
    surface : Surface, optional
        What draws the molecule (default: QPainter in the canvas).
    extensions : Extensions, optional
        What plugins draw and add to the right-click menu.
    tasks : TaskRunner, optional
        Where slow commands run on large structures, shown in the busy bar
        (default: they run at once).
    parent : QWidget, optional
        Parent widget.

    Signals
    -------
    selectionChanged()
        The selected atoms changed.
    focused()
        The canvas was clicked.
    nudge(steps)
        Option-scroll over a selection of 2-4 atoms, in dial ticks.
    playRequested()
        Space was pressed: start or stop playing the frames.
    documentChanged(kind)
        The document changed, after the canvas followed it; ``kind`` is a
        :class:`~miles.ui.document.Change`. Listening here rather than to
        the document lets the connection go with the canvas when its tab
        closes.
    documentRenamed()
        The document's name, file or modified state changed.
    """

    selectionChanged = Signal()
    focused = Signal()
    nudge = Signal(float)
    playRequested = Signal()
    documentChanged = Signal(str)
    documentRenamed = Signal()

    def __init__(self, doc: Document, builder: Builder, surface: Surface | None = None,
                 extensions: Extensions | None = None, tasks: TaskRunner | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.doc = doc
        self.builder = builder
        self.surface: Surface = surface if surface is not None else PainterSurface()
        self.extensions = extensions if extensions is not None else Extensions()
        self.tasks = tasks
        self.theme: Theme = LIGHT
        self.display = Style()
        self.hover: Hit | None = None
        self.gesture: Gesture | None = None
        self.info_text = ""
        self.toast_text = ""
        self.toast_top = 16        # moved down by the stage when the sequence strip is shown
        self.view_insets = (0, 0)   # pixels at the top and bottom that panels cover, set by the stage
        self.placed_last = False   # the last build click placed a free fragment
        self.plugin_state: dict[str, dict[str, Any]] = copy.deepcopy(doc.info.view.get("plugins") or {})
        self._selection_seen: list[int] = list(doc.selected)   # to tell when edits or undo move it
        self._stereo_key: tuple[int, int] | None = None
        self._stereo: tuple[dict[int, stereo.Chirality], dict[tuple[int, int], stereo.BondStereo]] = ({}, {})

        self.view = ViewController(self)
        self.scenes = SceneCache()
        self.drawing = Drawing(self)
        self.input = InputController(self)
        self.commands = Commands(self)

        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumSize(200, 150)
        if self.surface.draws_in_canvas:
            self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent)
        else:   # the molecule is drawn underneath: only the overlay is painted here
            self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        self._toast_timer = QTimer(self)
        self._toast_timer.setSingleShot(True)
        self._toast_timer.timeout.connect(self._clear_toast)
        self._info_timer = QTimer(self)
        self._info_timer.setSingleShot(True)
        self._info_timer.setInterval(INFO_DELAY_MS)
        self._info_timer.timeout.connect(self.update_status)
        builder.toolChanged.connect(self._tool_changed)
        builder.fragmentChanged.connect(self.redraw)
        builder.mouseModeChanged.connect(self._mouse_mode_changed)
        doc.changed.connect(self._doc_changed)
        doc.changed.connect(self.documentChanged)
        doc.renamed.connect(self.documentRenamed)
        self.update_status()

    # -- what is shown ------------------------------------------------------------------------

    @property
    def mol(self) -> Molecule:
        """The document's molecule."""
        return self.doc.mol

    @property
    def camera(self) -> Camera:
        """The view shown now (kept by :attr:`view`)."""
        return self.view.camera

    @camera.setter
    def camera(self, camera: Camera) -> None:
        self.view.stop()
        self.view.camera = camera

    @property
    def tool(self) -> CanvasTool:
        """Behaviour of the active tool."""
        return TOOLS[self.builder.tool]

    def redraw(self) -> None:
        """Draw the molecule and the overlay again soon."""
        self.update()
        self.surface.redraw(self)

    def set_theme(self, theme: Theme) -> None:
        """Switch colours."""
        self.theme = theme
        self.scenes.invalidate()
        self.redraw()

    def set_style(self, **changes: Any) -> None:
        """Change display options, e.g. ``set_style(show_hydrogens=False)``."""
        for name, value in changes.items():
            setattr(self.display, name, value)
        if changes.keys() & {"mode", "show_hydrogens"}:
            self.update_status()   # what plugins report about the atoms drawn may have changed
        self.redraw()

    def fit(self, animate: bool = False) -> None:
        """Zoom to show the whole structure (see :meth:`ViewController.fit`)."""
        self.view.fit(animate)

    def zoom_to(self, atoms: Iterable[int]) -> None:
        """Zoom smoothly so that some atoms fill the view (see :meth:`ViewController.zoom_to`)."""
        self.view.zoom_to(atoms)

    def render_export(self, scale: float, background: Background, ambient_occlusion: bool) -> QImage:
        """Picture of the current view for publication, drawn by the canvas's surface."""
        return self.surface.render_export(self, scale, background, ambient_occlusion)

    # -- messages and the status line ---------------------------------------------------------

    def toast(self, text: str, ms: int = 1600) -> None:
        """Show a short message at the top of the canvas."""
        self.toast_text = text
        self._toast_timer.start(ms)
        self.redraw()

    def _clear_toast(self) -> None:
        self.toast_text = ""
        self.redraw()

    def update_status(self) -> None:
        """Work out the status line again: formula, mass, point group and what plugins report."""
        mol = self.mol
        parts = []
        if len(mol):
            parts = [subscript_digits(mol.formula()), f"{mol.masses().sum():.2f} g/mol"]
            if len(mol) <= POINT_GROUP_LIMIT:
                parts.append(find_point_group(mol.symbols, mol.coords, 0.1).name)
            for layer in self.extensions.status:
                text = layer.function(self, None)
                if text:
                    parts.append(str(text))
        self.info_text = SEPARATOR.join(parts)
        self.redraw()

    # -- the selection ------------------------------------------------------------------------

    @property
    def selection(self) -> list[int]:
        """Selected atoms (order matters for measurements); do not change the list in place.

        The document keeps the selection by atom, so it stays on the same
        atoms when others are added or deleted and through undo (see
        :attr:`miles.core.document.Document.selected`). Assigning selects
        without announcing it; :meth:`set_selection` announces it.
        """
        return self.doc.selected

    @selection.setter
    def selection(self, atoms: list[int]) -> None:
        self.doc.select(atoms)
        self._selection_seen = list(self.doc.selected)

    def set_selection(self, atoms: list[int]) -> None:
        """Select atoms (order matters for measurements)."""
        self.selection = list(atoms)
        self.selectionChanged.emit()
        self.redraw()

    def clear_selection(self) -> None:
        """Select nothing."""
        if self.selection:
            self.set_selection([])

    def selected_or_all(self) -> list[int]:
        """Return the selection, or every atom when nothing is selected."""
        return list(self.selection) or list(range(len(self.mol)))

    def selected_with_hydrogens(self) -> list[int]:
        """Return the selection with the hydrogens bonded to it, or every atom when nothing is selected.

        Selecting the heavy atoms of a ligand then describes the whole
        ligand: formulas, identifiers and drawings need its hydrogens.
        """
        if not self.selection:
            return list(range(len(self.mol)))
        chosen = set(self.selection)
        for i in self.selection:
            chosen.update(j for j in self.mol.neighbors(i) if self.mol.atoms[j].symbol == "H")
        return sorted(chosen)

    # -- what can be picked -------------------------------------------------------------------

    def dragging(self) -> bool:
        """Whether atoms are being dragged."""
        return self.gesture is not None and self.gesture.atoms is not None and self.gesture.moved

    def shown_atoms(self) -> np.ndarray:
        """Boolean mask of the atoms shown: drawn, or standing for their residue in a cartoon (CA or P)."""
        return self.picker().pickable

    def drawn_atoms(self) -> np.ndarray:
        """Boolean mask of the atoms drawn as spheres (not hidden, not part of a cartoon)."""
        return visible_atoms(self.mol, self.display, display=self.doc.display)

    def picker(self) -> Picker:
        """Hit tester for the structure as it is drawn now."""
        return self.scenes.picker(self.mol, self.doc.layout_revision, self.display, self.doc.display)

    def pick(self, pos: QPointF) -> Hit | None:
        """Atom or bond under a point of the canvas."""
        return self.picker().pick(self.camera, self.width(), self.height(), pos.x(), pos.y())

    def stereo_labels(self) -> tuple[dict[int, stereo.Chirality], dict[tuple[int, int], stereo.BondStereo]]:
        """R/S labels of stereocentres and E/Z labels of double bonds, recomputed after edits."""
        key = (id(self.mol), self.doc.revision)
        if key != self._stereo_key:
            self._stereo_key = key
            self._stereo = (stereo.stereocenters(self.mol), stereo.double_bond_stereo(self.mol))
        return self._stereo

    # -- following the document and the builder -----------------------------------------------

    def _doc_changed(self, kind: str) -> None:
        selected = self.doc.selected   # follows its atoms through edits; undo restores it
        if selected != self._selection_seen:
            self._selection_seen = list(selected)
            self.selectionChanged.emit()
        if kind == Change.RESULTS and self.doc.stale_results():
            stale = self.doc.stale_results()
            what = ", ".join(stale[:-1]) + " and " + stale[-1] if len(stale) > 1 else stale[0]
            self.toast(f"The {what} from the file belong to the structure before this edit; "
                       "undo to use them again", 6000)
        if kind == Change.PROPERTIES:
            self.drawing.forget_layers()   # plugins may draw what depends on the charges
        self.hover = None
        self.view.fit_if_needed()
        self._info_timer.start()
        self.redraw()

    def _tool_changed(self) -> None:
        if self.builder.tool != Tool.SELECT:
            self.selection = []
        self.selectionChanged.emit()
        self.redraw()

    def _mouse_mode_changed(self) -> None:
        mode = self.builder.mouse_mode
        self.setCursor(MODE_CURSORS[mode])
        self.toast(f"{mode.capitalize()} mode")

    # -- Qt events: painting here, input to :attr:`input` -------------------------------------

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        self.view.fit_if_needed()
        p = QPainter(self)
        if self.surface.draws_in_canvas:
            self.surface.paint(p, self)
        paint_overlay(p, self)
        p.end()

    @override
    def resizeEvent(self, e: QResizeEvent) -> None:
        self.view.fit_if_needed()
        super().resizeEvent(e)

    @override
    def enterEvent(self, e: QEnterEvent) -> None:
        self.input.enter()
        super().enterEvent(e)

    @override
    def leaveEvent(self, e: QEvent) -> None:
        self.input.leave()

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        self.input.press(e)

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        self.input.move(e)

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        self.input.release(e)

    @override
    def mouseDoubleClickEvent(self, e: QMouseEvent) -> None:
        if self.input.double_click(e):
            self.input.press(e)

    @override
    def wheelEvent(self, e: QWheelEvent) -> None:
        self.input.wheel(e)

    @override
    def keyPressEvent(self, e: QKeyEvent) -> None:
        if not self.input.key(e):
            super().keyPressEvent(e)

    @override
    def event(self, e: QEvent) -> bool:
        if e.type() == QEvent.Type.ToolTip and isinstance(e, QHelpEvent):
            self.input.tooltip(e)
            return True
        if e.type() == QEvent.Type.NativeGesture and self.input.native_gesture(cast(QNativeGestureEvent, e)):
            return True
        return super().event(e)
