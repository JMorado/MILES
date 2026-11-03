"""The area below the tab bar: stacked canvases with panels floating above them."""

from typing import override

from PySide6.QtGui import QResizeEvent
from PySide6.QtWidgets import QStackedWidget, QWidget

from miles.render.theme import Theme
from miles.ui.builder import Builder
from miles.ui.busy_bar import BusyBar
from miles.ui.canvas import Canvas, GLView, make_surface
from miles.ui.console import ConsolePanel
from miles.ui.floating import FloatingPanel
from miles.ui.focus_bar import FocusBar
from miles.ui.rail import Rail
from miles.ui.selection_panel import SelectionPanel
from miles.ui.sequence_strip import SequenceStrip
from miles.ui.tasks import TaskRunner
from miles.ui.timeline import Timeline

MARGIN = 14
BOTTOM = 40   # space kept free for the status line
GAP = 8


class Stage(QWidget):
    """Canvases (one per tab, stacked) and the floating panels of the current one.

    The rail sits at the top left with the sequence strip beside it, the
    focus bar at the top right, and the timeline and selection panel at the
    bottom centre. With OpenGL, one :class:`~miles.ui.canvas.GLView` under
    the canvases draws the molecule of the canvas shown (:attr:`surface`).
    """

    def __init__(self, builder: Builder, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.tasks = TaskRunner()   # slow commands running in the background
        self.canvas: Canvas | None = None
        self.surface = make_surface(self)   # what draws the molecule, shared by every canvas
        self.stack = QStackedWidget(self)
        self.rail = Rail(builder, self)
        self.strip = SequenceStrip(self)
        self.focus = FocusBar(self)
        self.busy = BusyBar(self.tasks, self)
        self.timeline = Timeline(self)
        self.console = ConsolePanel(self)
        self.extra_panels: dict[str, FloatingPanel] = {}   # added by plugins
        self.bottom_panels: list[FloatingPanel] = []   # those of them centred above the timeline
        self._theme: Theme | None = None
        self.panel = SelectionPanel(self)
        self.panel.hide()
        for widget in (self.strip, self.focus, self.busy, self.timeline, self.console, self.panel):
            widget.relayout.connect(self.layout_overlays)

    def floating_panels(self) -> dict[str, FloatingPanel]:
        """Return the panels the user can move and fold, by the name their arrangement is saved under."""
        return {"python": self.console, **self.extra_panels}

    def add_panel(self, name: str, panel: FloatingPanel, place: str = "right") -> None:
        """Take a plugin's panel into the right-hand column, or above the timeline.

        Parameters
        ----------
        name : str
            Name its arrangement is saved under.
        panel : FloatingPanel
            The panel; it is reparented to the stage.
        place : str, optional
            ``"right"``: under the panels added before it on the right;
            ``"bottom"``: centred above the timeline (and the panels added
            there before it), where fitting the view keeps clear of it.
        """
        if place not in ("right", "bottom"):
            raise ValueError(f"A panel goes on the right or at the bottom, not {place!r}")
        panel.setParent(self)
        panel.hide()
        panel.relayout.connect(self.layout_overlays)
        self.extra_panels[name] = panel
        if place == "bottom":
            self.bottom_panels.append(panel)
        if self._theme is not None:
            panel.set_icon_color(self._theme.muted)

    def remove_panel(self, name: str) -> None:
        """Take a plugin's panel out again (when its plugin is switched off)."""
        panel = self.extra_panels.pop(name, None)
        if panel is None:
            return
        if panel in self.bottom_panels:
            self.bottom_panels.remove(panel)
        panel.hide()
        panel.relayout.disconnect(self.layout_overlays)
        panel.deleteLater()
        self.layout_overlays()

    def set_canvas(self, canvas: Canvas) -> None:
        """Show a canvas and point the panels at it."""
        self.canvas = canvas
        self.stack.setCurrentWidget(canvas)
        if isinstance(self.surface, GLView):
            self.surface.show_canvas(canvas)
        self.panel.set_canvas(canvas)
        self.panel.refresh()
        self.timeline.set_document(canvas.doc)
        self.strip.set_canvas(canvas)
        self.focus.set_canvas(canvas)
        self.layout_overlays()

    def doc_changed(self, kind: str) -> None:
        """Update the panels after an edit of the current document."""
        self.panel.doc_changed(kind)
        self.timeline.doc_changed(kind)
        self.strip.refresh()
        self.focus.refresh()

    def selection_changed(self) -> None:
        """Update the panels after the selection changed."""
        self.panel.selection_changed()
        self.strip.update()

    def tool_changed(self) -> None:
        """Show or hide the selection panel for the new tool."""
        self.panel.refresh()

    def set_theme(self, theme: Theme) -> None:
        """Recolour the panels."""
        self.rail.set_theme(theme)
        self.panel.set_theme(theme)
        self.timeline.set_theme(theme)
        self.console.set_theme(theme)
        self._theme = theme
        for extra in self.extra_panels.values():
            extra.set_icon_color(theme.muted)
        self.strip.theme = theme
        self.strip.update()

    @override
    def resizeEvent(self, event: QResizeEvent) -> None:
        self.stack.setGeometry(self.rect())
        if isinstance(self.surface, GLView):   # under the canvases, which only paint their overlay
            self.surface.setGeometry(self.rect())
            self.surface.lower()
        self.layout_overlays()
        super().resizeEvent(event)

    def layout_overlays(self) -> None:
        """Place the floating panels around the edges of the canvas."""
        w, h = self.width(), self.height()
        self.rail.adjustSize()
        self.rail.move(MARGIN, MARGIN)
        self.rail.raise_()

        top = MARGIN
        if self.strip.isVisibleTo(self):
            x = MARGIN + self.rail.width() + 12
            self.strip.setGeometry(x, MARGIN, max(100, w - x - MARGIN), self.strip.height())
            self.strip.raise_()
            top = MARGIN + self.strip.height() + GAP
        right_top = top
        if self.focus.isVisibleTo(self):
            self.focus.adjustSize()
            self.focus.move(w - self.focus.width() - MARGIN, top)
            self.focus.raise_()
            right_top += self.focus.height() + GAP
        if self.busy.isVisibleTo(self):
            self.busy.adjustSize()
            self.busy.move(w - self.busy.width() - MARGIN, right_top)
            self.busy.raise_()
            right_top += self.busy.height() + GAP
        for panel in self.extra_panels.values():   # stacked on the right
            if panel in self.bottom_panels:
                continue
            if not panel.isVisibleTo(self):
                continue
            panel.fit()
            if not panel.place_in(self.size()):   # unless the user moved it
                panel.move(w - panel.width() - MARGIN, right_top)
                right_top += panel.height() + GAP
            panel.raise_()
        if self.canvas is not None:
            self.canvas.toast_top = top

        bottom = h - BOTTOM
        if self.console.isVisibleTo(self):   # along the bottom, under everything else that sits there
            self.console.fit()
            if not self.console.place_in(self.size()):
                height = self.console.height()
                self.console.setGeometry(MARGIN, bottom - height, w - 2 * MARGIN, height)
                bottom -= height + GAP
            self.console.raise_()
        if self.timeline.isVisibleTo(self):
            self.timeline.adjustSize()
            width = min(max(self.timeline.sizeHint().width(), 560), w - 2 * MARGIN)
            self.timeline.setGeometry((w - width) // 2, bottom - self.timeline.height(), width,
                                      self.timeline.height())
            self.timeline.raise_()
            bottom -= self.timeline.height() + GAP
        for panel in self.bottom_panels:   # centred above the timeline
            if not panel.isVisibleTo(self):
                continue
            panel.fit()
            if not panel.place_in(self.size()):
                panel.move((w - panel.width()) // 2, bottom - panel.height())
                bottom -= panel.height() + GAP
            panel.raise_()
        if self.canvas is not None:   # fitting the view keeps the structure clear of the panels at the edges
            self.canvas.view_insets = (top, max(0, h - bottom - BOTTOM))
        if self.panel.isVisibleTo(self):
            self.panel.adjustSize()
            self.panel.move((w - self.panel.width()) // 2, bottom - self.panel.height())
            self.panel.raise_()
