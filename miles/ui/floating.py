"""Panels floating over the canvas that can be moved by dragging and collapsed to their top row.

:class:`FloatingPanel` is the panel itself; :class:`DocumentPanel` adds
what every panel about the tab in front needs: it follows that tab, opens
and closes on request, and keeps data for each tab.
"""

from typing import TYPE_CHECKING, Any, ClassVar, override

from PySide6.QtCore import QPoint, QSize, Qt, Signal
from PySide6.QtGui import QColor, QHideEvent, QMouseEvent, QShowEvent
from PySide6.QtWidgets import QFrame, QToolButton, QVBoxLayout, QWidget

from miles.core.document import Change
from miles.ui import icons

if TYPE_CHECKING:
    from miles.api import DocumentView
    from miles.plugins.context import PluginContext
    from miles.render.theme import Theme

EDGE = 4   # pixels a moved panel keeps from the edges of the window


class FloatingPanel(QFrame):
    """A floating panel with a body that folds away and a place the user can choose.

    Subclasses put their top row in their own layout, with
    :attr:`chevron` at its left, and everything else in :attr:`body`
    (through :attr:`body_layout`), then add :attr:`body` to their layout.
    The chevron collapses the panel to its top row. Dragging the panel by
    its background or labels moves it anywhere inside its parent;
    :attr:`placed` then holds its position, which the stage keeps instead
    of the usual place. A double click on the background puts it back.
    """

    relayout = Signal()
    arranged = Signal()   # the user moved, folded or reset the panel (to remember it)
    opened = Signal()     # the panel appeared
    closed = Signal()     # the panel was hidden

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Panel")
        self.placed: QPoint | None = None
        self.collapsed = False
        self._shown = False   # what opened and closed last said
        self._grab: QPoint | None = None
        self._moved = False
        self._icon_color = QColor(128, 128, 128)
        self.chevron = QToolButton()
        self.chevron.setObjectName("Order")
        self.chevron.setFixedSize(22, 22)
        self.chevron.setIconSize(QSize(16, 16))
        self.chevron.setCursor(Qt.CursorShape.PointingHandCursor)
        self.chevron.clicked.connect(self.toggle_collapsed)
        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 0, 0, 0)
        self._update_chevron()

    def set_icon_color(self, color: QColor) -> None:
        """Colour of the chevron, from the theme."""
        self._icon_color = QColor(color)
        self._update_chevron()

    def toggle_collapsed(self) -> None:
        """Fold the panel to its top row, or open it again."""
        self.set_collapsed(not self.collapsed)

    def set_collapsed(self, collapsed: bool) -> None:
        """Show only the top row (True) or the whole panel (False)."""
        changed = collapsed != self.collapsed
        self.collapsed = collapsed
        self.body.setVisible(not collapsed)
        self._update_chevron()
        self.fit()
        self.relayout.emit()
        if changed:
            self.arranged.emit()

    def fit(self) -> None:
        """Take the size the contents need."""
        self.adjustSize()

    def open_panel(self) -> bool:
        """Show the panel; returns whether it is shown (a subclass may have nothing to show)."""
        self.setVisible(True)
        self.fit()
        self.relayout.emit()
        return True

    def close_panel(self) -> None:
        """Hide the panel until it is opened again."""
        self.hide()
        self.relayout.emit()

    def toggle(self) -> bool:
        """Open the panel, or close it when it is open; returns whether it is open now."""
        if self.isVisible():
            self.close_panel()
            return False
        return self.open_panel()

    @override
    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        if not event.spontaneous() and not self._shown:
            self._shown = True
            self.opened.emit()

    @override
    def hideEvent(self, event: QHideEvent) -> None:
        super().hideEvent(event)
        if not event.spontaneous() and self._shown and not self.isVisible():
            self._shown = False
            self.closed.emit()

    def place_in(self, area: QSize) -> bool:
        """Move to the place chosen by dragging, kept inside ``area``; False when there is none."""
        if self.placed is None:
            return False
        self.move(self._clamped(self.placed, area))
        return True

    def _clamped(self, point: QPoint, area: QSize) -> QPoint:
        x = max(EDGE, min(point.x(), area.width() - self.width() - EDGE))
        y = max(EDGE, min(point.y(), area.height() - self.height() - EDGE))
        return QPoint(x, y)

    def _update_chevron(self) -> None:
        self.chevron.setIcon(icons.icon("expand" if self.collapsed else "collapse", self._icon_color))
        self.chevron.setToolTip("Show the whole panel" if self.collapsed else "Fold the panel to this row")

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if e.button() == Qt.MouseButton.LeftButton:
            self._grab = e.position().toPoint()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            e.accept()
        else:
            super().mousePressEvent(e)

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        parent = self.parentWidget()
        if self._grab is None or parent is None:
            super().mouseMoveEvent(e)
            return
        self.placed = self._clamped(self.mapToParent(e.position().toPoint()) - self._grab, parent.size())
        self.move(self.placed)
        self.raise_()
        self._moved = True

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        if self._grab is not None:
            self._grab = None
            self.unsetCursor()
            if self._moved:
                self._moved = False
                self.arranged.emit()
        super().mouseReleaseEvent(e)

    @override
    def mouseDoubleClickEvent(self, e: QMouseEvent) -> None:
        if e.button() == Qt.MouseButton.LeftButton and self.placed is not None:
            self.placed = None   # back to the usual place
            self.relayout.emit()
            self.arranged.emit()
        super().mouseDoubleClickEvent(e)


class DocumentPanel(FloatingPanel):
    """A floating panel about the tab in front, such as the Ramachandran plot or the vibrations.

    Register it with :meth:`~miles.plugins.PluginContext.add_document_panel`,
    which makes it follow the tab in front, its edits and (when
    :attr:`follows_selection` is set) its selection, and adds a command
    that opens and closes it. A subclass says whether a tab has anything to
    show (:meth:`has_content`) and fills the panel (:meth:`show_contents`);
    the panel shows itself only while it is open and the tab has content.

    Attributes
    ----------
    document : DocumentView or None
        The tab the panel shows.
    wanted : bool
        Whether the user opened the panel (it may still be hidden when the
        tab in front has nothing to show).
    context : PluginContext or None
        The plugin's context, set when the panel is registered.
    refresh_on : frozenset of str
        Kinds of change of the tab (see :class:`~miles.plugin.Change`)
        after which the panel is filled again.
    follows_selection : bool
        Whether the panel is filled again when the selection changes.
    unavailable : str
        What the command says when the tab in front has nothing to show.
    """

    refresh_on: ClassVar[frozenset[str]] = frozenset({Change.STRUCTURE, Change.BONDS})
    follows_selection: ClassVar[bool] = False
    unavailable: ClassVar[str] = "Nothing to show for this structure"

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.document: DocumentView | None = None
        self.wanted = False
        self.context: PluginContext | None = None
        self.hide()

    # -- what a subclass gives ------------------------------------------------------------------

    def has_content(self, document: "DocumentView") -> bool:
        """Whether a tab has anything to show (default: any tab)."""
        return True

    def show_contents(self, document: "DocumentView") -> None:
        """Fill the panel for a tab; called while the panel is open and the tab has content."""

    def closing(self) -> None:
        """React to the user closing the panel, e.g. stop an animation it started (default: nothing)."""

    def set_theme(self, theme: "Theme") -> None:
        """Take the colours of the interface (default: the chevron's)."""
        self.set_icon_color(theme.muted)

    # -- per-tab data ---------------------------------------------------------------------------

    def store(self, document: "DocumentView") -> dict[str, Any]:
        """Return a dictionary for working data about a tab, dropped when the tab closes."""
        if self.context is None:
            raise RuntimeError("The panel is not registered with add_document_panel")
        return self.context.document_cache(document)

    def saved(self, document: "DocumentView") -> dict[str, Any]:
        """Return a dictionary for data about a tab that session files keep (JSON values only)."""
        if self.context is None:
            raise RuntimeError("The panel is not registered with add_document_panel")
        return self.context.document_state(document)

    # -- following the tab in front ----------------------------------------------------------

    def set_document(self, document: "DocumentView | None") -> None:
        """Show another tab (MILES calls this when the tab in front changes)."""
        self.document = document
        self.refresh()

    def document_changed(self, kind: str) -> None:
        """React to a change of the tab shown (default: fill again for the kinds in :attr:`refresh_on`)."""
        if kind in self.refresh_on:
            self.refresh()

    def selection_changed(self) -> None:
        """React to a new selection in the tab shown (default: fill again if :attr:`follows_selection`)."""
        if self.follows_selection:
            self.refresh()

    def refresh(self) -> None:
        """Fill the panel for the tab shown and show it; hide it when it is closed or has nothing to show."""
        doc = self.document
        visible = self.wanted and doc is not None and self.has_content(doc)
        if visible and doc is not None:
            self.show_contents(doc)
        if visible != self.isVisible():
            self.setVisible(visible)
        if visible:
            self.fit()
        self.relayout.emit()

    @override
    def open_panel(self) -> bool:
        """Open the panel; returns False, leaving it closed, when the tab in front has nothing to show."""
        if self.document is None or not self.has_content(self.document):
            return False
        self.wanted = True
        self.refresh()
        return True

    @override
    def close_panel(self) -> None:
        """Close the panel until it is opened again."""
        if self.wanted:
            self.closing()
        self.wanted = False
        self.refresh()

    @override
    def toggle(self) -> bool:
        """Open the panel, or close it when it is open; returns whether it is open now."""
        if self.wanted and self.isVisible():
            self.close_panel()
            return False
        return self.open_panel()
