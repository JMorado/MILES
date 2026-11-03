"""The open documents of a window, in tab order, and which one is shown."""

from PySide6.QtCore import QObject, Signal
from PySide6.QtWidgets import QStackedWidget, QTabBar

from miles.ui.canvas import Canvas


class DocumentTabs(QObject):
    """The tabs of a window: the tab bar, with each tab's canvas kept in its data.

    The tab bar is the one record of which documents are open and in what
    order, so tabs can be dragged into another order; the stack of canvases
    only shows the current one.

    Parameters
    ----------
    bar : QTabBar
        The window's tab bar; it is hidden while only one document is open.
    stack : QStackedWidget
        Where the canvases live.

    Signals
    -------
    currentChanged(canvas)
        Another tab came to the front, once per change.
    """

    currentChanged = Signal(object)

    def __init__(self, bar: QTabBar, stack: QStackedWidget, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.bar = bar
        self.stack = stack
        self._shown: Canvas | None = None
        bar.setMovable(True)
        bar.currentChanged.connect(lambda _index: self._follow_bar())

    @property
    def canvases(self) -> list[Canvas]:
        """Canvases of the open tabs, in tab order."""
        return [self.bar.tabData(i) for i in range(self.bar.count())]

    @property
    def current(self) -> Canvas | None:
        """Canvas of the tab in front."""
        index = self.bar.currentIndex()
        return self.bar.tabData(index) if index >= 0 else None

    def index_of(self, canvas: Canvas) -> int:
        """Position of a canvas's tab, or -1 if it has none."""
        return next((i for i in range(self.bar.count()) if self.bar.tabData(i) is canvas), -1)

    def add(self, canvas: Canvas, title: str) -> None:
        """Open a tab for a canvas and bring it to the front."""
        self.stack.addWidget(canvas)
        index = self.bar.addTab(title)
        self.bar.setTabData(index, canvas)
        self.bar.setCurrentIndex(index)
        self._follow_bar()
        self.bar.setVisible(self.bar.count() > 1)

    def remove(self, canvas: Canvas) -> None:
        """Close a canvas's tab; the canvas is deleted."""
        index = self.index_of(canvas)
        if index < 0:
            return
        self.stack.removeWidget(canvas)
        self.bar.removeTab(index)
        self._follow_bar()
        self.bar.setVisible(self.bar.count() > 1)
        canvas.deleteLater()

    def set_current(self, canvas: Canvas) -> None:
        """Bring a canvas's tab to the front."""
        index = self.index_of(canvas)
        if index >= 0:
            self.bar.setCurrentIndex(index)

    def set_title(self, canvas: Canvas, title: str) -> None:
        """Change the text of a canvas's tab."""
        index = self.index_of(canvas)
        if index >= 0:
            self.bar.setTabText(index, title)

    def _follow_bar(self) -> None:
        canvas = self.current
        if canvas is None or canvas is self._shown:   # no tab yet, or only the order changed
            if canvas is None:
                self._shown = None
            return
        self._shown = canvas
        self.stack.setCurrentWidget(canvas)
        self.currentChanged.emit(canvas)
