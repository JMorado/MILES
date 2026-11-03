"""Command palette: type to find and run any menu command."""

from dataclasses import dataclass
from typing import override

from PySide6.QtCore import QEvent, QModelIndex, QObject, QPersistentModelIndex, QPoint, QRect, QSize, Qt
from PySide6.QtGui import QAction, QKeyEvent, QPainter
from PySide6.QtWidgets import (
    QFrame,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QMenuBar,
    QStyle,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QVBoxLayout,
    QWidget,
)

from miles.render.theme import Theme

MAX_ROWS = 12
SHORTCUT_ROLE = Qt.ItemDataRole.UserRole + 1


@dataclass(frozen=True)
class Entry:
    """A menu command: where it lives in the menus and its action."""

    path: str
    action: QAction

    @property
    def shortcut(self) -> str:
        """Its keyboard shortcut in the platform's notation, or ""."""
        return self.action.shortcut().toString(self.action.shortcut().SequenceFormat.NativeText)


def menu_entries(bar: QMenuBar) -> list[Entry]:
    """Every command of a menu bar, in menu order, with its menu path."""
    entries: list[Entry] = []

    def walk(menu: QMenu, prefix: str) -> None:
        for action in menu.actions():
            if action.isSeparator() or not action.text():
                continue
            submenu = action.menu()
            name = action.text().replace("&", "").removesuffix("...")
            if isinstance(submenu, QMenu):
                walk(submenu, f"{prefix}{name} > ")
            else:
                entries.append(Entry(prefix + name, action))

    for top in bar.actions():
        menu = top.menu()
        if isinstance(menu, QMenu):
            walk(menu, f"{top.text().replace('&', '')} > ")
    return entries


NAME_BONUS = 5.0
"""Added when the query matches the command's own name, not just its menu path."""


def match_score(query: str, path: str) -> float | None:
    """How well ``query`` matches a command (higher is better), or None if it does not.

    Every query character must appear in order. Matches at word starts and
    runs of consecutive characters score higher, and so do shorter texts.
    A match within the command's name beats a match spread over its menus.
    """
    name = path.rsplit(" > ", 1)[-1]
    in_name = _subsequence_score(query, name)
    if in_name is not None:
        return in_name + NAME_BONUS
    return _subsequence_score(query, path)


def _subsequence_score(query: str, text: str) -> float | None:
    query, lowered = query.lower().replace(" ", ""), text.lower()
    if not query:
        return 0.0
    score, position, previous = 0.0, 0, -2
    for ch in query:
        found = lowered.find(ch, position)
        if found < 0:
            return None
        if found == previous + 1:
            score += 3.0
        if found == 0 or not lowered[found - 1].isalnum():
            score += 2.0
        score += 1.0
        previous, position = found, found + 1
    return score - 0.01 * len(text)


class _Delegate(QStyledItemDelegate):
    """Paints the command on the left and its shortcut, muted, on the right."""

    def __init__(self, palette: "CommandPalette") -> None:
        super().__init__(palette)
        self.palette_widget = palette

    @override
    def paint(self, painter: QPainter, option: QStyleOptionViewItem,
              index: QModelIndex | QPersistentModelIndex) -> None:
        theme = self.palette_widget.theme
        if theme is None:
            return
        rect: QRect = option.rect
        selected = bool(option.state & QStyle.StateFlag.State_Selected)
        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        if selected:
            color = theme.accent.lighter(170) if not theme.dark else theme.accent.darker(170)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawRoundedRect(rect.adjusted(2, 1, -2, -1), 7, 7)
        inner = rect.adjusted(12, 0, -12, 0)
        painter.setPen(theme.muted)
        painter.drawText(inner, Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignRight,
                         str(index.data(SHORTCUT_ROLE) or ""))
        painter.setPen(theme.text)
        painter.drawText(inner.adjusted(0, 0, -110, 0), Qt.AlignmentFlag.AlignVCenter, str(index.data()))
        painter.restore()

    @override
    def sizeHint(self, option: QStyleOptionViewItem, index: QModelIndex | QPersistentModelIndex) -> QSize:
        return QSize(560, 30)


class CommandPalette(QWidget):
    """Popover with a search field and the matching menu commands."""

    def __init__(self, bar: QMenuBar, parent: QWidget | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.bar = bar
        self.theme: Theme | None = None
        self._entries: list[Entry] = []
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        frame = QFrame()
        frame.setObjectName("Popover")
        outer.addWidget(frame)
        lay = QVBoxLayout(frame)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(6)
        self.search = QLineEdit()
        self.search.setObjectName("Search")
        self.search.setPlaceholderText("Type a command, e.g. surface, hydrogens, export")
        self.search.textChanged.connect(self._filter)
        self.search.installEventFilter(self)
        lay.addWidget(self.search)
        self.results = QListWidget()
        self.results.setObjectName("Palette")
        self.results.setItemDelegate(_Delegate(self))
        self.results.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.results.itemActivated.connect(self._run)
        self.results.itemClicked.connect(self._run)
        lay.addWidget(self.results)
        self.setFixedWidth(600)

    def set_theme(self, theme: Theme) -> None:
        """Recolour the rows."""
        self.theme = theme

    def open_at(self, pos: QPoint) -> None:
        """Show the palette with every command, the search field focused."""
        self._entries = [e for e in menu_entries(self.bar) if e.action.isEnabled()]
        self.search.clear()
        self._filter("")
        self.move(pos)
        self.show()
        self.search.setFocus()

    def _filter(self, query: str) -> None:
        scored = [(score, n, e) for n, e in enumerate(self._entries)
                  if (score := match_score(query, e.path)) is not None]
        scored.sort(key=lambda t: (-t[0], t[1]) if query else (t[1], 0))
        self.results.clear()
        for _, _, entry in scored:
            item = QListWidgetItem(entry.path)
            item.setData(SHORTCUT_ROLE, entry.shortcut)
            item.setData(Qt.ItemDataRole.UserRole, entry)
            self.results.addItem(item)
        if self.results.count():
            self.results.setCurrentRow(0)
        rows = min(MAX_ROWS, max(1, self.results.count()))
        self.results.setFixedHeight(rows * 30 + 4)
        self.adjustSize()

    def _run(self, item: QListWidgetItem) -> None:
        entry = item.data(Qt.ItemDataRole.UserRole)
        self.close()
        if isinstance(entry, Entry):
            entry.action.trigger()

    @override
    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        if watched is self.search and isinstance(event, QKeyEvent) and event.type() == QEvent.Type.KeyPress:
            key = event.key()
            if key in (Qt.Key.Key_Down, Qt.Key.Key_Up):
                step = 1 if key == Qt.Key.Key_Down else -1
                row = max(0, min(self.results.count() - 1, self.results.currentRow() + step))
                self.results.setCurrentRow(row)
                return True
            if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
                item = self.results.currentItem()
                if item is not None:
                    self._run(item)
                return True
            if key == Qt.Key.Key_Escape:
                self.close()
                return True
        return False
