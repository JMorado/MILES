"""Python console inside the window: run code on the open structures, as in PyMOL."""

import code
import contextlib
import os
import re
import rlcompleter
import sys
import tokenize
import traceback
from collections.abc import Iterator
from types import TracebackType
from typing import Any, override

from PySide6.QtCore import QEvent, QPoint, QSettings, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFontDatabase,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QTextCharFormat,
    QTextCursor,
)
from PySide6.QtWidgets import QApplication, QHBoxLayout, QLineEdit, QPlainTextEdit, QVBoxLayout, QWidget

from miles.api import App
from miles.render.theme import LIGHT, Theme
from miles.scripting import BANNER, namespace, refresh, settle, snapshot
from miles.ui.controls import button, label
from miles.ui.floating import FloatingPanel

HISTORY_LIMIT = 500        # commands kept between sessions
OUTPUT_LINES = 5000        # lines kept in the output area
DEFAULT_HEIGHT = 230       # pixels
MIN_HEIGHT = 120
PROMPT, MORE = ">>>", "..."
CLEAR_MODIFIER = (Qt.KeyboardModifier.MetaModifier if sys.platform == "darwin"
                  else Qt.KeyboardModifier.ControlModifier)   # the Control key, as in terminals


class _Interpreter(code.InteractiveConsole):
    """Interactive interpreter that writes errors to the console, never to the error dialog."""

    def __init__(self, names: dict[str, Any], console: "ConsolePanel") -> None:
        super().__init__(names, filename="<console>")
        self.console = console

    @override
    def write(self, data: str) -> None:
        self.console.write(data, "error")

    @override
    def showtraceback(self) -> None:
        kind, value, tb = sys.exc_info()
        if kind is None or value is None:
            return
        sys.last_type, sys.last_value, sys.last_traceback = kind, value, tb
        self.write("".join(traceback.format_exception(kind, value, tb.tb_next if tb else None)))

    @override
    def showsyntaxerror(self, filename: str | None = None, **kwargs: Any) -> None:
        kind, value, _ = sys.exc_info()
        if kind is not None and value is not None:
            self.write("".join(traceback.format_exception_only(kind, value)))


class _Stream:
    """Stand-in for ``sys.stdout`` or ``sys.stderr`` that writes into the console."""

    encoding = "utf-8"

    def __init__(self, console: "ConsolePanel", kind: str) -> None:
        self.console = console
        self.kind = kind

    def write(self, text: str) -> int:
        self.console.write(text, self.kind)
        return len(text)

    def flush(self) -> None:
        """Nothing is buffered."""

    def isatty(self) -> bool:
        """Not a terminal, so tools such as ``help`` print plainly."""
        return False


class ConsoleInput(QLineEdit):
    """One line of input with history, completion and multi-line paste."""

    submitted = Signal(str)
    completion = Signal()
    cleared = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.history: list[str] = []
        self._position = 0      # index into history while browsing; len(history) means the new line
        self._draft = ""

    def set_history(self, lines: list[str]) -> None:
        """Commands from earlier sessions, oldest first."""
        self.history = lines[-HISTORY_LIMIT:]
        self._position = len(self.history)

    def remember(self, line: str) -> None:
        """Add a command to the history (not blank lines or repeats of the last one)."""
        if line.strip() and (not self.history or self.history[-1] != line):
            self.history.append(line)
            del self.history[:-HISTORY_LIMIT]
        self._position = len(self.history)

    def _browse(self, step: int) -> None:
        if not self.history:
            return
        if self._position == len(self.history):
            self._draft = self.text()
        self._position = max(0, min(len(self.history), self._position + step))
        self.setText(self._draft if self._position == len(self.history) else self.history[self._position])

    def _is_clear(self, e: QKeyEvent) -> bool:
        return e.key() == Qt.Key.Key_L and bool(e.modifiers() & CLEAR_MODIFIER)

    @override
    def event(self, e: QEvent) -> bool:
        if e.type() == QEvent.Type.ShortcutOverride and isinstance(e, QKeyEvent) and self._is_clear(e):
            e.accept()   # Control-L clears here, before any menu shortcut sees it
            return True
        if e.type() == QEvent.Type.KeyPress and isinstance(e, QKeyEvent) and e.key() == Qt.Key.Key_Tab:
            self.completion.emit()   # Tab completes instead of moving the focus
            return True
        return super().event(e)

    @override
    def keyPressEvent(self, e: QKeyEvent) -> None:
        if e.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            line = self.text()
            self.clear()
            self.submitted.emit(line)
        elif e.key() == Qt.Key.Key_Up:
            self._browse(-1)
        elif e.key() == Qt.Key.Key_Down:
            self._browse(1)
        elif self._is_clear(e):
            self.cleared.emit()
        elif e.matches(QKeySequence.StandardKey.Paste):
            text = QApplication.clipboard().text()
            if "\n" in text.rstrip("\n"):   # a block: run it line by line, as a terminal would
                lines = (self.text() + text).split("\n")
                self.clear()
                for line in lines:
                    self.submitted.emit(line)
            else:
                super().keyPressEvent(e)
        else:
            super().keyPressEvent(e)


class _Grip(QWidget):
    """Strip along the top edge of the console; dragging it changes the console's height."""

    def __init__(self, console: "ConsolePanel") -> None:
        super().__init__(console)
        self.console = console
        self.setFixedHeight(7)
        self.setCursor(Qt.CursorShape.SizeVerCursor)
        self.setToolTip("Drag to resize the console")
        self._start: tuple[float, int] | None = None

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        color = QColor(self.console.theme.muted)
        color.setAlphaF(0.5)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(color)
        middle = self.width() / 2
        p.drawRoundedRect(int(middle - 18), 2, 36, 3, 1.5, 1.5)
        p.end()

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if e.button() == Qt.MouseButton.LeftButton:
            self._start = (e.globalPosition().y(), self.console.height())
            e.accept()

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        if self._start is not None:
            y, height = self._start
            self.console.set_console_height(round(height + y - e.globalPosition().y()))

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        if self._start is not None:
            self._start = None
            self.console.save_height()


class ConsolePanel(FloatingPanel):
    """Python console docked along the bottom of the window.

    It runs code on the interface thread with the semantics of Python's
    interactive prompt (:class:`code.InteractiveConsole`): statements and
    expressions, multi-line blocks with a ``...`` prompt, results echoed,
    output and tracebacks shown in the console. The namespace comes from
    :func:`miles.scripting.namespace`; ``doc`` and ``mol`` follow the
    current tab before every command, and the view is redrawn after it.

    The console spans the bottom of the window above the timeline; drag its
    top edge to change its height. Like the other panels it can be folded
    and moved. History (the last :data:`HISTORY_LIMIT` commands), height and
    whether it is open are kept in the settings once :meth:`attach` gives it
    the application.
    """

    closed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.theme: Theme = LIGHT
        self.app: App | None = None
        self.names: dict[str, Any] = {}
        self.wanted_height = DEFAULT_HEIGHT
        self._settings: QSettings | None = None
        self._interpreter: _Interpreter | None = None
        self._chunks: list[tuple[str, str]] = []   # output as (text, kind), to recolour on theme changes
        self._more = False
        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 0, 10, 8)
        outer.setSpacing(2)
        self.grip = _Grip(self)
        outer.addWidget(self.grip)
        top = QHBoxLayout()
        top.addWidget(self.chevron)
        top.addWidget(label("Python", "Title"))
        top.addStretch(1)
        clear = button("Link", "Clear", "Clear the output (Control-L)")
        clear.clicked.connect(self.clear)
        top.addWidget(clear)
        close = button("Link", "Close")
        close.clicked.connect(self.close_console)
        top.addWidget(close)
        outer.addLayout(top)
        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        font.setPointSizeF(max(font.pointSizeF(), 11.0))
        self.output = QPlainTextEdit()
        self.output.setObjectName("ConsoleOutput")
        self.output.setReadOnly(True)
        self.output.setMaximumBlockCount(OUTPUT_LINES)
        self.output.setFont(font)
        self.output.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.body_layout.setSpacing(4)
        self.body_layout.addWidget(self.output, 1)
        row = QHBoxLayout()
        row.setSpacing(6)
        self.prompt = label(PROMPT, "Muted")
        self.prompt.setFont(font)
        row.addWidget(self.prompt)
        self.input = ConsoleInput()
        self.input.setObjectName("ConsoleInput")
        self.input.setFont(font)
        self.input.setPlaceholderText("Python: try select(\"protein and name CA\") or show(\"tube\")")
        self.input.submitted.connect(self.execute)
        self.input.completion.connect(self.complete)
        self.input.cleared.connect(self.clear)
        row.addWidget(self.input, 1)
        self.body_layout.addLayout(row)
        outer.addWidget(self.body, 1)
        self.resize(600, self.wanted_height)
        self.hide()

    # -- setup ----------------------------------------------------------------------------------

    def attach(self, app: App, settings: QSettings | None = None) -> None:
        """Connect to the application; restore history and height from ``settings``."""
        self.app = app
        self.names = namespace(app)
        self._interpreter = _Interpreter(self.names, self)
        self._settings = settings
        if settings is not None:
            history: Any = settings.value("console/history", [])
            lines = [history] if isinstance(history, str) else list(history or [])
            self.input.set_history([str(line) for line in lines])
            height: Any = settings.value("console/height", DEFAULT_HEIGHT)
            try:
                self.wanted_height = max(MIN_HEIGHT, int(height))
            except (TypeError, ValueError):
                self.wanted_height = DEFAULT_HEIGHT
        self.write(BANNER, "info")

    def set_theme(self, theme: Theme) -> None:
        """Recolour the console and the text already in it."""
        self.theme = theme
        self.set_icon_color(theme.muted)
        border = QColor(theme.text)
        border.setAlphaF(0.14)
        field = QColor(theme.text)
        field.setAlphaF(0.04)
        self.output.setStyleSheet(
            f"QPlainTextEdit {{ background: transparent; border: none; color: {theme.text.name()};"
            f" selection-background-color: {theme.accent.name()}; }}")
        self.input.setStyleSheet(
            f"QLineEdit {{ border: 1px solid {border.name(QColor.NameFormat.HexArgb)}; border-radius: 6px;"
            f" padding: 3px 6px; color: {theme.text.name()};"
            f" background: {field.name(QColor.NameFormat.HexArgb)};"
            f" selection-background-color: {theme.accent.name()}; }}")
        chunks, self._chunks = self._chunks, []
        self.output.clear()
        for text, kind in chunks:
            self.write(text, kind)
        self.grip.update()

    def set_open(self, shown: bool, remember: bool = True) -> None:
        """Show the console with the input focused, or hide it.

        ``remember`` keeps the choice for the next start; only the user's own
        choice is kept, not a console opened to run a script.
        """
        self.setVisible(shown)
        if shown:
            self.input.setFocus()
        if remember and self._settings is not None:
            self._settings.setValue("console/open", "true" if shown else "false")
        self.relayout.emit()

    def close_console(self) -> None:
        """Hide the console (the Close link)."""
        self.set_open(False)
        self.closed.emit()

    # -- size -----------------------------------------------------------------------------------

    @override
    def fit(self) -> None:
        if self.collapsed:
            self.resize(self.width(), self.minimumSizeHint().height())
        else:
            self.resize(self.width(), self.wanted_height)

    @override
    def set_collapsed(self, collapsed: bool) -> None:
        self.grip.setVisible(not collapsed)
        super().set_collapsed(collapsed)

    def set_console_height(self, height: int) -> None:
        """Change the height, keeping the bottom edge in place."""
        parent = self.parentWidget()
        limit = parent.height() - 60 if parent is not None else 2000
        height = max(MIN_HEIGHT, min(height, limit))
        bottom = self.geometry().bottom()
        self.wanted_height = height
        self.fit()
        if self.placed is not None:
            self.placed = QPoint(self.placed.x(), bottom - height + 1)
        self.relayout.emit()

    def save_height(self) -> None:
        """Remember the height for the next session."""
        if self._settings is not None:
            self._settings.setValue("console/height", self.wanted_height)
        self.arranged.emit()

    # -- output ---------------------------------------------------------------------------------

    def write(self, text: str, kind: str = "output") -> None:
        """Add text to the output area.

        Parameters
        ----------
        text : str
            Text to add, with its own newlines.
        kind : {"output", "error", "input", "info"}, optional
            Sets the colour.
        """
        if not text:
            return
        self._chunks.append((text, kind))
        del self._chunks[:-OUTPUT_LINES]
        colors = {"output": self.theme.text, "error": self.theme.danger, "input": self.theme.muted,
                  "info": self.theme.muted}
        fmt = QTextCharFormat()
        fmt.setForeground(colors.get(kind, self.theme.text))
        cursor = self.output.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        cursor.insertText(text, fmt)
        self.output.setTextCursor(cursor)
        self.output.ensureCursorVisible()

    def clear(self) -> None:
        """Empty the output area."""
        self._chunks = []
        self.output.clear()

    def text(self) -> str:
        """Everything in the output area."""
        return self.output.toPlainText()

    @contextlib.contextmanager
    def _capturing(self) -> Iterator[None]:
        """Send ``print`` output and errors to the console while code runs; redraw afterwards."""
        if self.app is None:
            raise RuntimeError("The console is not attached to the application")
        refresh(self.names, self.app)
        before = snapshot(self.app)
        out, err = _Stream(self, "output"), _Stream(self, "error")
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                yield
            finally:
                settle(self.app, before)
                refresh(self.names, self.app)

    # -- running code ---------------------------------------------------------------------------

    def execute(self, line: str) -> bool:
        """Run one line typed at the prompt; return whether more lines are expected.

        Parameters
        ----------
        line : str
            The line, as typed.
        """
        if self._interpreter is None:
            return False
        self.write(f"{MORE if self._more else PROMPT} {line}\n", "input")
        self.input.remember(line)
        self._save_history()
        with self._capturing():
            try:
                self._more = self._interpreter.push(line)
            except SystemExit:
                self._more = False
                self.write("The console stays open; close the window to quit MILES.\n", "info")
        self.prompt.setText(MORE if self._more else PROMPT)
        return self._more

    def run_source(self, source: str, filename: str = "<script>") -> bool:
        """Run a whole program in the console's namespace; return whether it ran without error.

        Errors are shown as a traceback in the console and passed to the
        application's error handler (the error dialog).

        Parameters
        ----------
        source : str
            Python code.
        filename : str, optional
            Name in tracebacks.
        """
        if self.app is None:
            return False
        self.write(f"# running {filename}\n", "info")
        name = self.names.get("__name__")
        self.names["__name__"] = "__main__"
        try:
            with self._capturing():
                try:
                    exec(compile(source, filename, "exec"), self.names)
                except SystemExit:
                    self.write(f"# {os.path.basename(filename)} stopped\n", "info")
                except Exception:
                    kind, value, tb = sys.exc_info()
                    self._report(kind, value, tb)
                    return False
        finally:
            self.names["__name__"] = name
        return True

    def run_file(self, path: str) -> bool:
        """Run a Python file in the console's namespace; return whether it ran without error.

        ``__file__`` is the file's path and ``__name__`` is ``"__main__"``
        while it runs; what it defines stays available in the console.
        """
        with tokenize.open(path) as f:
            source = f.read()
        previous = self.names.get("__file__")
        self.names["__file__"] = os.path.abspath(path)
        try:
            return self.run_source(source, path)
        finally:
            if previous is None:
                self.names.pop("__file__", None)
            else:
                self.names["__file__"] = previous

    def _report(self, kind: type[BaseException] | None, value: BaseException | None,
                tb: TracebackType | None) -> None:
        if kind is None or value is None:
            return
        self.write("".join(traceback.format_exception(kind, value, tb.tb_next if tb else None)), "error")
        if sys.excepthook is not sys.__excepthook__:
            sys.excepthook(kind, value, tb)   # the error dialog

    def _save_history(self) -> None:
        if self._settings is not None:
            self._settings.setValue("console/history", self.input.history)

    # -- completion -----------------------------------------------------------------------------

    def complete(self) -> None:
        """Complete the name before the cursor; with several candidates, list them."""
        text = self.input.text()[:self.input.cursorPosition()]
        token = re.search(r"[A-Za-z_][\w.]*$", text)
        if token is None:
            return
        word = token.group(0)
        completer = rlcompleter.Completer(self.names)
        matches: list[str] = []
        while len(matches) < 200:
            match = completer.complete(word, len(matches))
            if match is None:
                break
            matches.append(match)
        matches = sorted(set(matches))
        if not matches:
            return
        common = matches[0]
        for match in matches[1:]:   # the longest start all candidates share
            while not match.startswith(common):
                common = common[:-1]
        if len(common) > len(word):
            self.input.insert(common[len(word):])
        elif len(matches) > 1:
            self.write("  ".join(m.rsplit(".", 1)[-1] for m in matches) + "\n", "info")
