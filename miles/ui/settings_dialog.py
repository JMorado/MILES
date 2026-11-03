"""Settings: the appearance, autosave, how new documents are drawn, and the settings of plugins.

Every change applies at once and is kept with the other settings, so the
window has no OK or Cancel button.
"""

import contextlib
from collections.abc import Callable
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractSpinBox,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from miles.core.molecule import Representation
from miles.plugins import PluginState, Setting
from miles.render.cartoon import CartoonColor
from miles.render.scene import REPRESENTATION_NAMES, LabelMode, Style
from miles.render.theme import Appearance, Theme
from miles.ui.autosave import INTERVALS, interval_text
from miles.ui.controls import dialog_stylesheet, label, theme_palette

if TYPE_CHECKING:
    from miles.ui.mainwindow.window import MainWindow

FIELD_WIDTH = 220   # every field has this width, so the column lines up
TEXT_INSET = 7      # pixels before the text of a number field, as in the drop-down lists
IMAGES = Path(__file__).parent / "images"
DEFAULTS = "defaults"
"""Settings group of the display options new documents start with."""

APPEARANCES = [(Appearance.SYSTEM, "Match System"), (Appearance.LIGHT, "Light"), (Appearance.DARK, "Dark")]
STYLES = [(rep, REPRESENTATION_NAMES[rep]) for rep in REPRESENTATION_NAMES
          if rep not in (Representation.CARTOON, Representation.HIDDEN)]   # cartoons only suit polymers
LABELS = [(LabelMode.NONE, "None"), (LabelMode.SYMBOL, "Element"), (LabelMode.NUMBER, "Number"),
          (LabelMode.BOTH, "Element and Number")]
CARTOON_COLORS = [(CartoonColor.SECONDARY_STRUCTURE, "Secondary Structure"), (CartoonColor.CHAIN, "Chain"),
                  (CartoonColor.RAINBOW, "Rainbow")]

_CHOICES: dict[str, tuple[type[StrEnum], list[tuple[Any, str]]]] = {
    "mode": (Representation, STYLES),
    "labels": (LabelMode, LABELS),
    "cartoon_color": (CartoonColor, CARTOON_COLORS),
}
_SWITCHES = ("show_hydrogens", "depth_fade")


def new_document_style(settings: QSettings) -> dict[str, Any]:
    """Display options chosen in Settings for new documents, as keyword arguments of ``set_style``.

    Only options the user changed are returned, each checked against the
    choices offered, so a document without any keeps the built-in defaults.
    """
    style: dict[str, Any] = {}
    for name, (kind, options) in _CHOICES.items():
        raw = settings.value(f"{DEFAULTS}/{name}", None)
        if raw is None:
            continue
        with contextlib.suppress(ValueError):
            value = kind(str(raw))
            if any(value == option for option, _ in options):
                style[name] = value
    for name in _SWITCHES:
        raw = settings.value(f"{DEFAULTS}/{name}", None)
        if raw is not None:
            style[name] = str(raw).lower() in ("true", "1")
    return style


def setting_editor(setting: Setting) -> QWidget:
    """Editor for a plugin's setting: a check box, a number field or a text field, saving each change."""
    value = setting.value
    if isinstance(setting.default, bool):
        box = QCheckBox(setting.title)
        box.setChecked(bool(value))
        box.toggled.connect(setting.set)
        return box
    if isinstance(setting.default, (int, float)):
        field: QSpinBox | QDoubleSpinBox
        if isinstance(setting.default, int):
            field = QSpinBox()
            field.setRange(int(setting.minimum if setting.minimum is not None else -10**9),
                           int(setting.maximum if setting.maximum is not None else 10**9))
            field.setSingleStep(int(setting.step or 1))
        else:
            field = QDoubleSpinBox()
            step = setting.step or 0.1
            field.setDecimals(max(1, len(f"{step:g}".partition(".")[2])))
            field.setRange(setting.minimum if setting.minimum is not None else -1e9,
                           setting.maximum if setting.maximum is not None else 1e9)
            field.setSingleStep(step)
        field.setValue(value)
        if setting.unit:
            field.setSuffix(f" {setting.unit}")
        field.setKeyboardTracking(False)   # one change when typing ends, not one per key
        field.valueChanged.connect(setting.set)
        return field
    line = QLineEdit(str(value))
    line.setObjectName("Field")
    line.editingFinished.connect(lambda: setting.set(line.text()))
    return line


def _rgba(color: QColor, alpha: float) -> str:
    return f"rgba({color.red()}, {color.green()}, {color.blue()}, {alpha:.3f})"


def fields_stylesheet(theme: Theme) -> str:
    """Drop-down lists, number fields and text fields drawn alike, so they line up in one column.

    Native drop-down lists are drawn inside their frame by a margin that
    differs between platforms, while styled fields fill it; styling them all
    the same way keeps the column straight and readable in either theme.
    """
    field = _rgba(theme.text, 0.18)
    down = (IMAGES / "chevron_down.svg").as_posix()
    up = (IMAGES / "chevron_up.svg").as_posix()
    return f"""
    QScrollArea, QWidget#SettingsBody {{ background: {theme.panel.name()}; }}
    QDialog QComboBox, QDialog QAbstractSpinBox, QDialog QLineEdit#Field {{
        border: 1px solid {field}; border-radius: 6px; padding: 4px 24px 4px 8px; min-height: 18px;
        background: {theme.panel.name()}; color: {theme.text.name()};
        selection-background-color: {theme.accent.name()}; }}
    QDialog QLineEdit#Field {{ padding-right: 8px; }}
    QDialog QComboBox:hover, QDialog QAbstractSpinBox:hover {{ border-color: {_rgba(theme.text, 0.32)}; }}
    QDialog QComboBox::drop-down {{ subcontrol-origin: padding; subcontrol-position: center right;
        width: 22px; border: none; }}
    QDialog QComboBox::down-arrow {{ image: url({down}); width: 10px; height: 10px; }}
    QDialog QComboBox QAbstractItemView {{ background: {theme.panel.name()}; color: {theme.text.name()};
        border: 1px solid {field}; outline: none; padding: 4px;
        selection-background-color: {_rgba(theme.accent, 0.28)}; selection-color: {theme.text.name()}; }}
    QDialog QAbstractSpinBox::up-button, QDialog QAbstractSpinBox::down-button {{
        subcontrol-origin: border; width: 22px; border: none; background: transparent; }}
    QDialog QAbstractSpinBox::up-button {{ subcontrol-position: top right; }}
    QDialog QAbstractSpinBox::down-button {{ subcontrol-position: bottom right; }}
    QDialog QAbstractSpinBox::up-arrow {{ image: url({up}); width: 10px; height: 10px; }}
    QDialog QAbstractSpinBox::down-arrow {{ image: url({down}); width: 10px; height: 10px; }}
    """


class SettingsDialog(QDialog):
    """The Settings window: one form, in sections, whose changes apply at once.

    Parameters
    ----------
    window : MainWindow
        The window whose settings are shown and changed.
    """

    def __init__(self, window: "MainWindow") -> None:
        super().__init__(window)
        self.main = window
        self.setWindowTitle("Settings")
        self.sections: dict[str, QLabel] = {}   # section titles by plugin key, to scroll to
        self._combos: list[QComboBox] = []
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 14)
        lay.setSpacing(10)
        self.area = QScrollArea()
        self.area.setWidgetResizable(True)
        self.area.setFrameShape(QFrame.Shape.NoFrame)
        self.area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        lay.addWidget(self.area, 1)
        buttons = QHBoxLayout()
        buttons.setContentsMargins(22, 0, 22, 0)
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.setObjectName("Primary")
        close.setDefault(True)
        close.clicked.connect(self.close)
        buttons.addWidget(close)
        lay.addLayout(buttons)
        self.rebuild()
        self.set_theme(window.theme)
        window.themeChanged.connect(self.set_theme)
        self.fit()

    def fit(self) -> None:
        """Take the height of the whole form, as far as the screen allows."""
        body = self.area.widget()
        if body is None:
            return
        wanted = body.sizeHint().height() + 74   # and the row of buttons below
        screen = self.screen().availableGeometry().height() if self.screen() is not None else 800
        self.resize(560, min(wanted, int(screen * 0.85)))

    # -- building ---------------------------------------------------------------------------

    def rebuild(self) -> None:
        """Lay the form out again, e.g. after a plugin was switched on or off."""
        body = QWidget()
        body.setObjectName("SettingsBody")
        form = QFormLayout(body)
        form.setContentsMargins(22, 16, 22, 8)
        form.setHorizontalSpacing(14)
        form.setVerticalSpacing(9)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
        self.form = form
        self.sections.clear()
        self._combos.clear()
        settings = self.main.settings

        self._section("General")
        self.appearance = self._choice("Appearance:", APPEARANCES, self.main.appearance, self._set_appearance)
        autosave = self.main.autosave
        intervals = [*sorted({*INTERVALS, autosave.interval} - {0}), 0]   # one set by a script too
        self.autosave = self._choice("Autosave:", [(m, interval_text(m)) for m in intervals],
                                     autosave.interval, autosave.set_interval)
        self._note("Copies of documents with unsaved changes, offered again if MILES stops "
                   "before you save them.")

        self._section("New Documents")
        chosen = new_document_style(settings)
        builtin = Style()
        for name, title in (("mode", "Style:"), ("labels", "Labels:"), ("cartoon_color", "Cartoon colors:")):
            _, options = _CHOICES[name]
            current = chosen.get(name, getattr(builtin, name))
            self._choice(title, options, current, self._default_setter(name))
        for name, text in (("show_hydrogens", "Show hydrogens"), ("depth_fade", "Fade distant atoms")):
            box = QCheckBox(text)
            box.setChecked(bool(chosen.get(name, getattr(builtin, name))))
            box.toggled.connect(lambda on, n=name: self._set_default(n, on))
            form.addRow("", box)
        self._note("Used for documents opened from now on. The View menu changes the current document.")

        manager = getattr(self.main, "plugins", None)   # not there yet while the window starts
        for info in [] if manager is None else manager.plugins:
            if manager is None or info.state != PluginState.LOADED:
                continue
            plugin_settings = manager.settings.get(info.key, [])
            if not plugin_settings:
                continue
            self.sections[info.key] = self._section(info.name)
            for setting in plugin_settings:
                editor = setting_editor(setting)
                if isinstance(editor, QCheckBox):
                    form.addRow("", editor)   # its title is its text
                    continue
                editor.setFixedWidth(FIELD_WIDTH)
                inner = editor.findChild(QLineEdit) if isinstance(editor, QAbstractSpinBox) else None
                if inner is not None:   # spin boxes ignore horizontal padding on macOS
                    inner.setTextMargins(TEXT_INSET, 0, 0, 0)
                form.addRow(f"{setting.title}:", editor)
        self.area.setWidget(body)   # the old form is deleted with its widget
        if hasattr(self, "theme"):
            self.set_theme(self.theme)

    def _section(self, title: str) -> QLabel:
        heading = label(title, "Title")
        if self.form.rowCount():
            heading.setContentsMargins(0, 12, 0, 0)
        self.form.addRow(heading)
        return heading

    def _note(self, text: str) -> None:
        note = label(text, "Muted")
        note.setWordWrap(True)
        self.form.addRow("", note)

    def _choice(self, title: str, options: list[tuple[Any, str]], current: Any,
                chosen: Callable[[Any], None]) -> QComboBox:
        combo = QComboBox()
        combo.setFixedWidth(FIELD_WIDTH)
        for _, text in options:
            combo.addItem(text)
        index = next((k for k, (value, _) in enumerate(options) if value == current), 0)
        combo.setCurrentIndex(index)
        combo.currentIndexChanged.connect(lambda k: chosen(options[k][0]))
        self._combos.append(combo)
        self.form.addRow(title, combo)
        return combo

    # -- changes ----------------------------------------------------------------------------

    def _set_appearance(self, appearance: str) -> None:
        if appearance != self.main.appearance:
            self.main.set_appearance(appearance)

    def _default_setter(self, name: str) -> Callable[[Any], None]:
        return lambda value: self._set_default(name, value)

    def _set_default(self, name: str, value: Any) -> None:
        self.main.settings.setValue(f"{DEFAULTS}/{name}", str(value).lower() if isinstance(value, bool)
                                    else str(value))

    def show_plugin(self, key: str) -> None:
        """Scroll to the settings of a plugin."""
        heading = self.sections.get(key)
        if heading is not None:
            QTimer.singleShot(0, lambda: self.area.ensureWidgetVisible(heading, 0, 40))

    def set_theme(self, theme: Theme) -> None:
        """Draw the window in the app's theme and show the appearance in use."""
        self.theme = theme
        self.setStyleSheet(dialog_stylesheet(theme) + fields_stylesheet(theme))
        for combo in self._combos:
            theme_palette(combo, theme)
        index = next((k for k, (value, _) in enumerate(APPEARANCES) if value == self.main.appearance), -1)
        if index >= 0 and index != self.appearance.currentIndex():
            self.appearance.blockSignals(True)
            self.appearance.setCurrentIndex(index)
            self.appearance.blockSignals(False)
