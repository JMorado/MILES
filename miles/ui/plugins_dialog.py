"""Help > Plugins: the plugins found, what they added, their errors, and switching them on or off.

The settings of plugins are edited in the Settings window
(:mod:`miles.ui.settings_dialog`); the Settings button here opens it
at the selected plugin.
"""

from PySide6.QtCore import Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFontDatabase
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QPlainTextEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from miles.plugins import PluginInfo, PluginManager, PluginState, plugin_folder
from miles.render.theme import Theme
from miles.ui.controls import dialog_stylesheet, label

COLUMNS = ["On", "Plugin", "Version", "From", "State"]


def plugin_details(info: PluginInfo) -> str:
    """Plain-text description of a plugin, its commands, its warnings and its last error."""
    lines = [info.name + (f" {info.version}" if info.version else "")]
    if info.description:
        lines.append(info.description)
    lines += ["", f"From: {info.source}"]
    if info.path:
        lines.append(f"File: {info.path}")
    lines.append(f"State: {info.state}")
    if info.commands:
        lines += ["", "Commands:", *(f"  {path}" for path in info.commands)]
    if info.extras:
        lines += ["", "Also adds:", *(f"  {what}" for what in info.extras)]
    if info.warnings:
        lines += ["", "Warnings:", *(f"  {line}" for line in info.warnings.splitlines() if line.strip())]
    if info.error:
        lines += ["", "Error:", info.error.rstrip()]
    return "\n".join(lines)


class PluginsDialog(QDialog):
    """Table of plugins with a switch for each and the details of the one selected.

    Parameters
    ----------
    manager : PluginManager
        The application's plugins.
    theme : Theme
        Colours.
    parent : QWidget, optional
        Window to centre on.
    """

    settingsRequested = Signal(str)   # key of the plugin whose settings to show; the dialog closes

    def __init__(self, manager: PluginManager, theme: Theme, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.manager = manager
        self.theme = theme
        self.setWindowTitle("Plugins")
        self.resize(720, 640)
        self.setStyleSheet(dialog_stylesheet(theme))
        lay = QVBoxLayout(self)
        lay.setContentsMargins(18, 16, 18, 16)
        lay.setSpacing(10)
        lay.addWidget(label("Plugins", "Title"))
        self.table = QTableWidget(len(manager.plugins), len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().hide()
        self.table.setShowGrid(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        disabled = manager.disabled
        for row, info in enumerate(manager.plugins):
            switch = QTableWidgetItem()
            flags = Qt.ItemFlag
            switch.setFlags(flags.ItemIsUserCheckable | flags.ItemIsEnabled | flags.ItemIsSelectable)
            switch.setCheckState(Qt.CheckState.Unchecked if info.key in disabled else Qt.CheckState.Checked)
            switch.setToolTip("Switch the plugin on or off, now and at later starts")
            self.table.setItem(row, 0, switch)
            for column in range(1, len(COLUMNS)):
                self.table.setItem(row, column, QTableWidgetItem())
            self._fill_row(row)
        self.table.itemChanged.connect(self._switched)
        self.table.currentCellChanged.connect(lambda row, *_: self._show(row))
        lay.addWidget(self.table, 1)
        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.details.setFixedHeight(150)
        lay.addWidget(self.details)
        note = label(f"Switching a plugin on or off takes effect at once. Reload reads a changed plugin "
                     f"file again. The settings of plugins are in Settings. Plugin files go in "
                     f"{plugin_folder()}", "Muted")
        note.setWordWrap(True)
        note.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(note)
        buttons = QHBoxLayout()
        folder = QPushButton("Open Plugin Folder")
        folder.clicked.connect(self.open_folder)
        buttons.addWidget(folder)
        self.reload_button = QPushButton("Reload")
        self.reload_button.setToolTip("Unload the selected plugin and load it again, reading its file anew")
        self.reload_button.clicked.connect(self.reload_selected)
        buttons.addWidget(self.reload_button)
        self.settings_button = QPushButton("Settings...")
        self.settings_button.setToolTip("Change the settings of the selected plugin")
        self.settings_button.clicked.connect(self.open_settings)
        buttons.addWidget(self.settings_button)
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.setObjectName("Primary")
        close.setDefault(True)
        close.clicked.connect(self.accept)
        buttons.addWidget(close)
        lay.addLayout(buttons)
        failed = next((r for r, info in enumerate(manager.plugins) if info.state == PluginState.FAILED), 0)
        if manager.plugins:
            self.table.selectRow(failed)
            self._show(failed)

    def _fill_row(self, row: int) -> None:
        """Write a plugin's name, version, source and state into its row."""
        info = self.manager.plugins[row]
        texts = (info.name, info.version, info.source, str(info.state))
        for column, text in enumerate(texts, start=1):
            item = self.table.item(row, column)
            if item is None:
                continue
            item.setText(text)
            item.setToolTip(info.description if column == 1 else info.path if column == 3 else "")
            if column == 4:
                item.setForeground(self.theme.danger if info.state == PluginState.FAILED else self.theme.text)

    def _show(self, row: int) -> None:
        if not 0 <= row < len(self.manager.plugins):
            return
        info = self.manager.plugins[row]
        self.details.setPlainText(plugin_details(info))
        on = info.state != PluginState.SKIPPED and info.key not in self.manager.disabled
        self.reload_button.setEnabled(on)
        self.settings_button.setEnabled(self._has_settings(info))

    def _has_settings(self, info: PluginInfo) -> bool:
        return info.state == PluginState.LOADED and bool(self.manager.settings.get(info.key))

    def open_settings(self) -> None:
        """Close the dialog and ask for the Settings window at the selected plugin."""
        row = self.table.currentRow()
        if not 0 <= row < len(self.manager.plugins):
            return
        self.accept()
        self.settingsRequested.emit(self.manager.plugins[row].key)

    def _switched(self, item: QTableWidgetItem) -> None:
        if item.column() == 0:
            row = item.row()
            info = self.manager.plugins[row]
            self.table.blockSignals(True)
            self.manager.set_enabled(info.key, item.checkState() == Qt.CheckState.Checked)
            self._fill_row(row)
            self.table.blockSignals(False)
            if row == self.table.currentRow():
                self._show(row)

    def reload_selected(self) -> None:
        """Unload the selected plugin and load it again (a plugin file is read again)."""
        row = self.table.currentRow()
        if not 0 <= row < len(self.manager.plugins):
            return
        self.manager.reload(self.manager.plugins[row])
        self.table.blockSignals(True)
        self._fill_row(row)
        self.table.blockSignals(False)
        self._show(row)

    def open_folder(self) -> None:
        """Create the plugin folder if needed and show it in the file browser."""
        folder = plugin_folder()
        folder.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))
