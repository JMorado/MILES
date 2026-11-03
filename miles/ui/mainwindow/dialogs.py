"""The popovers and dialogs a window opens, and the offer to restore documents after a crash."""

from typing import TYPE_CHECKING

from PySide6.QtCore import QPoint, QTimer
from PySide6.QtWidgets import QDialog, QWidget

from miles.ui.builder import Tool
from miles.ui.dialogs import AtomListDialog, ShortcutsDialog
from miles.ui.document import Document
from miles.ui.recovery_dialog import DISCARD, RecoveryDialog
from miles.ui.settings_dialog import SettingsDialog

if TYPE_CHECKING:
    from miles.ui.mainwindow.window import MainWindow


class WindowDialogs:
    """Fragments, command palette, symmetry, atom lists, shortcuts, plugins, settings and recovery.

    Parameters
    ----------
    window : MainWindow
        The window the dialogs open over.

    Attributes
    ----------
    settings_dialog : SettingsDialog or None
        The Settings window, once opened.
    recovery_dialog : RecoveryDialog or None
        The offer to restore documents, while it is open.
    """

    def __init__(self, window: "MainWindow") -> None:
        self.window = window
        self.settings_dialog: SettingsDialog | None = None
        self.recovery_dialog: RecoveryDialog | None = None
        self._atom_lists: dict[int, AtomListDialog] = {}

    def _popup_pos(self, anchor: QWidget) -> QPoint:
        return anchor.mapToGlobal(QPoint(anchor.width() + 10, -6))

    def show_fragments(self, search: bool = False) -> None:
        """Open the fragment popover next to the build button."""
        window = self.window
        window.builder.set_tool(Tool.BUILD)
        window.fragment_popover.show_at(self._popup_pos(window.rail.build))
        if search:
            window.fragment_popover.focus_search()

    def search_fragments(self) -> None:
        """Open the fragment popover with the cursor in the search box."""
        self.show_fragments(search=True)

    def look_up(self, text: str) -> None:
        """Look a text up in the fragment window, as if typed in its search box.

        A whole structure (e.g. a PDB entry) opens in a new tab; a molecule
        becomes the fragment to place.
        """
        self.show_fragments(search=True)
        self.window.fragment_popover.search_for(text)

    def show_command_palette(self) -> None:
        """Find and run any menu command by typing part of its name."""
        window = self.window
        top = window.mapToGlobal(QPoint((window.width() - window.palette_popover.width()) // 2, 70))
        window.palette_popover.open_at(top)

    def show_symmetry(self) -> None:
        """Open the point group popover."""
        doc = self.window.current_doc()
        if doc is not None:
            self.window.symmetry_popover.show_for(doc, self._popup_pos(self.window.rail.symmetry))

    def show_atom_list(self) -> None:
        """Open (or raise) the atom table of the current document."""
        doc = self.window.current_doc()
        if doc is None:
            return
        dialog = self._atom_lists.get(id(doc))
        if dialog is None:
            dialog = self._atom_lists[id(doc)] = AtomListDialog(doc, self.window, self.window.theme)
        dialog.show()
        dialog.raise_()

    def close_atom_list(self, doc: Document) -> None:
        """Close the atom table of a document whose tab closes."""
        dialog = self._atom_lists.pop(id(doc), None)
        if dialog is not None:
            dialog.close()

    def set_theme(self) -> None:
        """Recolour the open atom tables."""
        for dialog in self._atom_lists.values():
            dialog.set_theme(self.window.theme)

    def show_shortcuts(self) -> None:
        """Show the shortcut reference."""
        ShortcutsDialog(self.window, self.window.theme, self.window.menuBar()).exec()

    def show_plugins(self) -> None:
        """List the plugins, what they added and any errors; switch them on or off."""
        from miles.ui.plugins_dialog import PluginsDialog

        window = self.window
        dialog = PluginsDialog(window.plugins, window.theme, window)
        dialog.settingsRequested.connect(
            lambda key: QTimer.singleShot(0, window, lambda: self.show_settings(key)))   # once it has closed
        dialog.exec()
        if self.settings_dialog is not None and self.settings_dialog.isVisible():
            self.settings_dialog.rebuild()   # plugins may have been switched on or off

    def show_settings(self, plugin: str = "") -> None:
        """Open the Settings window, scrolled to a plugin's settings if one is named by its key."""
        if self.settings_dialog is None:
            self.settings_dialog = SettingsDialog(self.window)
        else:
            self.settings_dialog.rebuild()
        self.settings_dialog.show()
        self.settings_dialog.raise_()
        self.settings_dialog.activateWindow()
        if plugin:
            self.settings_dialog.show_plugin(plugin)

    def offer_recovery(self) -> None:
        """Offer the documents autosaved by a MILES that stopped before they were saved."""
        if self.recovery_dialog is not None:
            return
        window = self.window
        leftovers = window.autosave.find_leftovers()
        if not leftovers:
            return
        dialog = RecoveryDialog(leftovers, window.theme, window)
        dialog.finished.connect(lambda result: self._recovery_answered(dialog, result))
        self.recovery_dialog = dialog
        dialog.open()

    def _recovery_answered(self, dialog: RecoveryDialog, result: int) -> None:
        self.recovery_dialog = None
        autosave = self.window.autosave
        if result == QDialog.DialogCode.Accepted:
            self.window.files.restore_documents(dialog.leftovers)
        elif result == DISCARD:
            autosave.discard_leftovers()
        else:
            autosave.release_leftovers()   # asked again at the next start
        dialog.deleteLater()
