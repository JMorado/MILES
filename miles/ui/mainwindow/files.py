"""Opening and saving documents, recent files, exported pictures and data, and fragment libraries."""

import os
import sys
from typing import TYPE_CHECKING

from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QMenu,
    QMessageBox,
)

from miles.core.fragments import Fragment, read_library
from miles.core.molecule import Molecule
from miles.core.protein import default_representations
from miles.files import session
from miles.files.common import FileFormatError, JobInfo, extension
from miles.files.formats import read_file, read_filter, write_file, write_filter
from miles.ui.autosave import Leftover, read_leftover
from miles.ui.builder import Tool
from miles.ui.canvas.view_state import apply_view_state, view_state
from miles.ui.document import Document
from miles.ui.export_dialog import ExportDialog
from miles.ui.fragment_popover import fragment_folder
from miles.ui.plot_export import COPY, CURVE, IMAGE, PlotSource, save_plot_image
from miles.ui.text import SEPARATOR, plural

if TYPE_CHECKING:
    from miles.ui.mainwindow.window import MainWindow

STRUCTURE_EXTENSIONS = {"xyz", "mol", "sdf", "pdb", "cif", "gjf", "com", "log", "out", "miles"}
"""Extensions that make a command-line item a file name rather than a name or SMILES string."""

READ_ONLY_EXTENSIONS = {"log", "out"}
"""Program output that Save writes elsewhere instead of overwriting."""

RECENT_LIMIT = 10
"""Files kept in File > Open Recent."""


class FileController:
    """What a window reads and writes: documents, recent files, exports and fragment libraries.

    Parameters
    ----------
    window : MainWindow
        The window whose documents these are; dialogs open over it and new
        documents become its tabs.
    """

    def __init__(self, window: "MainWindow") -> None:
        self.window = window
        self.settings = window.settings
        self._recent_menu: QMenu | None = None   # File > Open Recent, once the menus fill it

    def last_dir(self) -> str:
        """Folder of the file opened or saved last, where file dialogs start."""
        return str(self.settings.value("last_dir", os.path.expanduser("~")))

    def remember_dir(self, path: str, folder: bool = False) -> None:
        """Start the next file dialog in the folder of ``path``, or in ``path`` itself if it is a folder."""
        self.settings.setValue("last_dir", os.path.abspath(path) if folder else
                               os.path.dirname(os.path.abspath(path)))

    def open_dialog(self) -> None:
        """Ask for files and open them."""
        paths, _ = QFileDialog.getOpenFileNames(self.window, "Open", self.last_dir(), read_filter())
        for path in paths:
            self.open_path(path)

    def export_plot(self, source: PlotSource, kind: str) -> None:
        """Copy a plot's data, or ask where to save it (CSV) or its picture (PNG or SVG).

        Parameters
        ----------
        source : PlotSource
            Panel whose plot was right-clicked.
        kind : str
            One of the kinds in :mod:`miles.ui.plot_export`.
        """
        window = self.window
        canvas = window.current_canvas()
        if kind == IMAGE:
            plot = source.plot_widget()
            if plot is None:
                return
            start = os.path.join(self.last_dir(), source.export_name() + ".png")
            path, _ = QFileDialog.getSaveFileName(window, "Export Image", start,
                                                  "PNG image (*.png);;SVG drawing (*.svg)")
            if not path:
                return
            try:
                save_plot_image(plot, path)
            except OSError as exc:
                QMessageBox.warning(window, "Could not export", str(exc))
                return
        else:
            table = source.export_table(kind)
            if table is None:
                return
            if kind == COPY:
                QApplication.clipboard().setText(table.tsv())
                if canvas is not None:
                    canvas.toast(f"Copied {len(table.rows)} rows for a spreadsheet", 2000)
                return
            suffix = "_spectrum" if kind == CURVE else ""
            start = os.path.join(self.last_dir(), source.export_name() + suffix + ".csv")
            path, _ = QFileDialog.getSaveFileName(window, "Export Data", start, "CSV files (*.csv)")
            if not path:
                return
            try:
                table.write_csv(path)
            except OSError as exc:
                QMessageBox.warning(window, "Could not export", str(exc))
                return
        self.remember_dir(path)
        if canvas is not None:
            canvas.toast(f"Exported {os.path.basename(path)}", 2000)

    def open_path(self, path: str) -> None:
        """Open a structure file in a new tab."""
        try:
            mol, info = read_file(path)
        except (OSError, FileFormatError, ValueError, IndexError) as exc:
            QMessageBox.warning(self.window, "Could not open file", f"{os.path.basename(path)}\n\n{exc}")
            return
        self.remember_dir(path)
        self._add_recent(path)
        self.open_structure(mol, info, path)

    def recent_files(self) -> list[str]:
        """Files opened or saved lately, most recent first."""
        value = self.settings.value("recent_files", [])
        if isinstance(value, str):   # QSettings gives a single entry back as a plain string
            value = [value]
        return [str(path) for path in value if path] if isinstance(value, list) else []

    def _set_recent(self, paths: list[str]) -> None:
        self.settings.setValue("recent_files", paths[:RECENT_LIMIT])
        if self._recent_menu is not None:
            self.fill_recent_menu(self._recent_menu)

    def _add_recent(self, path: str) -> None:
        path = os.path.abspath(path)
        self._set_recent([path, *(p for p in self.recent_files() if p != path)])

    def fill_recent_menu(self, menu: QMenu) -> None:
        """List the recent files in File > Open Recent."""
        self._recent_menu = menu
        menu.clear()
        paths = self.recent_files()
        names = [os.path.basename(p) for p in paths]
        for path, name in zip(paths, names, strict=True):
            folder = os.path.basename(os.path.dirname(path))
            text = f"{name}   ({folder})" if names.count(name) > 1 and folder else name
            action = menu.addAction(text)
            action.setToolTip(path)
            action.triggered.connect(lambda _=False, p=path: self.open_recent(p))
        if not paths:
            menu.addAction("No Recent Files").setEnabled(False)
        menu.addSeparator()
        clear = menu.addAction("Clear Menu")
        clear.setEnabled(bool(paths))
        clear.triggered.connect(self.clear_recent)

    def open_recent(self, path: str) -> None:
        """Open a file from File > Open Recent, forgetting it if it is gone."""
        if not os.path.exists(path):
            self._set_recent([p for p in self.recent_files() if p != path])
            QMessageBox.warning(self.window, "File not found",
                                f"{path}\n\nIt was moved or deleted, so it was removed from the list.")
            return
        self.open_path(path)

    def clear_recent(self) -> None:
        """Empty File > Open Recent."""
        self._set_recent([])

    def open_item(self, item: str) -> None:
        """Open a file, or download or build a PDB code, name or SMILES string (command line)."""
        if os.path.exists(item):
            self.open_path(item)
        elif extension(item) in STRUCTURE_EXTENSIONS:
            QMessageBox.warning(self.window, "File not found", item)
        else:
            self.window.lookup.start(item)

    def open_structure(self, mol: Molecule, info: JobInfo, path: str | None = None) -> None:
        """Show a structure in a new tab, replacing an empty untitled one."""
        window = self.window
        reps = default_representations(mol)
        polymer = reps is not None
        current = window.current_canvas()
        replace = current is not None and not len(current.doc.mol) and not current.doc.modified
        view = info.view
        saved_display = info.display is not None   # a session keeps its own representations
        doc = Document(mol, info, path)
        if reps is not None and not saved_display:
            doc.display = doc.display.with_rep(reps)
        canvas = window.add_document(doc)
        canvas.view.fit()
        if view:
            apply_view_state(canvas, view)
            window.sync_view_menu(canvas)
        if replace and current is not None:
            window.close_canvas(current)
        if polymer:
            window.builder.set_tool(Tool.SELECT)
        name = os.path.basename(path) if path else info.title
        extra = ""
        if info.vibrations:
            extra = f"{SEPARATOR}{len(info.vibrations)} vibrational modes"
        elif info.orbitals is not None:
            extra = f"{SEPARATOR}{plural(len(info.orbitals.orbitals), 'orbital')}"
        elif info.volumes:
            extra = f"{SEPARATOR}{plural(len(info.volumes), 'grid')}"
        message = f"Opened {name}{SEPARATOR}{len(mol)} atoms{extra}"
        if info.notes:   # what the reader had to leave out, e.g. orbitals in an unknown convention
            canvas.toast(f"{message}. {' '.join(info.notes)}", 8000)
        else:
            canvas.toast(message, 2600)

    def open_found(self, mol: Molecule, attach: int | None, smiles: str, text: str) -> None:
        """Open a molecule found by name or SMILES on the command line."""
        self.open_structure(mol, JobInfo(title=text))

    def load_fragment_libraries(self) -> None:
        """Add the fragment libraries of the user's fragment folder to the fragment window."""
        folder = fragment_folder()
        if not folder.is_dir():
            return
        for path in sorted(folder.glob("*.json")):
            try:
                self.add_fragment_library(str(path))
            except (OSError, ValueError) as exc:
                print(f"Fragment library {path} could not be read: {exc}", file=sys.stderr)

    def add_fragment_library(self, path: str) -> list[tuple[str, list[Fragment]]]:
        """Add the fragments of a library file (the format of ``data/fragments.json``) to the fragment window.

        Returns
        -------
        list of (str, list of Fragment)
            What was added to each category, for :meth:`remove_fragments`.

        Raises
        ------
        OSError, ValueError
            If the file cannot be read or is not a library.
        """
        added = list(read_library(path).items())
        for category, fragments in added:
            self.window.fragment_popover.add_fragments(category, fragments)
        return added

    def remove_fragments(self, added: list[tuple[str, list[Fragment]]]) -> None:
        """Take fragments that :meth:`add_fragment_library` added out of the fragment window again."""
        for category, fragments in added:
            self.window.fragment_popover.remove_fragments(category, fragments)

    def save(self) -> None:
        """Save the current document, asking for a name if it has none."""
        doc = self.window.current_doc()
        if doc is not None:
            self.save_document(doc)

    def save_as(self) -> None:
        """Save the current document under a new name."""
        doc = self.window.current_doc()
        if doc is not None:
            self.save_document(doc, ask=True)

    def save_document(self, doc: Document, ask: bool = False) -> bool:
        """Write a document to its file (or a chosen one); returns whether it was saved."""
        path = doc.path
        if ask or not path or extension(path) in READ_ONLY_EXTENSIONS:
            base = os.path.splitext(doc.path)[0] if doc.path else os.path.join(self.last_dir(), "molecule")
            path, _ = QFileDialog.getSaveFileName(self.window, "Save", base + ".xyz", write_filter())
            if not path:
                return False
            if not os.path.splitext(path)[1]:
                path += ".xyz"
        doc.settle()   # the document's own coordinates, not a blend shown during playback
        info = doc.info
        canvas = next((c for c in self.window.canvases if c.doc is doc), None)
        if extension(path) == session.EXTENSION and canvas is not None:
            info = doc.session_info(view_state(canvas))
        try:
            write_file(path, doc.mol, info)
        except (OSError, FileFormatError) as exc:
            QMessageBox.warning(self.window, "Could not save", str(exc))
            return False
        self.remember_dir(path)
        self._add_recent(path)
        doc.set_path(path)
        canvas = self.window.current_canvas()
        if canvas is not None and canvas.doc is doc:
            canvas.toast(f"Saved {os.path.basename(path)}")
        return True

    def export_image(self) -> None:
        """Open the publication picture dialog."""
        window = self.window
        canvas = window.current_canvas()
        if canvas is None or not len(canvas.mol):
            return
        dialog = ExportDialog(canvas, window.theme, self.settings, window)
        if dialog.exec() and dialog.saved_path:
            canvas.toast(f"Exported {os.path.basename(dialog.saved_path)}")

    def confirm_discard(self, doc: Document) -> bool:
        """Ask to save a changed document; returns False if the user cancelled."""
        if not doc.modified or not len(doc.mol):
            return True
        box = QMessageBox(self.window)
        box.setText(f"Save changes to {doc.name}?")
        box.setInformativeText("Your changes will be lost if you don't save them.")
        buttons = QMessageBox.StandardButton
        box.setStandardButtons(buttons.Save | buttons.Discard | buttons.Cancel)
        box.setDefaultButton(buttons.Save)
        answer = box.exec()
        if answer == buttons.Save:
            return self.save_document(doc)
        return answer == buttons.Discard

    def restore_documents(self, leftovers: list[Leftover]) -> int:
        """Open autosaved documents as modified tabs, then delete the autosaved copies.

        Each document gets back the file it came from, if any, so Save
        writes there; it is marked as changed, so it still has to be saved.

        Returns
        -------
        int
            Number of documents restored.
        """
        window = self.window
        restored = 0
        problems: list[str] = []
        unreadable: list[Leftover] = []
        for leftover in reversed(leftovers):   # the newest ends up in front
            try:
                mol, info = read_leftover(leftover)
            except (OSError, FileFormatError, ValueError, KeyError) as exc:
                problems.append(f"{leftover.name}: {exc}")
                unreadable.append(leftover)
                continue
            self.open_structure(mol, info, leftover.path)
            canvas = window.current_canvas()
            if canvas is None:
                continue
            canvas.doc.modified = True
            canvas.doc.renamed.emit()   # the tab and the title show the unsaved changes
            window.autosave.mark_changed(canvas.doc)
            restored += 1
        window.autosave.discard_leftovers(keep=unreadable)   # what could not be read is offered again
        if problems:
            QMessageBox.warning(window, "Could not restore", "\n\n".join(problems))
        window.autosave.save_now()   # the copies now belong to this window
        if restored:
            window.toast(f"Restored {plural(restored, 'document')}", 2600)
        return restored
