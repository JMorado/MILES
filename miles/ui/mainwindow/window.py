"""The main window: its parts, the documents it shows, its panels, and closing it."""

import os
from collections.abc import Callable
from typing import override

from PySide6.QtCore import QPoint, QSettings, QTimer, Signal
from PySide6.QtGui import QCloseEvent, QDragEnterEvent, QDropEvent, QGuiApplication
from PySide6.QtWidgets import QFileDialog, QMainWindow, QTabBar, QVBoxLayout, QWidget

from miles.api import App
from miles.core.fragments import Fragment
from miles.core.molecule import Molecule
from miles.files.common import JobInfo
from miles.plugins import NO_PLUGINS_ENV, PluginManager
from miles.render.theme import Appearance, theme_for
from miles.ui import icons
from miles.ui.autosave import Autosave
from miles.ui.builder import Builder
from miles.ui.canvas import Canvas
from miles.ui.canvas.extensions import Extensions
from miles.ui.command_palette import CommandPalette
from miles.ui.controls import stylesheet
from miles.ui.document import Change, Document
from miles.ui.floating import FloatingPanel
from miles.ui.fragment_popover import FragmentPopover
from miles.ui.lookup import StructureLookup, StructureSource
from miles.ui.mainwindow.dialogs import WindowDialogs
from miles.ui.mainwindow.files import FileController
from miles.ui.mainwindow.tabs import DocumentTabs
from miles.ui.menus import MenuActions, build_menus
from miles.ui.plot_export import PlotSource
from miles.ui.settings_dialog import new_document_style
from miles.ui.stage import Stage
from miles.ui.symmetry_popover import SymmetryPopover
from miles.ui.text import BULLET

APP_NAME = "MILES"

CanvasAction = Callable[[Canvas], object]
Action = Callable[[], object]


class MainWindow(QMainWindow):
    """The application window. Each open document has a tab with its own canvas.

    The window puts its parts together and connects them:

    * :attr:`tabs` (:class:`~miles.ui.mainwindow.tabs.DocumentTabs`): the
      open documents in tab order and the one in front;
    * :attr:`files` (:class:`~miles.ui.mainwindow.files.FileController`):
      opening, saving, recent files, exports and fragment libraries;
    * :attr:`dialogs` (:class:`~miles.ui.mainwindow.dialogs.WindowDialogs`):
      the popovers and dialogs;
    * :attr:`stage` (:class:`~miles.ui.stage.Stage`): the canvases and the
      panels floating over them.

    The menus call these parts, or the canvas in front for commands on a
    document (:mod:`miles.ui.menus`). The methods of the window itself are
    what scripts, plugins and the other parts use: the documents, opening
    and saving, the theme and the Python console.

    Parameters
    ----------
    settings : QSettings, optional
        Where preferences are kept (default: the user's settings).
    plugins : bool, optional
        Load plugins of other packages and of the plugin folder as well as
        the built-in ones (False for ``--no-plugins``).

    Signals
    -------
    documentOpened(canvas)
        A tab was opened.
    documentClosed(canvas)
        A tab is about to close.
    documentSelectionChanged(canvas)
        The selection of a document changed.
    documentStructureChanged(canvas)
        The atoms, bonds or coordinates of a document changed.
    currentDocumentChanged(canvas)
        Another tab came to the front.
    documentChanged(canvas, kind)
        A document changed; ``kind`` is a :class:`miles.ui.document.Change`.
    themeChanged(theme)
        The light or dark theme changed.
    """

    documentOpened = Signal(object)
    documentClosed = Signal(object)
    documentSelectionChanged = Signal(object)
    documentStructureChanged = Signal(object)
    currentDocumentChanged = Signal(object)
    documentChanged = Signal(object, str)
    themeChanged = Signal(object)

    def __init__(self, settings: QSettings | None = None, plugins: bool = True) -> None:
        super().__init__()
        self.setWindowIcon(icons.app_icon())
        self.resize(1280, 820)
        self.setAcceptDrops(True)
        self.settings = settings or QSettings("miles", "miles")
        self.appearance = str(self.settings.value("appearance", Appearance.SYSTEM))
        self.theme = theme_for(self.appearance)
        self.builder = Builder()
        self.extensions = Extensions()   # what plugins draw and add to the canvases' menus
        self.structure_sources: list[StructureSource] = []   # added by plugins, shared by both lookups
        self.autosave = Autosave(lambda: self.canvases, self.settings)
        self.autosave.failed.connect(lambda message: self.toast(message, 4000))

        central = QWidget()
        lay = QVBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        self.tabbar = QTabBar()
        self.tabbar.setDocumentMode(True)
        self.tabbar.setTabsClosable(True)
        self.tabbar.setExpanding(False)
        self.tabbar.setDrawBase(False)
        self.tabbar.hide()
        lay.addWidget(self.tabbar)
        self.stage = Stage(self.builder)
        lay.addWidget(self.stage, 1)
        self.setCentralWidget(central)
        self.rail = self.stage.rail
        self.tabs = DocumentTabs(self.tabbar, self.stage.stack, self)
        self.tabs.currentChanged.connect(self._current_changed)
        self.tabbar.tabCloseRequested.connect(lambda index: self.close_canvas(self.tabbar.tabData(index)))

        self.files = FileController(self)
        self.dialogs = WindowDialogs(self)
        self.fragment_popover = FragmentPopover(self.builder, self, self.structure_sources)
        self.files.load_fragment_libraries()
        self.symmetry_popover = SymmetryPopover(self)
        self.lookup = StructureLookup(self.structure_sources)
        self.lookup.structure.connect(self.open_structure)
        self.lookup.fragment.connect(self.files.open_found)
        self.lookup.failed.connect(lambda message: self.toast(message, 4000))
        self.rail.fragmentRequested.connect(self.dialogs.show_fragments)
        self.rail.toolRequested.connect(self.builder.set_tool)
        self.rail.cleanRequested.connect(lambda: self._on_canvas(lambda c: c.commands.clean_geometry()))
        self.rail.symmetryRequested.connect(self.dialogs.show_symmetry)
        self.fragment_popover.placeRequested.connect(
            lambda: self._on_canvas(lambda c: c.commands.place_current_fragment()))
        self.fragment_popover.structureLoaded.connect(self.open_structure)
        self.builder.toolChanged.connect(self._tool_changed)
        self._restore_panels()
        for name, panel in self.stage.floating_panels().items():
            panel.arranged.connect(lambda n=name, p=panel: self._save_panel(n, p))
        self.stage.timeline.exportRequested.connect(lambda kind: self.export_plot(self.stage.timeline, kind))

        self.menu_actions: MenuActions = build_menus(self, self.menuBar())
        self.palette_popover = CommandPalette(self.menuBar(), self)
        self.menu_actions.check("tool", self.builder.tool)
        self.menu_actions.check("appearance", self.appearance)
        self.builder.mouseModeChanged.connect(
            lambda: self.menu_actions.check("mouse_mode", self.builder.mouse_mode))
        QGuiApplication.styleHints().colorSchemeChanged.connect(lambda *_: self.apply_theme())
        self.apply_theme()
        self.new_document()
        self.api: App = App(self)   # what scripts, plugins and the Python console use
        self.stage.console.attach(self.api, self.settings)
        console_toggle = self.menu_actions.toggles["python_console"]
        self.stage.console.closed.connect(lambda: console_toggle.setChecked(False))
        if str(self.settings.value("console/open", "false")).lower() in ("true", "1"):
            self.show_console(True, remember=False)   # already remembered
        self.plugins = PluginManager(self)
        self.plugins.load(third_party=plugins and not os.environ.get(NO_PLUGINS_ENV))
        QTimer.singleShot(0, self, self.dialogs.offer_recovery)   # once the window is up

    # -- documents ----------------------------------------------------------------------------

    @property
    def canvases(self) -> list[Canvas]:
        """Canvases of the open tabs, in tab order."""
        return self.tabs.canvases

    def current_canvas(self) -> Canvas | None:
        """Canvas of the visible tab."""
        return self.tabs.current

    def current_doc(self) -> Document | None:
        """Document of the visible tab."""
        canvas = self.tabs.current
        return canvas.doc if canvas is not None else None

    def _on_canvas(self, act: CanvasAction) -> None:
        canvas = self.tabs.current
        if canvas is not None:
            act(canvas)

    def toast(self, text: str, ms: int = 1600) -> None:
        """Show a short message at the top of the visible tab."""
        self._on_canvas(lambda canvas: canvas.toast(text, ms))

    def new_document(self) -> Canvas:
        """Open an empty tab."""
        return self.add_document(Document())

    def add_document(self, doc: Document) -> Canvas:
        """Open a tab for a document and bring it to the front."""
        canvas = Canvas(doc, self.builder, self.stage.surface, self.extensions, self.stage.tasks)
        canvas.set_theme(self.theme)
        style = new_document_style(self.settings)
        if style:
            canvas.set_style(**style)
        self.autosave.track(doc)
        canvas.selectionChanged.connect(lambda c=canvas: self._selection_changed(c))
        canvas.nudge.connect(
            lambda steps, c=canvas: self._if_current(c, lambda: self.stage.panel.nudge(steps)))
        canvas.playRequested.connect(lambda c=canvas: self._if_current(c, self.stage.timeline.toggle_playing))
        canvas.documentChanged.connect(lambda kind, c=canvas: self._doc_changed(c, kind))
        canvas.documentRenamed.connect(lambda c=canvas: self._refresh_title(c))
        self.tabs.add(canvas, doc.name)
        canvas.setFocus()
        self.documentOpened.emit(canvas)
        return canvas

    def close_canvas(self, canvas: Canvas) -> None:
        """Close a tab, asking to save changes first."""
        if canvas not in self.canvases or not self.files.confirm_discard(canvas.doc):
            return
        self.autosave.forget(canvas.doc)
        self.documentClosed.emit(canvas)
        self.dialogs.close_atom_list(canvas.doc)
        self.tabs.remove(canvas)
        if not self.canvases:
            self.new_document()

    def close_tab(self, index: int) -> None:
        """Close the tab at a position, asking to save changes first."""
        if 0 <= index < len(self.canvases):
            self.close_canvas(self.canvases[index])

    def close_current_tab(self) -> None:
        """Close the visible tab."""
        self._on_canvas(self.close_canvas)

    def _if_current(self, canvas: Canvas, act: Action) -> None:
        if canvas is self.tabs.current:
            act()

    def _current_changed(self, canvas: Canvas) -> None:
        self.stage.set_canvas(canvas)
        self._refresh_title(canvas)
        self.sync_view_menu(canvas)
        canvas.setFocus()
        self.currentDocumentChanged.emit(canvas)

    def _refresh_title(self, canvas: Canvas) -> None:
        doc = canvas.doc
        self.tabs.set_title(canvas, doc.name + (f" {BULLET}" if doc.modified else ""))
        if canvas is self.tabs.current:
            self.setWindowTitle(f"{doc.name}[*]")
            self.setWindowModified(doc.modified)
            self.setWindowFilePath(doc.path or "")
            self._update_undo()

    def _doc_changed(self, canvas: Canvas, kind: str) -> None:
        if kind in (Change.STRUCTURE, Change.BONDS, Change.COORDS):
            self.documentStructureChanged.emit(canvas)
        self.documentChanged.emit(canvas, kind)
        if canvas is self.tabs.current:
            self.stage.doc_changed(kind)
            self._update_undo()

    def _selection_changed(self, canvas: Canvas) -> None:
        self.documentSelectionChanged.emit(canvas)
        if canvas is self.tabs.current:
            self.stage.selection_changed()

    def _tool_changed(self) -> None:
        self.menu_actions.check("tool", self.builder.tool)
        self.stage.tool_changed()

    def _update_undo(self) -> None:
        doc = self.current_doc()
        self.menu_actions.commands["undo"].setEnabled(doc is not None and doc.can_undo())
        self.menu_actions.commands["redo"].setEnabled(doc is not None and doc.can_redo())

    def sync_view_menu(self, canvas: Canvas) -> None:
        """Check the View menu's items that match a canvas's display options."""
        display = canvas.display
        self.menu_actions.check("mode", display.mode)
        self.menu_actions.check("labels", display.labels)
        self.menu_actions.check("cartoon_color", display.cartoon_color)
        self.menu_actions.toggles["show_hydrogens"].setChecked(display.show_hydrogens)
        self.menu_actions.toggles["depth_fade"].setChecked(display.depth_fade)
        self.menu_actions.toggles["stereo"].setChecked(display.stereo)

    # -- files: what scripts, plugins and the command line use (see FileController) -----------

    def open_path(self, path: str) -> None:
        """Open a structure file in a new tab."""
        self.files.open_path(path)

    def open_item(self, item: str) -> None:
        """Open a file, or download or build a PDB code, name or SMILES string (command line)."""
        self.files.open_item(item)

    def open_structure(self, mol: Molecule, info: JobInfo, path: str | None = None) -> None:
        """Show a structure in a new tab, replacing an empty untitled one."""
        self.files.open_structure(mol, info, path)

    def save_document(self, doc: Document, ask: bool = False) -> bool:
        """Write a document to its file (or a chosen one); returns whether it was saved."""
        return self.files.save_document(doc, ask)

    def last_dir(self) -> str:
        """Folder of the file opened or saved last, where file dialogs start."""
        return self.files.last_dir()

    def remember_dir(self, path: str, folder: bool = False) -> None:
        """Start the next file dialog in the folder of ``path``, or in ``path`` itself if it is a folder."""
        self.files.remember_dir(path, folder)

    def export_plot(self, source: PlotSource, kind: str) -> None:
        """Copy a plot's data, or save it or its picture (see :meth:`FileController.export_plot`)."""
        self.files.export_plot(source, kind)

    def add_fragment_library(self, path: str) -> list[tuple[str, list[Fragment]]]:
        """Add the fragments of a library file to the fragment window (see :class:`FileController`)."""
        return self.files.add_fragment_library(path)

    def remove_fragments(self, added: list[tuple[str, list[Fragment]]]) -> None:
        """Take fragments that :meth:`add_fragment_library` added out of the fragment window again."""
        self.files.remove_fragments(added)

    def look_up(self, text: str) -> None:
        """Look a text up in the fragment window, as if typed in its search box."""
        self.dialogs.look_up(text)

    @override
    def dragEnterEvent(self, event: QDragEnterEvent) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    @override
    def dropEvent(self, event: QDropEvent) -> None:
        for url in event.mimeData().urls():
            if url.isLocalFile():
                self.files.open_path(url.toLocalFile())

    # -- appearance ---------------------------------------------------------------------------

    def set_appearance(self, appearance: str) -> None:
        """Use the system, light or dark appearance and remember the choice."""
        self.appearance = appearance
        self.settings.setValue("appearance", appearance)
        self.menu_actions.check("appearance", appearance)
        self.apply_theme()

    def apply_theme(self) -> None:
        """Restyle every panel and canvas for the current appearance."""
        self.theme = theme_for(self.appearance)
        css = stylesheet(self.theme)
        for widget in (self.stage, self.tabbar, self.fragment_popover, self.symmetry_popover,
                       self.palette_popover):
            widget.setStyleSheet(css)
        self.palette_popover.set_theme(self.theme)
        self.stage.set_theme(self.theme)
        self.fragment_popover.set_theme(self.theme)
        for canvas in self.canvases:
            canvas.set_theme(self.theme)
        self.dialogs.set_theme()
        self.themeChanged.emit(self.theme)

    # -- the Python console -------------------------------------------------------------------

    def show_console(self, shown: bool, remember: bool = True) -> None:
        """Open or close the Python console; ``remember`` keeps the choice for the next start."""
        self.menu_actions.toggles["python_console"].setChecked(shown)
        self.stage.console.set_open(shown, remember)

    def run_script_dialog(self) -> None:
        """Ask for a Python file and run it in the console."""
        path, _ = QFileDialog.getOpenFileName(self, "Run Script", self.files.last_dir(),
                                              "Python scripts (*.py)")
        if path:
            self.files.remember_dir(path)
            self.run_script(path)

    def run_script(self, path: str) -> bool:
        """Run a Python file in the console's namespace (opening the console); return whether it succeeded."""
        self.show_console(True, remember=False)   # the next start opens the console only if the user chose it
        try:
            return self.stage.console.run_file(path)
        except OSError as exc:
            self.stage.console.write(f"Could not read {path}: {exc}\n", "error")
            return False

    def run_code(self, source: str) -> bool:
        """Run Python code in the console's namespace (opening the console); return whether it succeeded."""
        self.show_console(True, remember=False)
        return self.stage.console.run_source(source, "<command line>")

    # -- panels and closing -------------------------------------------------------------------

    def _restore_panels(self) -> None:
        """Put the floating panels where they were left last time, folded or not."""
        for name, panel in self.stage.floating_panels().items():
            self._restore_panel(name, panel)

    def _restore_panel(self, name: str, panel: FloatingPanel) -> None:
        position = str(self.settings.value(f"panels/{name}/position", ""))
        x, _, y = position.partition(",")
        if x.strip().lstrip("-").isdigit() and y.strip().lstrip("-").isdigit():
            panel.placed = QPoint(int(x), int(y))
        if str(self.settings.value(f"panels/{name}/collapsed", "false")).lower() in ("true", "1"):
            panel.set_collapsed(True)

    def track_panel(self, name: str, panel: FloatingPanel) -> None:
        """Restore a panel added later (by a plugin) and remember its arrangement from now on."""
        self._restore_panel(name, panel)
        panel.arranged.connect(lambda: self._save_panel(name, panel))

    def _save_panel(self, name: str, panel: FloatingPanel) -> None:
        """Remember where a panel is and whether it is folded."""
        if panel.placed is None:
            self.settings.remove(f"panels/{name}/position")
        else:
            self.settings.setValue(f"panels/{name}/position", f"{panel.placed.x()},{panel.placed.y()}")
        self.settings.setValue(f"panels/{name}/collapsed", "true" if panel.collapsed else "false")

    @override
    def closeEvent(self, event: QCloseEvent) -> None:
        for canvas in list(self.canvases):
            if not self.files.confirm_discard(canvas.doc):
                event.ignore()
                return
            canvas.doc.modified = False
        self.autosave.close()
        event.accept()
