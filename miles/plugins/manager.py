"""Loading plugins and connecting what they register to the main window (:class:`PluginManager`)."""

import functools
import importlib
import os
import sys
from collections.abc import Callable, Iterable
from dataclasses import replace
from types import ModuleType
from typing import TYPE_CHECKING, Any

from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import QMenu, QWidget

from miles.api import App, DocumentView, document_view
from miles.core.document import Change
from miles.core.molecule import Molecule
from miles.files.common import FileFormatError, JobInfo
from miles.files.formats import REGISTRY
from miles.plugins.colors import named_colors
from miles.plugins.context import MenuAnchor, PluginContext
from miles.plugins.discovery import (
    PLUGIN_API,
    Candidate,
    PluginInfo,
    PluginState,
    discover,
    forget_modules,
    plugin_traceback,
)
from miles.plugins.layers import PluginLayer
from miles.plugins.settings import Choice, Setting, Toggle
from miles.render.layers import Annotations, DrawContext, SceneItems
from miles.ui.canvas.extensions import SECTIONS, Layer, MenuItem, MenuTarget
from miles.ui.lookup import StructureSource

if TYPE_CHECKING:
    from miles.ui.canvas import Canvas
    from miles.ui.floating import FloatingPanel
    from miles.ui.mainwindow import MainWindow

DISABLED_KEY = "plugins/disabled"
"""Settings key listing the plugins switched off, separated by ``|``."""

DEFAULT_MENU = "Plugins"
"""Top-level menu for commands whose path names no menu."""


def _remove_action(menu: QWidget, action: QAction) -> None:
    menu.removeAction(action)
    action.setShortcut(QKeySequence())   # free the shortcut now: the action itself goes only later
    action.deleteLater()


def _remove_if_empty(parent: QWidget, menu: QMenu) -> None:
    """Remove a menu a plugin created once nothing is left in it."""
    if not [a for a in menu.actions() if not a.isSeparator()]:
        parent.removeAction(menu.menuAction())
        menu.deleteLater()


def _extensions(extensions: str | Iterable[str]) -> list[str]:
    items = [extensions] if isinstance(extensions, str) else list(extensions)
    return [e.lower().lstrip(".") for e in items if e.strip(". ")]


class PluginManager:
    """Loads plugins and connects what they register to the main window.

    Parameters
    ----------
    window : MainWindow
        The window to extend; it must already have its ``api``.
    """

    def __init__(self, window: "MainWindow") -> None:
        self.window = window
        self.app: App = window.api
        self.plugins: list[PluginInfo] = []
        self.toggles: dict[str, Toggle] = {}
        self.choices: dict[str, Choice] = {}
        self.settings: dict[str, list[Setting]] = {}   # by plugin key
        self._checks: dict[QAction, Callable[[], bool]] = {}   # whether a plugin command can run now
        self._layers: list[PluginLayer] = []   # drawings, worked out again when the selection changes
        self._candidates: dict[str, Candidate] = {}
        self._cleanups: dict[str, list[Callable[[], None]]] = {}   # undo what each plugin added
        self._front: Canvas | None = None   # tab last reported to on_current_changed
        self._listeners: dict[str, list[tuple[PluginInfo, Callable[..., Any]]]] = {
            event: [] for event in ("opened", "closed", "current", "selection", "structure", "frame",
                                    "changed", "theme")}
        window.documentOpened.connect(lambda canvas: self._notify("opened", canvas))
        window.documentClosed.connect(lambda canvas: self._notify("closed", canvas))
        window.currentDocumentChanged.connect(lambda canvas: self._notify("current", canvas))
        window.documentSelectionChanged.connect(self._selection_changed)
        window.documentStructureChanged.connect(lambda canvas: self._notify("structure", canvas))
        window.documentChanged.connect(self._changed)
        window.themeChanged.connect(lambda theme: self._call("theme", theme))
        window.currentDocumentChanged.connect(self._sync_toggles)
        window.documentOpened.connect(self._settle_toggles)

    def _selection_changed(self, canvas: "Canvas") -> None:
        """Draw the plugins' layers again, since they may show the selection, then tell the plugins."""
        drawn = [handle for handle in self._layers if handle.kind != "status"]
        for handle in drawn:
            handle.layer.invalidate()
        if drawn and canvas in self.window.canvases:
            canvas.redraw()
        self._notify("selection", canvas)

    # -- settings ---------------------------------------------------------------------------

    @property
    def disabled(self) -> set[str]:
        """Keys of the plugins switched off in Help > Plugins."""
        value = str(self.window.settings.value(DISABLED_KEY, "") or "")
        return {key for key in value.split("|") if key}

    def set_enabled(self, key: str, enabled: bool, now: bool = True) -> None:
        """Switch a plugin on or off, at once and at later starts.

        Parameters
        ----------
        key : str
            :attr:`PluginInfo.key` of the plugin.
        enabled : bool
            Switch it on (True) or off.
        now : bool, optional
            Also load or unload it now (False: only from the next start).
        """
        keys = self.disabled
        if enabled:
            keys.discard(key)
        else:
            keys.add(key)
        self.window.settings.setValue(DISABLED_KEY, "|".join(sorted(keys)))
        info = next((i for i in self.plugins if i.key == key), None)
        if not now or info is None or info.state == PluginState.SKIPPED:
            return
        if enabled and info.state in (PluginState.OFF, PluginState.FAILED):
            self.reload(info)
        elif not enabled and info.state in (PluginState.LOADED, PluginState.FAILED):
            self.unload(info)
            info.state = PluginState.OFF

    def unload(self, info: PluginInfo) -> None:
        """Take out of the window everything a plugin added: commands, panels, drawings, sources ..."""
        for cleanup in reversed(self._cleanups.pop(info.key, [])):
            try:
                cleanup()
            except Exception:   # undoing must go on even if one step fails
                info.error += plugin_traceback()
        for document in self.app.documents:
            document.forget_cache(info.key)
        info.commands.clear()
        info.extras.clear()
        for canvas in self.window.canvases:
            canvas.update_status()
            canvas.redraw()

    def reload(self, info: PluginInfo) -> None:
        """Unload a plugin, import it again (so edits to its file count) and register it.

        A plugin from the plugin folder is read from its file again; for the
        others the module is reloaded with :func:`importlib.reload`.
        """
        candidate = self._candidates.get(info.key)
        if candidate is None:
            return
        self.unload(info)
        info.error = ""
        info.warnings = ""
        if not info.path:
            if info.builtin:
                module_name: str | None = info.key
            elif info.key.startswith("module:"):
                module_name = info.key.removeprefix("module:")
            else:
                value = getattr(candidate.load, "__self__", None)
                module_name = str(getattr(value, "value", "")).split(":")[0] or None
            loaded = sys.modules.get(module_name) if module_name else None
            if module_name and loaded is not None and getattr(loaded, "__spec__", None) is not None:
                forget_modules(module_name)   # the package and its own modules are imported again
                if info.key.startswith("module:"):
                    try:
                        module = importlib.import_module(module_name)
                    except Exception:
                        info.state = PluginState.FAILED
                        info.error = plugin_traceback()
                        return
                    candidate = Candidate(info, lambda: module)
                    self._candidates[info.key] = candidate
        self._register(candidate)

    def _on_unload(self, info: PluginInfo, cleanup: Callable[[], Any]) -> None:
        self._cleanups.setdefault(info.key, []).append(cleanup)

    def on_unload(self, info: PluginInfo, callback: Callable[[], Any]) -> None:
        """Call a plugin's own clean-up when it is switched off or reloaded.

        See :meth:`PluginContext.on_unload`.
        """
        def run() -> None:
            try:
                callback()
            except Exception:
                info.error += plugin_traceback()
        self._on_unload(info, run)

    # -- loading ----------------------------------------------------------------------------

    def load(self, third_party: bool = True) -> list[PluginInfo]:
        """Find and register the plugins.

        Parameters
        ----------
        third_party : bool, optional
            Load plugins of other packages and of the plugin folder too
            (False for ``--no-plugins``; they are then listed as skipped).

        Returns
        -------
        list of PluginInfo
            Every plugin found, with its state.
        """
        disabled = self.disabled
        candidates = discover(third_party=True)
        for candidate in candidates:
            info = candidate.info
            self.plugins.append(info)
            self._candidates[info.key] = candidate
            if not info.builtin and not third_party:
                info.state = PluginState.SKIPPED
            elif info.key in disabled:
                info.state = PluginState.OFF
            else:
                self._register(candidate)
        return self.plugins

    def add_module(self, module: ModuleType) -> PluginInfo:
        """Register a module as a plugin now, e.g. a plugin under test (see :mod:`miles.testing`).

        Parameters
        ----------
        module : module
            The plugin's module, with its ``register(context)`` function.

        Returns
        -------
        PluginInfo
            The plugin, loaded or failed with its error.
        """
        info = PluginInfo(f"module:{module.__name__}", module.__name__.rsplit(".", 1)[-1], "added")
        candidate = Candidate(info, lambda: module)
        self.plugins.append(info)
        self._candidates[info.key] = candidate
        self._register(candidate)
        return info

    def _register(self, candidate: Candidate) -> None:
        info = candidate.info
        try:
            loaded = candidate.load()
            if isinstance(loaded, ModuleType):   # a module with register(context)
                module: ModuleType | None = loaded
                register = getattr(loaded, "register", None)
            else:   # an entry point naming the register function itself
                module = sys.modules.get(getattr(loaded, "__module__", "") or "")
                register = loaded
            if module is not None:
                info.name = str(getattr(module, "NAME", "") or info.name)
                info.version = str(getattr(module, "VERSION", "") or "")
                summary = (module.__doc__ or "").strip().split("\n")[0]
                info.description = str(getattr(module, "DESCRIPTION", "") or summary)
            if not callable(register):
                raise TypeError("The plugin has no register(context) function")
            required = int(getattr(module, "REQUIRES_API", 1) or 1) if module is not None else 1
            if required > PLUGIN_API:
                raise RuntimeError(f"The plugin needs version {required} of the plugin interface; "
                                   f"this MILES has version {PLUGIN_API}. Update MILES to use it.")
            register(PluginContext(self, info))
        except Exception:
            info.state = PluginState.FAILED
            info.error = plugin_traceback()
            print(f"Plugin {info.name} ({info.path or info.source}) failed to load:\n{info.error}",
                  file=sys.stderr)
            return
        info.state = PluginState.LOADED

    # -- guarded calls ----------------------------------------------------------------------

    def failed(self, info: PluginInfo) -> None:
        """Record the exception being handled and report it through the error dialog."""
        self._failed(info)

    def _failed(self, info: PluginInfo) -> None:
        """Record the exception being handled and report it through the error dialog."""
        info.error = plugin_traceback()
        sys.excepthook(*sys.exc_info())   # the error dialog shows each distinct error once

    def guard_io(self, info: PluginInfo, func: Callable[..., Any]) -> Callable[..., Any]:
        """Wrap a plugin's reader or writer so unexpected errors become a message, not a crash."""
        def guarded(*args: Any) -> Any:
            try:
                return func(*args)
            except (FileFormatError, OSError, ValueError):
                raise
            except Exception as exc:
                info.error = plugin_traceback()
                raise FileFormatError(f"The {info.name} plugin could not handle this file: {exc}") from exc
        return guarded

    # -- commands ---------------------------------------------------------------------------

    def menu_for(self, titles: list[str], info: PluginInfo | None = None) -> QMenu:
        """Find or create the menu at a path of menu titles (the first one top-level).

        Menus created for a plugin (``info``) are removed with it when they are empty.
        """
        bar = self.window.menuBar()
        top: QMenu | None = None
        for action in bar.actions():
            menu = action.menu()
            if isinstance(menu, QMenu) and action.text().replace("&", "") == titles[0]:
                top = menu
                break
        if top is None:
            top = QMenu(titles[0], self.window)
            help_action = next((a for a in bar.actions() if a.text().replace("&", "") == "Help"), None)
            if help_action is not None:
                bar.insertMenu(help_action, top)
            else:
                bar.addMenu(top)
            if info is not None:
                self._on_unload(info, functools.partial(_remove_if_empty, bar, top))
        menu = top
        for title in titles[1:]:
            sub = next((a.menu() for a in menu.actions()
                        if isinstance(a.menu(), QMenu) and a.text().replace("&", "") == title), None)
            if isinstance(sub, QMenu):
                menu = sub
                continue
            created = menu.addMenu(title)
            if info is not None:
                self._on_unload(info, functools.partial(_remove_if_empty, menu, created))
            menu = created
        return menu

    def _shortcut_taken(self, sequence: QKeySequence) -> bool:
        return any(a.shortcut() == sequence for a in self.window.findChildren(QAction)
                   if not a.shortcut().isEmpty())

    @staticmethod
    def _parts(path: str) -> list[str]:
        parts = [p.strip() for p in path.split(">") if p.strip()]
        if not parts:
            raise ValueError("A menu item needs a title")
        return [DEFAULT_MENU, *parts] if len(parts) == 1 else parts

    @staticmethod
    def _place(menu: QMenu, action: QAction, after: MenuAnchor, before: MenuAnchor) -> None:
        """Insert an action after the first ``after`` anchor found, else before the first ``before``."""
        actions = menu.actions()

        def find(anchor: MenuAnchor) -> int | None:
            names = [anchor] if isinstance(anchor, str) else list(anchor or [])
            for name in names:
                for k, a in enumerate(actions):
                    if a.text().replace("&", "") == name or (a.isSeparator() and a.objectName() == name):
                        return k
            return None

        k = find(after)
        if k is not None:
            if k + 1 < len(actions):
                menu.insertAction(actions[k + 1], action)
            else:
                menu.addAction(action)
            return
        k = find(before)
        if k is not None:
            menu.insertAction(actions[k], action)
        else:
            menu.addAction(action)

    def add_menu(self, path: str, after: MenuAnchor, before: MenuAnchor,
                 info: PluginInfo | None = None) -> QMenu:
        """Create a submenu at a place (see :meth:`PluginContext.add_menu`)."""
        parts = [p.strip() for p in path.split(">") if p.strip()]
        if len(parts) < 2:
            return self.menu_for(parts or [DEFAULT_MENU], info)
        parent = self.menu_for(parts[:-1], info)
        existing = next((a.menu() for a in parent.actions()
                         if isinstance(a.menu(), QMenu) and a.text().replace("&", "") == parts[-1]), None)
        if isinstance(existing, QMenu):
            return existing
        menu = QMenu(parts[-1], parent)
        self._place(parent, menu.menuAction(), after, before)
        if info is not None:
            self._on_unload(info, lambda: _remove_if_empty(parent, menu))
        return menu

    def add_separator(self, path: str, after: MenuAnchor, before: MenuAnchor,
                      info: PluginInfo | None = None) -> None:
        """Add a separator to a menu (see :meth:`PluginContext.add_separator`)."""
        menu = self.menu_for([p.strip() for p in path.split(">") if p.strip()] or [DEFAULT_MENU], info)
        separator = QAction(menu)
        separator.setSeparator(True)
        self._place(menu, separator, after, before)
        if info is not None:
            self._on_unload(info, lambda: _remove_action(menu, separator))

    def add_command(self, info: PluginInfo, path: str, callback: Callable[[App], Any], shortcut: str | None,
                    enabled: Callable[[App], bool] | None, after: MenuAnchor,
                    before: MenuAnchor = None) -> QAction:
        """Add a plugin's menu command (see :meth:`PluginContext.add_command`)."""
        parts = self._parts(path)
        menu = self.menu_for(parts[:-1], info)
        title = parts[-1]
        action = QAction(title, self.window)
        if shortcut:
            sequence = QKeySequence(shortcut)
            if self._shortcut_taken(sequence):
                note = f"Shortcut {shortcut} of {title} is already used; the command has none.\n"
                info.warnings += note
                print(f"Plugin {info.name}: {note}", file=sys.stderr, end="")
            else:
                action.setShortcut(sequence)

        def allowed() -> bool:
            if enabled is None:
                return True
            try:
                return bool(enabled(self.app))
            except Exception:
                self._failed(info)
                return False

        def run() -> None:
            if not allowed():   # chosen from the menu or a shortcut while it could not run
                self.app.message(f"{title} is not available now")
                return
            try:
                callback(self.app)
            except Exception:
                self._failed(info)

        action.triggered.connect(lambda _=False: run())
        if enabled is not None:
            self._checks[action] = allowed
        self._place(menu, action, after, before)
        alive = [True]   # the menu's signals stay connected after the plugin is unloaded
        if enabled is not None:
            menu.aboutToShow.connect(lambda: action.setEnabled(allowed()) if alive[0] else None)
            menu.aboutToHide.connect(lambda: action.setEnabled(True) if alive[0] else None)   # run checks
        info.commands.append(" > ".join(parts))

        def remove() -> None:
            alive[0] = False
            self._checks.pop(action, None)
            _remove_action(menu, action)

        self._on_unload(info, remove)
        return action

    def available(self, action: QAction) -> bool:
        """Whether a menu command can run now: its own check for a plugin's command, else its state."""
        check = self._checks.get(action)
        return action.isEnabled() and (check is None or check())

    # -- panels -----------------------------------------------------------------------------

    def add_panel(self, info: PluginInfo, title: str, factory: Callable[[App], QWidget], menu: str,
                  shortcut: str | None) -> "FloatingPanel | None":
        """Add a plugin's panel (see :meth:`PluginContext.add_panel`)."""
        from miles.ui.plugin_panel import PluginPanel

        try:
            content = factory(self.app)
        except Exception:
            self._failed(info)
            return None
        panel = PluginPanel(title)
        panel.body_layout.addWidget(content)
        name = f"plugin-{info.key}-{title}".replace(" ", "_").replace(":", "_").lower()
        self.add_floating_panel(name, panel, "right", info)
        self.add_command(info, f"{menu} > {title}", lambda _app: panel.toggle(), shortcut, None, None)
        return panel

    def add_floating_panel(self, name: str, panel: "FloatingPanel", place: str = "right",
                           info: PluginInfo | None = None) -> None:
        """Put a plugin's own panel in the window (see :meth:`PluginContext.add_floating_panel`)."""
        if name in self.window.stage.floating_panels():
            raise ValueError(f"There is already a panel named {name!r}")
        self.window.stage.add_panel(name, panel, place)
        self.window.track_panel(name, panel)
        if info is not None:
            info.extras.append(f"panel {name}")
            self._on_unload(info, lambda: self.window.stage.remove_panel(name))

    def add_file_functions(self, info: PluginInfo, extensions: str | Iterable[str],
                           reader: Callable[[str], tuple[Molecule, JobInfo]] | None = None,
                           writer: Callable[[str, Molecule, JobInfo], None] | None = None) -> None:
        """Register a plugin's file reader or writer (see :meth:`PluginContext.add_file_reader`)."""
        guarded_reader = self.guard_io(info, reader) if reader is not None else None
        guarded_writer = self.guard_io(info, writer) if writer is not None else None
        for ext in _extensions(extensions):
            REGISTRY.register(ext, reader=guarded_reader, writer=guarded_writer)
            self._on_unload(info, functools.partial(REGISTRY.unregister, ext, guarded_reader, guarded_writer))

    def add_fragments(self, info: PluginInfo, path: str) -> None:
        """Add a plugin's fragment library (see :meth:`PluginContext.add_fragments`)."""
        added = self.window.add_fragment_library(path)
        info.extras.append(f"fragments from {os.path.basename(path)}")
        self._on_unload(info, lambda: self.window.remove_fragments(added))

    # -- per-tab state and toggles ----------------------------------------------------------

    def state(self, info: PluginInfo, document: DocumentView) -> dict[str, Any]:
        """Return the dictionary a plugin keeps for a tab (saved in sessions under ``plugins/<key>``)."""
        return document.canvas.plugin_state.setdefault(info.key, {})

    def add_toggle(self, info: PluginInfo, path: str, name: str,
                   changed: Callable[[DocumentView, bool], Any] | None,
                   default: bool | Callable[[DocumentView], bool], shortcut: str | None,
                   after: MenuAnchor, before: MenuAnchor) -> Toggle:
        """Add a per-tab toggle (see :meth:`PluginContext.add_toggle`)."""
        if name in self.toggles or name in self.choices:
            raise ValueError(f"There is already a toggle or choice named {name!r}")
        toggle: Toggle | None = None

        def switched(app: App) -> None:
            doc = app.current
            if toggle is None or doc is None:
                action.setChecked(False)
                return
            toggle.set(doc, action.isChecked())

        action = self.add_command(info, path, switched, shortcut, None, after, before)
        action.setCheckable(True)
        toggle = self.toggles[name] = Toggle(self, info, name, action, default, changed)
        self._on_unload(info, lambda: self.toggles.pop(name, None))
        for canvas in self.window.canvases:
            self._settle(toggle, document_view(self.window, canvas))
        current = self.app.current
        if current is not None:
            action.setChecked(toggle.value(current))
        return toggle

    def add_choice(self, info: PluginInfo, path: str, name: str, options: list[tuple[str, str]],
                   changed: Callable[[DocumentView, str], Any] | None,
                   default: str | Callable[[DocumentView], str] | None,
                   enabled: Callable[[DocumentView, str], bool] | None,
                   read: Callable[[DocumentView], str] | None,
                   write: Callable[[DocumentView, str], Any] | None,
                   after: MenuAnchor, before: MenuAnchor) -> Choice:
        """Add a per-tab choice between menu items (see :meth:`PluginContext.add_choice`)."""
        if name in self.choices or name in self.toggles:
            raise ValueError(f"There is already a toggle or choice named {name!r}")
        if not options:
            raise ValueError("A choice needs at least one option")
        choice = Choice(self, info, name, options, default if default is not None else options[0][0], changed,
                        enabled, read, write)
        group = QActionGroup(self.window)
        group.setExclusive(True)
        previous: MenuAnchor = after
        for k, (value, title) in enumerate(options):
            def chosen(app: App, value: str = value) -> None:
                if app.current is not None:
                    choice.set(app.current, value)

            action = self.add_command(info, f"{path} > {title}", chosen, None, None, previous,
                                      before if k == 0 else None)
            action.setCheckable(True)
            group.addAction(action)
            choice.actions[value] = action
            previous = title
        menu = self.menu_for(self._parts(f"{path} > x")[:-1], info)
        alive = [True]
        menu.aboutToShow.connect(lambda: choice.sync(self.app.current) if alive[0] else None)
        self.choices[name] = choice

        def remove() -> None:
            alive[0] = False
            self.choices.pop(name, None)
            group.setParent(None)   # out of the window now, deleted once nothing holds it
            group.deleteLater()

        self._on_unload(info, remove)
        if read is None:
            for canvas in self.window.canvases:
                self._settle_choice(choice, document_view(self.window, canvas))
        choice.sync(self.app.current)
        return choice

    def _settle_choice(self, choice: Choice, document: DocumentView) -> None:
        """Store a choice's default for a tab, so that sessions save it even if it was never chosen."""
        state = self.state(choice.info, document)
        if choice.name not in state:
            try:
                state[choice.name] = choice.value(document)
            except Exception:
                self._failed(choice.info)

    def add_setting(self, info: PluginInfo, name: str, title: str, default: bool | int | float | str,
                    minimum: float | None, maximum: float | None, step: float | None, unit: str,
                    changed: Callable[[Any], Any] | None) -> Setting:
        """Add a plugin's setting (see :meth:`PluginContext.add_setting`)."""
        settings = self.settings.setdefault(info.key, [])
        if any(s.name == name for s in settings):
            raise ValueError(f"The plugin already has a setting named {name!r}")
        setting = Setting(self, info, name, title, default, minimum, maximum, step, unit, changed)
        settings.append(setting)
        self._on_unload(info, lambda: self.settings.pop(info.key, None))
        return setting

    def setting_value(self, plugin: str, name: str, default: Any) -> Any:
        """Return the value of a plugin's setting (see :meth:`PluginContext.setting_value`)."""
        found = next((s for s in self.settings.get(plugin, []) if s.name == name), None)
        if found is not None:
            return found.value
        stored = Setting(self, PluginInfo(plugin, plugin, ""), name, name, default,
                         None, None, None, "", None)
        return stored.value

    def _settle(self, toggle: Toggle, document: DocumentView) -> None:
        """Store a toggle's default for a tab, so that sessions save it even if it was never switched."""
        state = self.state(toggle.info, document)
        if toggle.name not in state:
            try:
                state[toggle.name] = toggle.value(document)
            except Exception:
                self._failed(toggle.info)

    def _settle_toggles(self, canvas: "Canvas") -> None:
        if canvas in self.window.canvases:
            view = document_view(self.window, canvas)
            for toggle in self.toggles.values():
                self._settle(toggle, view)
            for choice in self.choices.values():
                if choice._read is None:
                    self._settle_choice(choice, view)

    def _sync_toggles(self, canvas: "Canvas") -> None:
        if canvas not in self.window.canvases:
            return
        view = document_view(self.window, canvas)
        for toggle in self.toggles.values():
            try:
                toggle.action.setChecked(toggle.value(view))
            except Exception:
                self._failed(toggle.info)
        for choice in self.choices.values():
            choice.sync(view)

    # -- drawing, status, menus, sources ------------------------------------------------------

    def add_layer(self, info: PluginInfo, kind: str, function: Callable[..., Any]) -> PluginLayer:
        """Add annotations, scene items or status text of a plugin to every canvas."""
        expected: type = {"annotations": Annotations, "scene": SceneItems, "status": str}[kind]
        handle = PluginLayer(self, kind, info)

        def call(canvas: Any, context: DrawContext | None) -> Any:
            state = (id(canvas.mol), canvas.doc.revision, handle.layer.generation)
            if handle.failed_at.get(id(canvas)) == state:
                return None
            document = document_view(self.window, canvas)
            try:
                result = function(document) if kind == "status" else function(document, context)
                if result is not None and not isinstance(result, expected):
                    raise TypeError(f"{kind} functions must return {expected.__name__} or None, "
                                    f"not {type(result).__name__}")
                if isinstance(result, Annotations):   # colours may be given in any form
                    result = Annotations([named_colors(x) for x in result.lines],
                                         [named_colors(x) for x in result.labels],
                                         [named_colors(x) for x in result.arcs])
                return result
            except Exception:
                handle.failed_at[id(canvas)] = state
                self._failed(info)
                return None

        handle.layer = Layer(call)
        info.extras.append({"annotations": "lines and labels over the molecule", "scene": "items drawn in 3D",
                            "status": "text in the status line"}[kind])
        extensions = self.window.extensions
        layers = {"annotations": extensions.annotations, "scene": extensions.scene_items,
                  "status": extensions.status}[kind]
        layers.append(handle.layer)
        self._layers.append(handle)

        def remove() -> None:
            layers.remove(handle.layer)
            self._layers.remove(handle)

        self._on_unload(info, remove)
        return handle

    def add_context_items(self, info: PluginInfo,
                          function: Callable[[DocumentView, MenuTarget], Iterable[MenuItem]],
                          section: str) -> None:
        """Add a plugin's right-click menu items (see :meth:`PluginContext.add_context_items`)."""
        if section not in SECTIONS:
            raise ValueError(f"The section must be one of {', '.join(SECTIONS)}")

        def guarded(item: MenuItem) -> MenuItem:
            action = item.action
            if action is not None:
                def run(act: Callable[[], Any] = action) -> None:
                    try:
                        act()
                    except Exception:
                        self._failed(info)
                return replace(item, action=run, items=[guarded(child) for child in item.items])
            return replace(item, items=[guarded(child) for child in item.items])

        def call(canvas: Any, target: MenuTarget) -> list[MenuItem]:
            try:
                items = list(function(document_view(self.window, canvas), target) or [])
            except Exception:
                self._failed(info)
                return []
            return [guarded(item) for item in items]

        functions = self.window.extensions.menu_items[section]
        functions.append(call)
        info.extras.append(f"right-click menu items ({section})")
        self._on_unload(info, lambda: functions.remove(call))

    def add_structure_source(self, info: PluginInfo, source: StructureSource) -> None:
        """Add a plugin's structure source (see :meth:`PluginContext.add_structure_source`)."""
        fetch = source.fetch

        def guarded(text: str) -> tuple[Molecule, JobInfo]:
            try:
                return fetch(text)
            except (FileFormatError, OSError, ValueError):
                raise
            except Exception as exc:
                info.error = plugin_traceback()
                raise FileFormatError(f"The {info.name} plugin could not fetch {text}: {exc}") from exc

        added = replace(source, fetch=guarded)
        self.window.structure_sources.append(added)
        info.extras.append(f"structure source {source.name}")
        self._on_unload(info, lambda: self.window.structure_sources.remove(added))

    # -- events -----------------------------------------------------------------------------

    def listen(self, event: str, info: PluginInfo, callback: Callable[..., Any]) -> None:
        """Call ``callback`` on an event; document events pass a :class:`~miles.api.DocumentView`."""
        entry = (info, callback)
        listeners = self._listeners[event]
        listeners.append(entry)
        self._on_unload(info, lambda: listeners.remove(entry))

    def _call(self, event: str, *args: Any) -> None:
        for info, callback in self._listeners[event]:
            try:
                callback(*args)
            except Exception:
                self._failed(info)

    def _notify(self, event: str, canvas: "Canvas", *extra: Any) -> None:
        if event == "current":
            if canvas is self._front:
                return
            self._front = canvas
        if self._listeners[event] and canvas in self.window.canvases:
            self._call(event, document_view(self.window, canvas), *extra)

    def _changed(self, canvas: "Canvas", change: str) -> None:
        kind = Change(change)
        if kind in (Change.FRAME, Change.FRAMES):
            self._notify("frame", canvas)
        self._notify("changed", canvas, kind)
