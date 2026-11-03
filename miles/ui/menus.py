"""The menu bar, declared as data.

Each entry names what it calls, as a dotted path (see :func:`resolve`):
a method of the main window or of one of its parts (``"files.save"``), or
of the canvas in front (``"canvas.commands.paste"``; ``"structure."``
only when that canvas has atoms). :func:`build_menus` turns the
declaration into Qt menus and returns the actions that need updating later
(undo state, checked choices).
"""

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import QMenu, QMenuBar, QWidget

from miles.core.molecule import Representation
from miles.render.cartoon import CartoonColor
from miles.render.scene import REPRESENTATION_NAMES, LabelMode
from miles.render.theme import Appearance
from miles.ui.builder import MouseMode, Tool

Key = QKeySequence.StandardKey
Shortcut = str | QKeySequence.StandardKey | None


@dataclass(frozen=True)
class Command:
    """A menu item that calls ``<slot>(*args)`` (see :func:`resolve`).

    ``role`` lets macOS move the item to the application menu, as it does
    with Settings.
    """

    text: str
    slot: str
    shortcut: Shortcut = None
    args: tuple[Any, ...] = ()
    name: str = ""   # keep the action under this name
    role: QAction.MenuRole | None = None


@dataclass(frozen=True)
class ToggleItem:
    """A checkable item that calls ``<slot>(checked)``, or ``<slot>(**{option: checked})``."""

    text: str
    slot: str
    name: str
    shortcut: Shortcut = None
    checked: bool = True
    option: str = ""


@dataclass(frozen=True)
class ChoiceItem:
    """Mutually exclusive items in a submenu that call ``<slot>(value)``, or ``<slot>(**{option: value})``."""

    title: str | None   # submenu title, or None to list the items in place
    name: str
    slot: str
    options: Sequence[tuple[Any, str, Shortcut]]
    option: str = ""


@dataclass(frozen=True)
class Separator:
    """A line between items; a named one marks where plugins put a group of commands."""

    name: str = ""


SEPARATOR = Separator()


@dataclass(frozen=True)
class Submenu:
    """A nested menu."""

    title: str
    items: list["Item"]


@dataclass(frozen=True)
class Dynamic:
    """A submenu whose items change, filled by ``<slot>(menu)`` each time it opens."""

    title: str
    slot: str
    name: str


Item = Command | ToggleItem | ChoiceItem | Separator | Submenu | Dynamic


@dataclass(frozen=True)
class Menu:
    """A top-level menu."""

    title: str
    items: list[Item]


@dataclass
class MenuActions:
    """Actions kept for later updates."""

    commands: dict[str, QAction] = field(default_factory=dict)
    toggles: dict[str, QAction] = field(default_factory=dict)
    choices: dict[str, dict[Any, QAction]] = field(default_factory=dict)
    menus: dict[str, QMenu] = field(default_factory=dict)

    def check(self, choice: str, value: Any) -> None:
        """Check the item of a choice that has ``value``."""
        action = self.choices[choice].get(value)
        if action is not None:
            action.setChecked(True)


_DRAW_MODES = [rep for rep in REPRESENTATION_NAMES if rep != Representation.HIDDEN]

MENUS = [
    Menu("File", [
        Command("New", "new_document", Key.New),
        Command("Open...", "files.open_dialog", Key.Open),
        Dynamic("Open Recent", "files.fill_recent_menu", "recent"),
        SEPARATOR,
        Command("Save", "files.save", Key.Save),
        Command("Save As...", "files.save_as", Key.SaveAs),
        Command("Export Image...", "files.export_image", "Ctrl+Shift+E"),
        Separator("scripts"),
        Command("Run Script...", "run_script_dialog"),
        SEPARATOR,
        Command("Close", "close_current_tab", Key.Close),
    ]),
    Menu("Edit", [
        Command("Undo", "canvas.commands.undo", Key.Undo, name="undo"),
        Command("Redo", "canvas.commands.redo", Key.Redo, name="redo"),
        SEPARATOR,
        Command("Cut", "canvas.commands.copy", Key.Cut, args=(True,)),
        Command("Copy", "canvas.commands.copy", Key.Copy),
        Command("Paste", "canvas.commands.paste", Key.Paste),
        Command("Duplicate", "canvas.commands.duplicate", "Ctrl+D"),
        Command("Copy SMILES", "canvas.commands.copy_smiles", "Ctrl+Shift+C"),
        SEPARATOR,
        Command("Select All", "canvas.commands.select_all", Key.SelectAll),
        SEPARATOR,
        Command("Clean Geometry", "canvas.commands.clean_geometry", "Ctrl+L"),
        Command("Point Group...", "dialogs.show_symmetry", "Ctrl+G"),
        Command("Add Hydrogens", "canvas.commands.add_hydrogens", "Ctrl+Alt+H"),
        Command("Rebond", "canvas.commands.rebond"),
        Command("Move to Origin", "canvas.commands.move_to_origin"),
        Command("Clear Measurements", "structure.commands.clear_measurements"),
        SEPARATOR,
        Command("Atom List...", "dialogs.show_atom_list", "Ctrl+Shift+A"),
        SEPARATOR,
        Command("Settings...", "dialogs.show_settings", Key.Preferences,
                role=QAction.MenuRole.PreferencesRole),
    ]),
    Menu("Tools", [
        Command("Command Palette...", "dialogs.show_command_palette", "Ctrl+Shift+P"),
        Command("Fragments and SMILES...", "dialogs.search_fragments", "Ctrl+K"),
        Separator("protein tools"),
        Separator("tool choice"),
        ChoiceItem(None, "tool", "builder.set_tool",
               [(Tool.SELECT, "Select", "Ctrl+1"), (Tool.BUILD, "Build", "Ctrl+2"),
                (Tool.DELETE, "Delete", "Ctrl+3")]),
    ]),
    Menu("View", [
        ChoiceItem("Style", "mode", "canvas.commands.set_display_mode",
               [(rep, REPRESENTATION_NAMES[rep], f"Ctrl+Alt+{n}")
                for n, rep in enumerate(_DRAW_MODES, start=1)]),
        ChoiceItem("Labels", "labels", "canvas.set_style",
               [(LabelMode.NONE, "None", None), (LabelMode.SYMBOL, "Element", None),
                (LabelMode.NUMBER, "Number", None), (LabelMode.BOTH, "Element and Number", None)],
               option="labels"),
        ChoiceItem("Cartoon Colors", "cartoon_color", "canvas.set_style",
               [(CartoonColor.SECONDARY_STRUCTURE, "Secondary Structure", None),
                (CartoonColor.CHAIN, "Chain", None), (CartoonColor.RAINBOW, "Rainbow", None)],
               option="cartoon_color"),
        ChoiceItem("Mouse Mode", "mouse_mode", "builder.set_mouse_mode",
               [(MouseMode.ROTATE, "Rotate (R)", None), (MouseMode.TRANSLATE, "Translate (T)", None),
                (MouseMode.SCALE, "Scale (S)", None)]),
        ToggleItem("Hydrogens", "canvas.set_style", "show_hydrogens", "Ctrl+Shift+H",
                   option="show_hydrogens"),
        ToggleItem("Stereo Labels", "canvas.set_style", "stereo", checked=False, option="stereo"),
        ToggleItem("Depth Fade", "canvas.set_style", "depth_fade", option="depth_fade"),
        SEPARATOR,
        Command("Show Only Selection", "structure.commands.show_only_selection", "Ctrl+Shift+F"),
        Command("Hide Selection", "structure.commands.hide_selection"),
        Command("Show Everything", "structure.commands.show_everything", "Ctrl+Shift+U"),
        ToggleItem("Python Console", "show_console", "python_console", "Ctrl+Alt+P", checked=False),
        ChoiceItem("Appearance", "appearance", "set_appearance",
               [(Appearance.SYSTEM, "Match System", None), (Appearance.LIGHT, "Light", None),
                (Appearance.DARK, "Dark", None)]),
        SEPARATOR,
        Command("Fit to Window", "canvas.view.fit", "Ctrl+0", args=(True,)),
        Command("Center", "canvas.view.center", "Ctrl+E"),
        Command("Reset Orientation", "canvas.view.reset", "Ctrl+Shift+R"),
        *(Command(f"Look Along {axis.upper()}", "canvas.view.look_along", args=(axis,)) for axis in "xyz"),
    ]),
    Menu("Help", [
        Command("Shortcuts", "dialogs.show_shortcuts", "Ctrl+/"),
        Command("Plugins...", "dialogs.show_plugins"),
    ]),
]


def _action(window: QWidget, menu: QMenu, text: str, shortcut: Shortcut, checkable: bool = False) -> QAction:
    action = QAction(text, window)
    if shortcut is not None:
        action.setShortcut(QKeySequence(shortcut))
    action.setCheckable(checkable)
    menu.addAction(action)
    return action


def resolve(window: QWidget, slot: str) -> Callable[..., Any]:
    """Return what a menu item calls.

    Parameters
    ----------
    window : QWidget
        The main window.
    slot : str
        A dotted path: from the window (``"new_document"``,
        ``"files.save"``), or from the canvas in front when it starts with
        ``"canvas."`` (``"canvas.commands.paste"``), looked up each time
        the item is chosen. ``"structure."`` is the canvas in front only
        when its structure has atoms; otherwise the item does nothing.

    Returns
    -------
    callable
        The method to call.
    """
    head, _, rest = slot.partition(".")
    if head not in ("canvas", "structure"):
        return _attribute(window, slot)

    def on_current_canvas(*args: Any, **kwargs: Any) -> Any:
        canvas = window.current_canvas()   # type: ignore[attr-defined]
        if canvas is None or (head == "structure" and not len(canvas.mol)):
            return None
        return _attribute(canvas, rest)(*args, **kwargs)

    return on_current_canvas


def _attribute(owner: object, path: str) -> Any:
    target = owner
    for name in path.split("."):
        target = getattr(target, name)
    return target


def build_menus(window: QWidget, bar: QMenuBar, menus: list[Menu] = MENUS) -> MenuActions:
    """Create the menus on a menu bar, connected to what their items call (:func:`resolve`).

    Parameters
    ----------
    window : QWidget
        The main window; it owns the actions.
    bar : QMenuBar
        Menu bar to fill.
    menus : list of Menu, optional
        Declaration.

    Returns
    -------
    MenuActions
        Named commands, toggles and choice groups.
    """
    actions = MenuActions()
    for spec in menus:
        _fill(window, bar.addMenu(spec.title), spec.items, actions)
    return actions


def _fill(window: QWidget, menu: QMenu, items: list[Item], actions: MenuActions) -> None:
    """Add items to a menu, recursing into submenus."""
    for item in items:
        match item:
            case Separator(name):
                menu.addSeparator().setObjectName(name)
            case Submenu(title, subitems):
                _fill(window, menu.addMenu(title), subitems, actions)
            case Dynamic(title, slot, name):
                submenu = menu.addMenu(title)
                fill = resolve(window, slot)
                fill(submenu)   # filled now as well, so the command palette sees the items
                submenu.aboutToShow.connect(lambda m=submenu, f=fill: f(m))
                actions.menus[name] = submenu
            case Command(text, slot, shortcut, args, name, role):
                action = _action(window, menu, text, shortcut)
                if role is not None:
                    action.setMenuRole(role)
                method = resolve(window, slot)
                action.triggered.connect(lambda _=False, f=method, a=args: f(*a))
                if name:
                    actions.commands[name] = action
            case ToggleItem(text, slot, name, shortcut, checked, option):
                action = _action(window, menu, text, shortcut, checkable=True)
                action.setChecked(checked)
                action.triggered.connect(_with_value(resolve(window, slot), option))
                actions.toggles[name] = action
            case ChoiceItem(title, name, slot, options, option):
                target = menu if title is None else menu.addMenu(title)
                group = QActionGroup(window)
                actions.choices[name] = {}
                choose = _with_value(resolve(window, slot), option)
                for value, text, shortcut in options:
                    action = _action(window, target, text, shortcut, checkable=True)
                    group.addAction(action)
                    action.triggered.connect(lambda _=False, f=choose, v=value: f(v))
                    actions.choices[name][value] = action
                actions.check(name, options[0][0])


def _with_value(method: Callable[..., Any], option: str) -> Callable[[Any], Any]:
    """Call ``method(value)``, or ``method(**{option: value})`` when an option is named."""
    if not option:
        return method
    return lambda value: method(**{option: value})
