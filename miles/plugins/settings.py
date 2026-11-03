"""Toggles and choices kept for each tab, and settings, that plugins add.

See :class:`Toggle`, :class:`Choice` and :class:`Setting`.
"""

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

from PySide6.QtGui import QAction

from miles.api import DocumentView
from miles.plugins.discovery import PluginInfo

if TYPE_CHECKING:
    from miles.plugins.manager import PluginManager


class Toggle:
    """A menu command that switches something on or off for each tab (see :meth:`PluginContext.add_toggle`).

    Attributes
    ----------
    name : str
        Name of the setting.
    action : QAction
        The checkable menu item.
    """

    def __init__(self, manager: "PluginManager", info: PluginInfo, name: str, action: QAction,
                 default: bool | Callable[[DocumentView], bool],
                 changed: Callable[[DocumentView, bool], Any] | None) -> None:
        self._manager = manager
        self.info = info
        self.name = name
        self.action = action
        self._default = default
        self._changed = changed

    def value(self, document: DocumentView) -> bool:
        """Whether the toggle is on for a tab.

        A tab never switched takes what an older session saved among its
        display options, else the default.
        """
        state = self._manager.state(self.info, document)
        if self.name in state:
            return bool(state[self.name])
        older = document.info.view.get("style")
        if isinstance(older, dict) and self.name in older:
            return bool(older[self.name])
        return bool(self._default(document) if callable(self._default) else self._default)

    def set(self, document: DocumentView, on: bool) -> None:
        """Switch the toggle on or off for a tab, as choosing the menu item does."""
        self._manager.state(self.info, document)[self.name] = bool(on)
        if document == self._manager.app.current:
            self.action.setChecked(bool(on))
        if self._changed is not None:
            try:
                self._changed(document, bool(on))
            except Exception:
                self._manager.failed(self.info)


class Choice:
    """Menu items of which one is checked for each tab, like radio buttons.

    See :meth:`PluginContext.add_choice`.

    Attributes
    ----------
    name : str
        Name of the choice.
    options : list of tuple of (str, str)
        Each option's value and menu title, in menu order.
    actions : dict of str to QAction
        The menu item of each value.
    """

    def __init__(self, manager: "PluginManager", info: PluginInfo, name: str, options: list[tuple[str, str]],
                 default: str | Callable[[DocumentView], str],
                 changed: Callable[[DocumentView, str], Any] | None,
                 enabled: Callable[[DocumentView, str], bool] | None,
                 read: Callable[[DocumentView], str] | None,
                 write: Callable[[DocumentView, str], Any] | None) -> None:
        self._manager = manager
        self.info = info
        self.name = name
        self.options = options
        self.actions: dict[str, QAction] = {}
        self._default = default
        self._changed = changed
        self._enabled = enabled
        self._read = read
        self._write = write

    def value(self, document: DocumentView) -> str:
        """Return the option chosen for a tab."""
        if self._read is not None:
            return str(self._read(document))
        state = self._manager.state(self.info, document)
        if self.name in state and str(state[self.name]) in self.actions:
            return str(state[self.name])
        return str(self._default(document) if callable(self._default) else self._default)

    def set(self, document: DocumentView, value: str) -> None:
        """Choose an option for a tab, as choosing its menu item does.

        Raises
        ------
        ValueError
            If the value is not one of the options.
        """
        if value not in self.actions:
            raise ValueError(f"{value!r} is not an option of {self.name}; "
                             f"the options are {', '.join(self.actions)}")
        if self._write is not None:
            self._write(document, value)
        else:
            self._manager.state(self.info, document)[self.name] = value
        if document == self._manager.app.current:
            self.sync(document)
        if self._changed is not None:
            try:
                self._changed(document, value)
            except Exception:
                self._manager.failed(self.info)

    def sync(self, document: DocumentView | None) -> None:
        """Check the option of a tab (the one in front) and grey out the options it cannot take."""
        if document is None:
            return
        try:
            chosen = self.value(document)
            for value, action in self.actions.items():
                action.setChecked(value == chosen)
                action.setEnabled(self._enabled is None or bool(self._enabled(document, value)))
        except Exception:
            self._manager.failed(self.info)


class Setting:
    """A setting of a plugin, edited in Settings and kept with the other settings.

    See :meth:`PluginContext.add_setting`.

    Attributes
    ----------
    name, title : str
        Key in the settings and text in the Settings window.
    default : bool, int, float or str
        Value until the user changes it; its type is the setting's type.
    minimum, maximum, step : float or None
        Range and step of a number.
    unit : str
        Shown after a number.
    """

    def __init__(self, manager: "PluginManager", info: PluginInfo, name: str, title: str,
                 default: bool | int | float | str, minimum: float | None, maximum: float | None,
                 step: float | None, unit: str, changed: Callable[[Any], Any] | None) -> None:
        self._manager = manager
        self.info = info
        self.name = name
        self.title = title
        self.default = default
        self.minimum, self.maximum, self.step = minimum, maximum, step
        self.unit = unit
        self._changed = changed

    @property
    def key(self) -> str:
        """Settings key, ``plugins/settings/<plugin key>/<name>``."""
        return f"plugins/settings/{self.info.key}/{self.name}"

    @property
    def value(self) -> Any:
        """The current value, of the type of :attr:`default`."""
        raw = self._manager.window.settings.value(self.key, None)
        if raw is None:
            return self.default
        try:
            if isinstance(self.default, bool):
                value: Any = str(raw).lower() in ("true", "1")
            elif isinstance(self.default, int):
                value = int(float(str(raw)))
            elif isinstance(self.default, float):
                value = float(str(raw))
            else:
                value = str(raw)
        except (TypeError, ValueError):
            return self.default
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if self.minimum is not None:
                value = max(value, type(self.default)(self.minimum))
            if self.maximum is not None:
                value = min(value, type(self.default)(self.maximum))
        return value

    def set(self, value: Any) -> None:
        """Store a new value and tell the plugin."""
        self._manager.window.settings.setValue(self.key, value)
        if self._changed is not None:
            try:
                self._changed(self.value)
            except Exception:
                self._manager.failed(self.info)

    def reset(self) -> None:
        """Go back to the default."""
        self._manager.window.settings.remove(self.key)
        if self._changed is not None:
            try:
                self._changed(self.default)
            except Exception:
                self._manager.failed(self.info)
