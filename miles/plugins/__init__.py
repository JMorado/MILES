"""Plugins: extending MILES with plain Python, without changing its code.

A plugin is a Python module with a ``register(context)`` function. MILES
calls it once at start-up with a :class:`PluginContext`, through which the
plugin adds menu commands and toggles, file readers and writers, panels,
annotations and 3D items drawn with the molecule, status text, right-click
menu items, structure sources, fragment libraries, and callbacks for
events. Optional module attributes describe it: ``NAME``, ``VERSION``,
``DESCRIPTION`` and ``REQUIRES_API``.

Plugins import what they use from :mod:`miles.plugin`, the stable surface
for plugin authors. This package is the machinery behind it:

:mod:`miles.plugins.discovery`
    Finding plugins and what is recorded about each (:class:`PluginInfo`).
:mod:`miles.plugins.context`
    What ``register(context)`` receives (:class:`PluginContext`).
:mod:`miles.plugins.settings`
    Per-tab toggles and choices, and settings (:class:`Toggle`,
    :class:`Choice`, :class:`Setting`).
:mod:`miles.plugins.layers`
    Drawings and status text worked out again on demand (:class:`PluginLayer`).
:mod:`miles.plugins.background`
    Slow results worked out on a worker thread for each tab
    (:class:`BackgroundResults`).
:mod:`miles.plugins.dialogs`
    Dialogs in the colours of the interface (:class:`PluginDialog`).
:mod:`miles.plugins.colors` and :mod:`miles.plugins.text`
    Colours in any form, and atom labels with their residue.
:mod:`miles.plugins.manager`
    Loading plugins and connecting them to the window (:class:`PluginManager`).

Plugins are found in three places, loaded in this order:

1. The built-in plugins (:data:`BUILTIN`), features of the application
   written against the same interface;
2. installed packages that declare an entry point in the
   :data:`ENTRY_POINT_GROUP` group;
3. ``.py`` files and packages in the user's plugin folder
   (:func:`plugin_folder`), so a single file can be dropped in without
   packaging it.

A plugin that fails to import or to register, or whose callback raises,
never stops the application: the error is recorded and shown under
Help > Plugins, where plugins can also be switched off and reloaded.
"""

from miles.plugins.context import MenuAnchor, PluginContext
from miles.plugins.discovery import (
    BUILTIN,
    ENTRY_POINT_GROUP,
    FOLDER_ENV,
    NO_PLUGINS_ENV,
    PLUGIN_API,
    PluginInfo,
    PluginState,
    discover,
    plugin_folder,
)
from miles.plugins.layers import PluginLayer
from miles.plugins.manager import DEFAULT_MENU, DISABLED_KEY, PluginManager
from miles.plugins.settings import Setting, Toggle
from miles.render.layers import Annotations, Arc, DrawContext, Label, Line, MeshItem, SceneItems
from miles.ui.canvas.extensions import MenuItem, MenuTarget
from miles.ui.lookup import StructureSource

__all__ = [
    "BUILTIN",
    "DEFAULT_MENU",
    "DISABLED_KEY",
    "ENTRY_POINT_GROUP",
    "FOLDER_ENV",
    "NO_PLUGINS_ENV",
    "PLUGIN_API",
    "Annotations",
    "Arc",
    "DrawContext",
    "Label",
    "Line",
    "MenuAnchor",
    "MenuItem",
    "MenuTarget",
    "MeshItem",
    "PluginContext",
    "PluginInfo",
    "PluginLayer",
    "PluginManager",
    "PluginState",
    "SceneItems",
    "Setting",
    "StructureSource",
    "Toggle",
    "discover",
    "plugin_folder",
]

# The classes and functions of the submodules are documented, and named in
# tracebacks and the documentation, as members of this package.
for _public in (PluginContext, PluginInfo, PluginState, PluginLayer, PluginManager, Setting, Toggle,
                discover, plugin_folder):
    _public.__module__ = __name__
del _public
