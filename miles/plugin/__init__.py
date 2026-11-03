"""Everything a plugin needs to work with the window: the stable plugin interface.

A plugin imports the window side from here and the chemistry it computes
with from :mod:`miles` itself::

    import miles
    from miles.plugin import App, PluginContext, plural

    def count(app: App) -> None:
        if app.current is not None:
            app.message(plural(len(app.current.molecule), "atom"))

    def register(context: PluginContext) -> None:
        context.add_command("Tools > Count Atoms", count)

This module and :mod:`miles` are the stable surface: a release that
changes them in a way that breaks plugins raises ``PLUGIN_API``, and a
plugin that sets ``REQUIRES_API`` fails with a clear message instead of
misbehaving. Two submodules are provisional and may still change between
releases: :mod:`miles.plugin.widgets` (controls, colours and plot menus for
panels) and :mod:`miles.plugin.render` (meshes, arrows and surfaces).
Everything else in MILES is internal.

What it provides, by purpose:

Registering
    ``PluginContext`` (what ``register(context)`` receives), ``PluginInfo``,
    ``PluginState``, ``Setting``, ``Toggle``, ``Choice``, ``PluginLayer``,
    ``BackgroundResults`` (slow drawings and panels), ``MenuAnchor`` and
    ``PLUGIN_API``.
The application and its documents
    ``App``, ``DocumentView``, ``Picture``, ``current_app``, ``Severity``
    (of messages), the kinds of change a document announces, ``Change``
    and ``ChangeEvent``, and ``AtomDisplay``, the display of a document's
    atoms that ``DocumentView.edit_display`` hands to its function.
Drawing with the molecule
    ``Annotations`` with ``Line``, ``Label`` and ``Arc``; ``SceneItems``
    with ``MeshItem``; ``DrawContext``; ``color``, which reads a colour in
    any of the forms the drawing types accept.
Menus and lookups
    ``MenuItem``, ``MenuTarget`` and ``StructureSource``.
Panels and dialogs
    ``FloatingPanel``, ``DocumentPanel`` (a panel about the tab in
    front), ``PluginDialog`` (a dialog in the colours of the interface),
    ``Theme`` and ``read_filter`` (the file types for an open dialog).
Text
    ``plural``, ``subscript_digits``, ``atom_name`` (``A/ASP25:OD1``),
    ``ANGSTROM``, ``DEGREE`` and ``SEPARATOR`` (the dot between parts of a
    message).
"""

from miles.api import App, DocumentView, Picture, Severity, current_app
from miles.core.display import AtomDisplay
from miles.core.document import Change, ChangeEvent
from miles.files.formats import read_filter
from miles.plugins import (
    PLUGIN_API,
    Annotations,
    Arc,
    DrawContext,
    Label,
    Line,
    MenuAnchor,
    MenuItem,
    MenuTarget,
    MeshItem,
    PluginContext,
    PluginInfo,
    PluginLayer,
    PluginState,
    SceneItems,
    Setting,
    StructureSource,
    Toggle,
)
from miles.plugins.background import BackgroundResults
from miles.plugins.colors import color
from miles.plugins.dialogs import PluginDialog
from miles.plugins.settings import Choice
from miles.plugins.text import atom_name
from miles.render.theme import Theme
from miles.ui.floating import DocumentPanel, FloatingPanel
from miles.ui.text import ANGSTROM, DEGREE, SEPARATOR, plural, subscript_digits

__all__ = [
    "ANGSTROM",
    "DEGREE",
    "PLUGIN_API",
    "SEPARATOR",
    "Annotations",
    "App",
    "Arc",
    "AtomDisplay",
    "BackgroundResults",
    "Change",
    "ChangeEvent",
    "Choice",
    "DocumentPanel",
    "DocumentView",
    "DrawContext",
    "FloatingPanel",
    "Label",
    "Line",
    "MenuAnchor",
    "MenuItem",
    "MenuTarget",
    "MeshItem",
    "Picture",
    "PluginContext",
    "PluginDialog",
    "PluginInfo",
    "PluginLayer",
    "PluginState",
    "SceneItems",
    "Setting",
    "Severity",
    "StructureSource",
    "Theme",
    "Toggle",
    "atom_name",
    "color",
    "current_app",
    "plural",
    "read_filter",
    "subscript_digits",
]
