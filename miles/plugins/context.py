"""What a plugin's ``register(context)`` receives: :class:`PluginContext`."""

import os
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

from PySide6.QtGui import QAction
from PySide6.QtWidgets import QMenu, QWidget

from miles.api import App, DocumentView
from miles.core.molecule import Molecule
from miles.files.common import JobInfo
from miles.plugins.discovery import PluginInfo
from miles.render.layers import Annotations, DrawContext, SceneItems
from miles.ui.canvas.extensions import MenuItem, MenuTarget
from miles.ui.lookup import StructureSource

if TYPE_CHECKING:
    from miles.core.document import Change
    from miles.plugins.background import BackgroundResults
    from miles.plugins.layers import PluginLayer
    from miles.plugins.manager import PluginManager
    from miles.plugins.settings import Choice, Setting, Toggle
    from miles.render.theme import Theme
    from miles.ui.floating import DocumentPanel, FloatingPanel

MenuAnchor = str | Sequence[str] | None
"""Where a menu item goes: the title of an item, or the name of one of the named separators of the menus,
or several of these tried in turn."""


class PluginContext:
    """What a plugin's ``register(context)`` receives to extend MILES.

    Attributes
    ----------
    app : miles.api.App
        The running application; callbacks receive it too.
    info : PluginInfo
        The plugin being registered.
    """

    def __init__(self, manager: "PluginManager", info: PluginInfo) -> None:
        self._manager = manager
        self.app = manager.app
        self.info = info

    def add_command(self, path: str, callback: Callable[[App], Any], shortcut: str | None = None,
                    enabled: Callable[[App], bool] | None = None, after: MenuAnchor = None,
                    before: MenuAnchor = None) -> QAction:
        """Add a menu command, also listed in the command palette and by :meth:`App.commands`.

        Parameters
        ----------
        path : str
            Where the command goes, menus separated by ``>``: ``"Tools >
            Charges > Compute"`` puts "Compute" in a "Charges" submenu of the
            Tools menu (created when missing). A path without a menu, just a
            title, goes to a "Plugins" menu; a new top-level menu goes before
            Help.
        callback : callable
            Called with the :class:`~miles.api.App` when the command
            runs, on the interface thread.
        shortcut : str, optional
            Key sequence such as ``"Ctrl+Shift+Y"`` (Ctrl is Cmd on macOS).
            Ignored, with a note in Help > Plugins, if another command has it.
        enabled : callable, optional
            Called with the App when the menu opens; the command is greyed
            out when it returns False, and does nothing then.
        after : str or sequence of str, optional
            Title of an item of the same menu to place the command after
            (default: at the end), or the name of a named separator, or
            several of these: the first one in the menu is used.
        before : str or sequence of str, optional
            Like ``after``, the item to place the command before; used when
            ``after`` is not given or not in the menu.

        Returns
        -------
        QAction
            The menu item.
        """
        return self._manager.add_command(self.info, path, callback, shortcut, enabled, after, before)

    def add_menu(self, path: str, after: MenuAnchor = None, before: MenuAnchor = None) -> QMenu:
        """Create a submenu at a chosen place (commands added to a path create missing menus at the end).

        Parameters
        ----------
        path : str
            Menus separated by ``>``, e.g. ``"View > Charges"``.
        after, before : str or sequence of str, optional
            Where the submenu goes, as for :meth:`add_command`.

        Returns
        -------
        QMenu
            The submenu (an existing one is returned as it is).
        """
        return self._manager.add_menu(path, after, before, self.info)

    def add_separator(self, path: str, after: MenuAnchor = None, before: MenuAnchor = None) -> None:
        """Add a separator line to a menu (at the end, or where ``after``/``before`` say)."""
        self._manager.add_separator(path, after, before, self.info)

    def add_toggle(self, path: str, name: str, changed: Callable[[DocumentView, bool], Any] | None = None,
                   default: bool | Callable[[DocumentView], bool] = False, shortcut: str | None = None,
                   after: MenuAnchor = None, before: MenuAnchor = None) -> "Toggle":
        """Add a checkable menu command whose state belongs to each tab, like View > Hydrogens.

        The check mark follows the tab in front, the state is saved in
        session files, and scripts can set it with
        ``doc.style(name=True)``.

        Parameters
        ----------
        path : str
            Menu path, as for :meth:`add_command`.
        name : str
            Name of the setting, unique among the toggles of all plugins.
        changed : callable, optional
            ``changed(document, on)`` called after the state of a tab changed.
        default : bool or callable, optional
            State of a tab that was never switched, or ``default(document)``
            working it out (e.g. from a session written by an older version).
        shortcut : str, optional
            Key sequence.
        after, before : str or sequence of str, optional
            Where the command goes.

        Returns
        -------
        Toggle
            Reads and sets the state of a tab.
        """
        return self._manager.add_toggle(self.info, path, name, changed, default, shortcut, after, before)

    def add_choice(self, menu: str, name: str, options: Mapping[str, str] | Sequence[str],
                   changed: Callable[[DocumentView, str], Any] | None = None,
                   default: str | Callable[[DocumentView], str] | None = None,
                   enabled: Callable[[DocumentView, str], bool] | None = None,
                   read: Callable[[DocumentView], str] | None = None,
                   write: Callable[[DocumentView, str], Any] | None = None,
                   after: MenuAnchor = None, before: MenuAnchor = None) -> "Choice":
        """Add menu items of which one is checked for each tab, like the colourings of View > Surface.

        By default the choice of each tab is kept with the plugin's state
        of the tab (saved in session files) and scripts can set it with
        ``doc.style(name="value")``; ``read`` and ``write`` keep it
        elsewhere instead, such as in the document itself.

        Parameters
        ----------
        menu : str
            Menu the items go in, e.g. ``"View > Surface"``.
        name : str
            Name of the choice, unique among the toggles and choices of all plugins.
        options : mapping or sequence of str
            The values with their menu titles (``{"gasteiger": "Gasteiger"}``),
            or values that are their own titles, in menu order.
        changed : callable, optional
            ``changed(document, value)`` called after a tab's choice changed.
        default : str or callable, optional
            Value of a tab never chosen, or ``default(document)``; the first
            option when not given.
        enabled : callable, optional
            ``enabled(document, value)``: whether an option can be chosen
            for the tab in front; checked when the menu opens.
        read, write : callable, optional
            ``read(document)`` returning the value and
            ``write(document, value)`` storing it, to keep it elsewhere.
        after, before : str or sequence of str, optional
            Where the first item goes, as for :meth:`add_command`.

        Returns
        -------
        Choice
            Reads and sets the choice of a tab.
        """
        items = list(options.items()) if isinstance(options, Mapping) else [(str(o), str(o)) for o in options]
        return self._manager.add_choice(self.info, menu, name, items, changed, default, enabled, read, write,
                                        after, before)

    def add_setting(self, name: str, title: str, default: bool | int | float | str,
                    minimum: float | None = None, maximum: float | None = None, step: float | None = None,
                    unit: str = "", changed: Callable[[Any], Any] | None = None) -> "Setting":
        """Add a setting that users change in Settings, such as a cutoff or a scale.

        The editor follows the type of ``default``: a check box for a bool,
        a number field for an int or a float, a text field for a string.
        The value is kept with the other settings.

        Parameters
        ----------
        name : str
            Key of the setting, unique within the plugin.
        title : str
            Text next to the editor.
        default : bool, int, float or str
            Value until the user changes it.
        minimum, maximum, step : float, optional
            Range and step of a number.
        unit : str, optional
            Shown after a number, e.g. ``"angstrom"``.
        changed : callable, optional
            ``changed(value)`` called after the user changed it.

        Returns
        -------
        Setting
            Its ``value`` is the current value.
        """
        return self._manager.add_setting(self.info, name, title, default, minimum, maximum, step, unit,
                                         changed)

    def setting_value(self, plugin: str, name: str, default: Any) -> Any:
        """Read the current value of a setting of another plugin, e.g. criteria that two plugins share.

        Parameters
        ----------
        plugin : str
            :attr:`~miles.plugins.PluginInfo.key` of the plugin, e.g.
            ``"miles.builtin_plugins.hydrogen_bonds"``.
        name : str
            Name of its setting.
        default : bool, int, float or str
            Value when the user never changed it (also its type).

        Returns
        -------
        object
            The value; ``default`` when that plugin is not loaded and the
            setting was never stored.
        """
        return self._manager.setting_value(plugin, name, default)

    def document_state(self, document: DocumentView) -> dict[str, Any]:
        """Return a dictionary this plugin keeps for a tab, saved in session files.

        Values must be JSON-compatible (numbers, strings, lists, dictionaries).
        The dictionary of a tab opened from a session holds what was saved.
        """
        return self._manager.state(self.info, document)

    def document_cache(self, document: DocumentView) -> dict[str, Any]:
        """Return a dictionary this plugin keeps for a tab, for working data that is not saved.

        Anything can go in it: arrays, meshes, results of slow work. It is
        dropped when the tab closes and when the plugin is switched off or
        reloaded, so the plugin needs no bookkeeping of its own.
        """
        return document.cache(self.info.key)

    def add_annotations(self, function: Callable[[DocumentView, DrawContext], Annotations | None]
                        ) -> "PluginLayer":
        """Draw lines and labels over the molecule, on screen, in exported pictures and in movies.

        Colours of lines, labels and arcs may be given in any form
        :func:`~miles.plugin.color` accepts: a name, a QColor or a tuple.

        Parameters
        ----------
        function : callable
            ``function(document, context)`` returning :class:`Annotations`
            for a tab, or None. ``context`` is a :class:`DrawContext`. On
            screen the result is kept until the structure, the selection,
            the frame, the display options or the theme change, or
            :meth:`PluginLayer.update` is called.

        Returns
        -------
        PluginLayer
            Call its ``update()`` when something the function depends on
            changed, such as a setting of the plugin.
        """
        return self._manager.add_layer(self.info, "annotations", function)

    def add_scene_items(self, function: Callable[[DocumentView, DrawContext], SceneItems | None]
                        ) -> "PluginLayer":
        """Draw arrows and meshes with the molecule, or colour its atoms, in both renderers.

        Parameters
        ----------
        function : callable
            ``function(document, context)`` returning :class:`SceneItems`
            for a tab, or None; called again under the same conditions as
            for :meth:`add_annotations`.

        Returns
        -------
        PluginLayer
            Call its ``update()`` when the items should be worked out again.
        """
        return self._manager.add_layer(self.info, "scene", function)

    def add_status(self, function: Callable[[DocumentView], str | None]) -> "PluginLayer":
        """Add text to the status line at the bottom left, after the formula, mass and point group.

        Parameters
        ----------
        function : callable
            ``function(document)`` returning the text, or an empty string or
            None for nothing. It is called a moment after the structure
            changes (not during playback or atom drags).

        Returns
        -------
        PluginLayer
            Call its ``update()`` to work the status line out again now.
        """
        return self._manager.add_layer(self.info, "status", function)

    def add_context_items(self, function: Callable[[DocumentView, MenuTarget], Iterable[MenuItem]],
                          section: str = "atom") -> None:
        """Add items to the canvas's right-click menu.

        Parameters
        ----------
        function : callable
            ``function(document, target)`` returning :class:`MenuItem`
            objects for what was clicked (:class:`MenuTarget`), possibly
            none; called each time the menu opens.
        section : str, optional
            ``"atom"`` (after the items for the atom clicked), ``"bond"``
            (after those for the bond clicked) or ``"selection"`` (after the
            measurement items, for the selection or the whole structure).
        """
        self._manager.add_context_items(self.info, function, section)

    def add_structure_source(self, source: StructureSource) -> None:
        """Let the fragment search, the command line, ``app.look_up`` and ``app.fetch`` fetch structures.

        See :class:`StructureSource`. Its ``fetch`` runs on a worker thread;
        errors other than ``FileFormatError`` are recorded as errors of the
        plugin.
        """
        self._manager.add_structure_source(self.info, source)

    def add_fragments(self, path: str | os.PathLike[str]) -> None:
        """Add the fragments of a JSON library (the format of ``core/data/fragments.json``) to the fragments.

        Categories ``"bio"`` and ``"ligand"`` add to the Biomolecules and
        Ligands tabs; any other category is the title of a tab of its own.

        Raises
        ------
        OSError, ValueError
            If the file cannot be read or is not a fragment library.
        """
        self._manager.add_fragments(self.info, os.fspath(path))

    def add_file_reader(self, extensions: str | Iterable[str],
                        reader: Callable[[str], tuple[Molecule, JobInfo]]) -> None:
        """Open files with these extensions with a reader of the plugin's own.

        Parameters
        ----------
        extensions : str or iterable of str
            Extensions without the dot, e.g. ``"toy"`` or ``["toy", "tyy"]``.
        reader : callable
            ``reader(path) -> (Molecule, JobInfo)``; raise
            :class:`~miles.files.common.FileFormatError` (or ValueError)
            with a message for the user when the file cannot be read.
        """
        self._manager.add_file_functions(self.info, extensions, reader=reader)

    def add_file_writer(self, extensions: str | Iterable[str],
                        writer: Callable[[str, Molecule, JobInfo], None]) -> None:
        """Save files with these extensions with a writer of the plugin's own.

        Parameters
        ----------
        extensions : str or iterable of str
            Extensions without the dot.
        writer : callable
            ``writer(path, molecule, info) -> None``.
        """
        self._manager.add_file_functions(self.info, extensions, writer=writer)

    def add_panel(self, title: str, factory: Callable[[App], QWidget], menu: str = "View",
                  shortcut: str | None = None) -> "FloatingPanel | None":
        """Add a floating panel on the right of the window, opened from a menu command.

        Parameters
        ----------
        title : str
            Title of the panel and of the command that shows or hides it.
        factory : callable
            ``factory(app) -> QWidget`` building the panel's contents once.
        menu : str, optional
            Menu path of the command (``"View"`` gives View > title).
        shortcut : str, optional
            Key sequence of the command.

        Returns
        -------
        FloatingPanel or None
            The panel (movable and foldable like the built-in panels), or None
            if the factory failed (the error is recorded). Its ``opened``
            and ``closed`` signals say when it appears and disappears, and
            ``open_panel()``, ``close_panel()`` and ``toggle()`` show and hide
            it from code. A panel about the tab in front is better made with
            :meth:`add_document_panel`.
        """
        return self._manager.add_panel(self.info, title, factory, menu, shortcut)

    def add_document_panel(self, name: str, panel: "DocumentPanel", command: str | None = None,
                           shortcut: str | None = None, place: str = "right", after: MenuAnchor = None,
                           before: MenuAnchor = None) -> "DocumentPanel":
        """Put a panel about the tab in front in the window, and keep it following that tab.

        The panel is told when another tab comes to the front
        (:meth:`~miles.ui.floating.DocumentPanel.set_document`), when that
        tab changes (:meth:`~miles.ui.floating.DocumentPanel.document_changed`),
        when its selection changes
        (:meth:`~miles.ui.floating.DocumentPanel.selection_changed`) and when
        the theme changes; its :meth:`~miles.ui.floating.DocumentPanel.store`
        and :meth:`~miles.ui.floating.DocumentPanel.saved` keep data for each
        tab. The command, when given, opens the panel or closes it when it
        is open, and says :attr:`~miles.ui.floating.DocumentPanel.unavailable`
        when the tab in front has nothing to show.

        Parameters
        ----------
        name : str
            Name its arrangement is saved under (see :meth:`add_floating_panel`).
        panel : DocumentPanel
            The panel.
        command : str, optional
            Menu path of the command that opens and closes it.
        shortcut : str, optional
            Key sequence of the command.
        place : str, optional
            ``"right"`` or ``"bottom"``, as for :meth:`add_floating_panel`.
        after, before : str or sequence of str, optional
            Where the command goes.

        Returns
        -------
        DocumentPanel
            The panel.
        """
        app = self.app
        panel.context = self
        self.add_floating_panel(name, panel, place)

        def follows(doc: DocumentView) -> bool:
            return doc == panel.document

        self.on_theme_changed(panel.set_theme)
        self.on_current_changed(panel.set_document)
        self.on_document_changed(lambda doc, change: panel.document_changed(change) if follows(doc) else None)
        self.on_selection_changed(lambda doc: panel.selection_changed() if follows(doc) else None)
        panel.set_document(app.current)
        if command is not None:
            def toggle(app: App) -> None:
                if app.current is None:
                    return
                was_open = panel.wanted and panel.isVisible()
                if not panel.toggle() and not was_open:   # it could not open: the tab has nothing to show
                    app.message(panel.unavailable)

            self.add_command(command, toggle, shortcut=shortcut, after=after, before=before)
        return panel

    def background_results(self, shows: "PluginLayer | Callable[[DocumentView], Any]", message: str,
                           name: str = "") -> "BackgroundResults":
        """Work out slow results for each tab on a worker thread, and show them when they arrive.

        See :class:`~miles.plugins.background.BackgroundResults`: its
        ``get(document, key, job, *args, now=drawing.export)`` returns the
        result for a tab, starting the work when the key is new.

        Parameters
        ----------
        shows : PluginLayer or callable
            The layer that draws the results (from :meth:`add_scene_items`
            or :meth:`add_annotations`), drawn again when a result arrives,
            or ``shows(document)`` to call then, e.g. to fill a panel.
        message : str
            What the busy bar says while the work runs, e.g. "Computing the surface".
        name : str, optional
            Tells several of these of one plugin apart.
        """
        from miles.plugins.background import BackgroundResults
        from miles.plugins.layers import PluginLayer

        ready = shows.update if isinstance(shows, PluginLayer) else shows
        return BackgroundResults(self, ready, f"{self.info.key}:{name}", message)

    def add_floating_panel(self, name: str, panel: "FloatingPanel", place: str = "right") -> None:
        """Put a panel of the plugin's own in the window, on the right below the others or at the bottom.

        This is for panels that need more than :meth:`add_panel` gives,
        such as buttons in the title row: a
        :class:`~miles.ui.floating.FloatingPanel` subclass that lays
        out its contents and shows and hides itself, opened by a command
        the plugin adds. It folds and moves like the other panels, and its
        place and fold are remembered under ``panels/<name>`` in the
        settings.

        Parameters
        ----------
        name : str
            Name its arrangement is saved under; no other panel may have it.
        panel : FloatingPanel
            The panel; it is hidden and reparented to the window.
        place : str, optional
            ``"right"`` (the default) stacks it on the right; ``"bottom"``
            centres it above the timeline, like the grid panel, and fitting
            the view keeps the structure clear of it.

        Raises
        ------
        ValueError
            If another panel has that name, or the place is unknown.
        """
        self._manager.add_floating_panel(name, panel, place, self.info)

    def on_document_opened(self, callback: Callable[[DocumentView], Any]) -> None:
        """Call ``callback(document)`` for each new tab (opened files, new documents).

        Tabs that were open before the plugin loaded, such as the empty
        first one, are in ``context.app.documents``.
        """
        self._manager.listen("opened", self.info, callback)

    def on_document_closed(self, callback: Callable[[DocumentView], Any]) -> None:
        """Call ``callback(document)`` just before a tab closes, e.g. to forget what was kept for it."""
        self._manager.listen("closed", self.info, callback)

    def on_current_changed(self, callback: Callable[[DocumentView], Any]) -> None:
        """Call ``callback(document)`` when another tab comes to the front, e.g. to refresh a panel."""
        self._manager.listen("current", self.info, callback)

    def on_selection_changed(self, callback: Callable[[DocumentView], Any]) -> None:
        """Call ``callback(document)`` when the selection of a tab changes."""
        self._manager.listen("selection", self.info, callback)

    def on_structure_changed(self, callback: Callable[[DocumentView], Any]) -> None:
        """Call ``callback(document)`` after edits of a tab's atoms, bonds or coordinates.

        Frame changes during playback and style changes do not count, so
        the callback is not called many times a second during animations;
        it is during atom drags, so keep it quick.
        """
        self._manager.listen("structure", self.info, callback)

    def on_frame_changed(self, callback: Callable[[DocumentView], Any]) -> None:
        """Call ``callback(document)`` when a tab shows another frame, or its frames are replaced or removed.

        During playback this happens many times a second, so keep the
        callback quick.
        """
        self._manager.listen("frame", self.info, callback)

    def on_document_changed(self, callback: Callable[[DocumentView, "Change"], Any]) -> None:
        """Call ``callback(document, change)`` after every change of a tab, saying what changed.

        ``change`` is a :class:`~miles.core.document.Change` (a string
        enumeration, so ``change == "coords"`` works too):

        - ``STRUCTURE``: atoms or bonds added or removed;
        - ``BONDS``: bond orders only;
        - ``COORDS``: positions only;
        - ``STYLE``: per-atom display (representation, hidden atoms) or measurements;
        - ``PROPERTIES``: per-atom properties such as partial charges, or the charge model;
        - ``FRAME``: another trajectory frame (or a point between two) is shown;
        - ``FRAMES``: the frames were replaced or removed, or an animation started or stopped;
        - ``NAMES``: named selections were defined or removed;
        - ``RESULTS``: grids were added or removed, or results read from the
          file (orbitals, modes, grids, charges) stopped or started to apply.

        The structure and frame events are filtered forms of this one; use
        it when a panel must tell changes apart, for instance to skip work
        that moving atoms does not affect. ``document.last_change`` says
        which atoms and properties changed, when the editor said so.
        """
        self._manager.listen("changed", self.info, callback)

    def on_unload(self, callback: Callable[[], Any]) -> None:
        """Call ``callback()`` when the plugin is switched off or reloaded, to free what MILES does not track.

        Commands, panels, layers, settings and listeners the plugin added
        are removed by MILES itself; this is for timers, threads, files or
        connections the plugin made on its own.
        """
        self._manager.on_unload(self.info, callback)

    def on_theme_changed(self, callback: Callable[["Theme"], Any]) -> None:
        """Call ``callback(theme)`` at once with the current theme, and again whenever the appearance changes.

        The :class:`~miles.render.theme.Theme` holds the colours that
        panels and plots drawn by the plugin should use.
        """
        callback(self.app.theme)
        self._manager.listen("theme", self.info, callback)
