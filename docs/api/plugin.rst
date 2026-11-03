The plugin module
#################

What a plugin needs to work with the window: ``from miles.plugin import ...``. :doc:`../developerguide/plugins` explains how to write one. This module is stable; its two submodules, :doc:`plugin_widgets` and :doc:`plugin_render`, are provisional (see :ref:`dg-stability`).

.. automodule:: miles.plugin
   :no-members:

Registering
***********

What ``register(context)`` receives and the settings, toggles, choices, drawing layers and slow results it creates.

.. autoclass:: miles.plugin.PluginContext
   :no-index:

.. autoclass:: miles.plugin.PluginInfo
   :no-index:

.. autoclass:: miles.plugin.PluginState
   :no-index:

.. autoclass:: miles.plugin.Setting
   :no-index:

.. autoclass:: miles.plugin.Toggle
   :no-index:

.. autoclass:: miles.plugin.Choice
   :no-index:

.. autoclass:: miles.plugin.PluginLayer
   :no-index:

.. autoclass:: miles.plugin.BackgroundResults
   :no-index:

.. data:: miles.plugin.MenuAnchor
   :no-index:

   Where a menu item goes: the title of an item, or the name of one of the named separators.

.. data:: miles.plugin.PLUGIN_API
   :no-index:

   Version of the plugin interface; a plugin that needs a newer one sets ``REQUIRES_API``.

The application and its documents
*********************************

The window as a script sees it, how it tells the user something, and what a document announces.

.. autoclass:: miles.plugin.App
   :no-index:

.. autoclass:: miles.plugin.DocumentView
   :no-index:

.. autoclass:: miles.plugin.Picture
   :no-index:

.. autofunction:: miles.plugin.current_app
   :no-index:

.. autoclass:: miles.plugin.Severity
   :no-index:

.. autoclass:: miles.plugin.Change
   :no-index:

.. autoclass:: miles.plugin.ChangeEvent
   :no-index:

.. autoclass:: miles.plugin.AtomDisplay
   :no-index:

Drawing with the molecule
*************************

Lines, labels, arcs and meshes drawn with the structure, and the colours they accept. The meshes themselves are built with :doc:`plugin_render`.

.. autoclass:: miles.plugin.Annotations
   :no-index:

.. autoclass:: miles.plugin.Line
   :no-index:

.. autoclass:: miles.plugin.Label
   :no-index:

.. autoclass:: miles.plugin.Arc
   :no-index:

.. autoclass:: miles.plugin.SceneItems
   :no-index:

.. autoclass:: miles.plugin.MeshItem
   :no-index:

.. autoclass:: miles.plugin.DrawContext
   :no-index:

.. autofunction:: miles.plugin.color
   :no-index:

Menus and lookups
*****************

Context menu items and structure sources.

.. autoclass:: miles.plugin.MenuItem
   :no-index:

.. autoclass:: miles.plugin.MenuTarget
   :no-index:

.. autoclass:: miles.plugin.StructureSource
   :no-index:

Panels and dialogs
******************

Floating panels, panels about the tab in front, dialogs in the colours of the interface, and the theme. The controls inside them are in :doc:`plugin_widgets`.

.. autoclass:: miles.plugin.FloatingPanel
   :no-index:

.. autoclass:: miles.plugin.DocumentPanel
   :no-index:

.. autoclass:: miles.plugin.PluginDialog
   :no-index:

.. autoclass:: miles.plugin.Theme
   :no-index:

.. autofunction:: miles.plugin.read_filter
   :no-index:

Text
****

Words and signs used in labels and messages.

.. autofunction:: miles.plugin.plural
   :no-index:

.. autofunction:: miles.plugin.subscript_digits
   :no-index:

.. autofunction:: miles.plugin.atom_name
   :no-index:

.. data:: miles.plugin.ANGSTROM
   :no-index:

   The angstrom sign, for labels and messages.

.. data:: miles.plugin.DEGREE
   :no-index:

   The degree sign.

.. data:: miles.plugin.SEPARATOR
   :no-index:

   The dot that separates the parts of a message.
