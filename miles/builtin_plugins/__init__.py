"""Features of MILES written as plugins.

Each module or package here is an ordinary plugin: a ``register(context)``
function that adds its commands, panels and drawings through
:class:`miles.plugin.PluginContext <miles.plugins.PluginContext>` and
imports only :mod:`miles` and :mod:`miles.plugin`. A package keeps the panels, dialogs and computations
that only its plugin uses next to the registration, so each feature is in
one place. They are loaded first, in
the order of :data:`miles.plugins.discovery.BUILTIN`, and are declared as entry
points of the ``miles.plugins`` group like any third-party plugin, so
they double as worked examples of the plugin interface.
"""
