API Reference
=============

Reference documentation generated from the docstrings, in two parts.

Public API
----------

What scripts and plugins use: the chemistry library ``miles`` and the plugin
module ``miles.plugin``, which are stable and versioned together, and the two
provisional submodules of ``miles.plugin`` for panels and meshes, and
``miles.testing`` for a plugin's tests, which may still change between
releases (see :ref:`dg-stability`).

.. toctree::
   :maxdepth: 2

   library
   plugin
   plugin_widgets
   plugin_render
   testing

Contributing: internals
-----------------------

Every other module, for people working on MILES itself. They may change in
any release. The application has four layers, and a layer only imports the
layers listed above it; plugins come last and use only the public API.

.. toctree::
   :maxdepth: 2

   core
   files
   render
   ui
   scripting
   plugins
