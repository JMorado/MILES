Meshes and surfaces (provisional)
#################################

What builds the meshes and arrows of :class:`~miles.plugin.SceneItems`: ``from miles.plugin.render import ...``. Provisional: these follow the renderer of MILES and may change between releases (see :ref:`dg-stability`).

.. automodule:: miles.plugin.render
   :no-members:

Meshes and arrows
*****************

Triangles and arrows drawn with the molecule.

.. autoclass:: miles.plugin.render.Mesh
   :no-index:

.. autofunction:: miles.plugin.render.arrow_rows
   :no-index:

.. autofunction:: miles.plugin.render.charge_colors
   :no-index:

.. autoclass:: miles.plugin.render.LabelMode
   :no-index:

.. data:: miles.plugin.render.RIBBON_ATOMS
   :no-index:

   Names of the backbone atoms that a cartoon replaces.

Surfaces and isosurfaces
************************

The molecular surface and the isosurfaces of grids.

.. autofunction:: miles.plugin.render.surface_mesh
   :no-index:

.. autoclass:: miles.plugin.render.SurfaceRequest
   :no-index:

.. autoclass:: miles.plugin.render.SurfaceColoring
   :no-index:

.. autofunction:: miles.plugin.render.ionisable_charges
   :no-index:

.. autofunction:: miles.plugin.render.volume_mesh
   :no-index:

.. autoclass:: miles.plugin.render.VolumeRequest
   :no-index:

.. data:: miles.plugin.render.MAPPED_COLORS
   :no-index:

   Colours of a surface coloured by the values of a second grid, from low to high.

.. autofunction:: miles.plugin.render.symmetric_range
   :no-index:
