Panel widgets (provisional)
###########################

Controls, colours and plot menus for the panels and dialogs of plugins: ``from miles.plugin.widgets import ...``. Provisional: these follow the panels of MILES and may change between releases (see :ref:`dg-stability`).

.. automodule:: miles.plugin.widgets
   :no-members:

Controls
********

The pieces the built-in panels are made of.

.. autofunction:: miles.plugin.widgets.button
   :no-index:

.. autofunction:: miles.plugin.widgets.label
   :no-index:

.. autofunction:: miles.plugin.widgets.muted
   :no-index:

.. autofunction:: miles.plugin.widgets.rule
   :no-index:

.. autoclass:: miles.plugin.widgets.Segmented
   :no-index:

.. autoclass:: miles.plugin.widgets.ScrubField
   :no-index:

.. autoclass:: miles.plugin.widgets.Dial
   :no-index:

Colours
*******

The two themes and the stylesheets of widgets.

.. data:: miles.plugin.widgets.LIGHT
   :no-index:

   The light theme.

.. data:: miles.plugin.widgets.DARK
   :no-index:

   The dark theme.

.. autofunction:: miles.plugin.widgets.rgb
   :no-index:

.. autofunction:: miles.plugin.widgets.dialog_stylesheet
   :no-index:

.. autofunction:: miles.plugin.widgets.theme_palette
   :no-index:

Plot menus
**********

The right-click menu of plots, and the data behind a plot (see :meth:`miles.plugin.App.export_plot`).

.. autofunction:: miles.plugin.widgets.show_plot_menu
   :no-index:

.. data:: miles.plugin.widgets.COPY
   :no-index:

   Plot menu entry that copies the plot as a picture.

.. data:: miles.plugin.widgets.DATA
   :no-index:

   Plot menu entry that saves the data behind a plot as a table.

.. data:: miles.plugin.widgets.CURVE
   :no-index:

   Plot menu entry that saves the curve behind a plot.

.. data:: miles.plugin.widgets.IMAGE
   :no-index:

   Plot menu entry that saves the plot as an image file.

.. autoclass:: miles.plugin.widgets.PlotSource
   :no-index:

.. autoclass:: miles.plugin.widgets.PaintablePlot
   :no-index:

.. autoclass:: miles.plugin.widgets.Table
   :no-index:

.. autofunction:: miles.plugin.widgets.format_number
   :no-index:
