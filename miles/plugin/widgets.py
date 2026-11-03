"""Controls, colours and plot menus for the panels and dialogs of plugins (provisional).

These are the pieces the built-in panels are made of. They are public, but
provisional: they may change between releases while the panels of MILES
change, so a plugin that uses them may need small updates. The stable
interface is :mod:`miles.plugin`.

Controls
    ``button``, ``label``, ``muted``, ``rule``, ``Segmented`` (a row of
    choices), ``ScrubField`` (a number dragged sideways) and ``Dial``.
Colours
    ``LIGHT`` and ``DARK`` (the two themes), ``rgb`` (a colour as three
    floats), ``dialog_stylesheet`` and ``theme_palette`` (for widgets drawn
    natively, such as combo boxes).
Plot menus
    ``show_plot_menu`` with ``COPY``, ``DATA``, ``CURVE`` and ``IMAGE``,
    ``PlotSource`` and ``PaintablePlot`` (see ``App.export_plot``), and
    ``Table`` with ``format_number`` for the data behind a plot.
"""

from miles.core.tables import Table, format_number
from miles.render.theme import DARK, LIGHT, rgb
from miles.ui.controls import (
    Dial,
    ScrubField,
    Segmented,
    button,
    dialog_stylesheet,
    label,
    rule,
    theme_palette,
)
from miles.ui.dialogs import muted
from miles.ui.plot_export import COPY, CURVE, DATA, IMAGE, PaintablePlot, PlotSource, show_plot_menu

__all__ = [
    "COPY",
    "CURVE",
    "DARK",
    "DATA",
    "IMAGE",
    "LIGHT",
    "Dial",
    "PaintablePlot",
    "PlotSource",
    "ScrubField",
    "Segmented",
    "Table",
    "button",
    "dialog_stylesheet",
    "format_number",
    "label",
    "muted",
    "rgb",
    "rule",
    "show_plot_menu",
    "theme_palette",
]
