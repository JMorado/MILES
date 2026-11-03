"""2D drawings: a skeletal formula of the selection, or of everything, as SVG or PNG.

Adds File > Export 2D Drawing... and Edit > Copy 2D Drawing.
"""

import os

from PySide6.QtWidgets import QApplication

from miles import ChemError
from miles.builtin_plugins.drawing.depiction import depiction_svg
from miles.builtin_plugins.drawing.files import drawing_mime, write_drawing
from miles.plugin import App, PluginContext, Severity

NAME = "2D Drawings"
VERSION = "1.0"
DESCRIPTION = "Export or copy a 2D drawing of the selection or the whole structure"


def drawing(app: App) -> str | None:
    """SVG drawing of the current selection or everything, or None after telling the user why not."""
    doc = app.current
    if doc is None or not len(doc.molecule):
        return None
    try:
        return depiction_svg(doc.molecule.subset(doc.selection_with_hydrogens()), doc.name)
    except ChemError as exc:
        app.message(str(exc), severity=Severity.WARNING)
        return None


def export_drawing(app: App) -> None:
    """Save a 2D drawing of the selection or everything as SVG or PNG."""
    svg = drawing(app)
    doc = app.current
    if svg is None or doc is None:
        return
    name = os.path.splitext(doc.name)[0] or "structure"
    path = app.ask_save_path("Export 2D Drawing", name + ".svg", "SVG drawing (*.svg);;PNG image (*.png)")
    if path is None:
        return
    try:
        write_drawing(svg, path)
    except OSError as exc:
        app.warn("Could not export", str(exc))
        return
    app.message(f"Exported {os.path.basename(path)}")


def copy_drawing(app: App) -> None:
    """Put a 2D drawing of the selection or everything on the clipboard."""
    svg = drawing(app)
    if svg is not None:
        QApplication.clipboard().setMimeData(drawing_mime(svg))
        app.message("Copied 2D drawing")


def register(context: PluginContext) -> None:
    """Add the export and copy commands next to the other export and copy commands."""
    context.add_command("File > Export 2D Drawing...", export_drawing, before="scripts")
    context.add_command("Edit > Copy 2D Drawing", copy_drawing, after="Copy SMILES")
