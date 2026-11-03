"""Export Movie: a turntable or the trajectory of the current view, as MP4, GIF or numbered PNG files.

Adds File > Export Movie... Every frame is drawn like an exported picture,
with what plugins add to it.
"""

import os

from miles.builtin_plugins.movie.dialog import MovieDialog
from miles.plugin import App, PluginContext

NAME = "Export Movie"
VERSION = "1.0"
DESCRIPTION = "Export a turntable or trajectory movie of the current view"


def export_movie(app: App) -> None:
    """Ask how to make the movie and write it."""
    doc = app.current
    if doc is None or not len(doc.molecule):
        return
    dialog = MovieDialog(app, doc)
    if dialog.run() and dialog.saved_path:
        app.message(f"Exported {os.path.basename(dialog.saved_path)}")


def register(context: PluginContext) -> None:
    """Add the command after File > Export Image..."""
    context.add_command("File > Export Movie...", export_movie, after="Export Image...")
