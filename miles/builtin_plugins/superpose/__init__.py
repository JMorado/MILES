"""Superpose: move the current structure onto the structure of another tab.

Adds Tools > Superpose...
"""

from miles.builtin_plugins.superpose.dialog import SuperposeDialog, describe, superpose_documents
from miles.plugin import App, PluginContext, Severity

NAME = "Superpose"
VERSION = "1.0"
DESCRIPTION = "Superpose the current structure onto the structure of another tab"


def superpose(app: App) -> None:
    """Ask for the other tab and superpose, as an undoable edit."""
    doc = app.current
    if doc is None:
        return
    targets = [other for other in app.documents if other != doc and len(other.molecule)]
    if not targets:
        app.message("Open the structure to superpose onto in another tab")
        return
    dialog = SuperposeDialog(app, doc, targets)
    if not dialog.run():
        return
    try:
        result = superpose_documents(doc, dialog.chosen, dialog.copy.isChecked())
    except ValueError as exc:
        app.message(str(exc), severity=Severity.WARNING)
        return
    app.message(describe(result), severity=Severity.WARNING)


def register(context: PluginContext) -> None:
    """Add the command after Tools > Align Frames."""
    context.add_command("Tools > Superpose...", superpose, before=("Vibrations", "protein tools"))
