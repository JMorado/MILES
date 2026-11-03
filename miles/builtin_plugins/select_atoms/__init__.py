"""Select atoms with an expression of the selection language, such as ``byres within 5 of ligand``.

Adds Edit > Select... (Cmd F). The language itself is part of the library
(:func:`miles.find <miles.core.selection.find>`); SMARTS patterns are part of it
(``smarts "[OX2H]"``) and the common ones are offered as examples.
"""

from miles.builtin_plugins.select_atoms.dialog import SelectionDialog
from miles.plugin import App, PluginContext, plural

NAME = "Select by Expression"
VERSION = "1.0"
DESCRIPTION = "Select atoms with expressions such as byres within 5 of ligand, or with SMARTS patterns"


def select_atoms(app: App) -> None:
    """Ask for an expression and select the atoms it picks."""
    doc = app.current
    if doc is None or not len(doc.molecule):
        return
    app.set_tool("select")
    dialog = SelectionDialog(app, doc)
    if dialog.run() and dialog.atoms:
        doc.selection = dialog.atoms
        app.message(f"Selected {plural(len(dialog.atoms), 'atom')}")


def register(context: PluginContext) -> None:
    """Add the command after Edit > Select All."""
    context.add_command("Edit > Select...", select_atoms, shortcut="Ctrl+F", after="Select All")
