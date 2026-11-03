"""Ramachandran Plot: backbone phi and psi of every residue, linked to the 3D selection.

The panel follows the current tab, edits, trajectory frames and the
selection; its right-click menu copies or exports the angles and the
picture. Adds Tools > Ramachandran Plot.
"""

from miles.builtin_plugins.ramachandran.panel import RamachandranPanel
from miles.plugin import PluginContext

NAME = "Ramachandran Plot"
VERSION = "1.0"
DESCRIPTION = "Plot the backbone phi and psi angles of a protein, linked to the selection"


def register(context: PluginContext) -> None:
    """Add the panel, and the command that opens and closes it before the tool choice in Tools."""
    app = context.app
    panel = RamachandranPanel()
    context.add_document_panel("ramachandran", panel, "Tools > Ramachandran Plot", shortcut="Ctrl+Alt+R",
                               before="tool choice")
    panel.exportRequested.connect(lambda kind: app.export_plot(panel, kind))
    panel.toolRequested.connect(app.set_tool)
