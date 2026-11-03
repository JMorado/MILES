"""Properties: composition, identifiers and drug-likeness of the structure or the selection.

The panel follows the current tab, the selection and edits of atoms and
bonds; clicking a value copies it. Adds Tools > Properties.
"""

from miles.builtin_plugins.properties.panel import PropertiesPanel
from miles.plugin import PluginContext

NAME = "Properties"
VERSION = "1.0"
DESCRIPTION = "Formula, masses, identifiers and drug-likeness descriptors, ready to copy"


def register(context: PluginContext) -> None:
    """Add the panel, and the command that opens and closes it after Tools > Superpose..."""
    context.add_document_panel("properties", PropertiesPanel(), "Tools > Properties", shortcut="Ctrl+I",
                               before=("Vibrations", "protein tools"))
