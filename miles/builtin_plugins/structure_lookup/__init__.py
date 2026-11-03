"""Structure Lookup: PDB entries by their code and compounds by name from PubChem.

Adds two structure sources, used by the fragment search, the command line,
``app.look_up`` and ``app.fetch``, and File > Open from PDB... A PDB entry
opens in a tab of its own; a compound becomes the fragment to place.
"""

from typing import override

from PySide6.QtWidgets import QLineEdit

from miles import JobInfo, Molecule
from miles.builtin_plugins.structure_lookup.remote import fetch_compound, fetch_pdb, is_pdb_code
from miles.plugin import App, PluginContext, PluginDialog, StructureSource
from miles.plugin.widgets import label

NAME = "Structure Lookup"
VERSION = "1.0"
DESCRIPTION = "Download PDB entries and look up compounds by name on PubChem"


def compound(name: str) -> tuple[Molecule, JobInfo]:
    """3D structure of a compound by name, from PubChem."""
    return fetch_compound(name), JobInfo(title=name)


class PdbDialog(PluginDialog):
    """A field for a PDB code; :attr:`code` holds it after :meth:`run` accepts."""

    def __init__(self, app: App) -> None:
        super().__init__(app, "Open from PDB")
        self.body.addWidget(label("PDB code", "Title"))
        self.field = QLineEdit()
        self.field.setObjectName("Field")
        self.field.setPlaceholderText("e.g. 1crn, 4hhb")
        self.field.setMinimumWidth(260)
        self.body.addWidget(self.field)
        self.message = label("The entry opens in a tab of its own.", "Muted")
        self.body.addWidget(self.message)
        self.button = self.buttons("Open")
        self.field.textChanged.connect(lambda text: self.button.setEnabled(is_pdb_code(text.strip())))
        self.button.setEnabled(False)
        self.code = ""

    @override
    def accept(self) -> None:
        """Close with a code that looks like a PDB code."""
        code = self.field.text().strip()
        if is_pdb_code(code):
            self.code = code
            super().accept()


def open_from_pdb(app: App) -> None:
    """Ask for a PDB code and download the entry."""
    dialog = PdbDialog(app)
    if dialog.run():
        app.look_up(dialog.code)


def register(context: PluginContext) -> None:
    """Add the PubChem and PDB sources and File > Open from PDB..."""
    context.add_structure_source(StructureSource("PubChem", compound, what="name", example="caffeine",
                                                 whole=False))
    context.add_structure_source(StructureSource(
        "PDB", fetch_pdb, matches=is_pdb_code, what="PDB code", example="1crn",
        message=lambda text: f"Fetching {text.upper()} from the PDB ..."))
    context.add_command("File > Open from PDB...", open_from_pdb, shortcut="Ctrl+Shift+O", after="Open...")
