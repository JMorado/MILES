"""Build from Sequence: a peptide, DNA or RNA built from its one-letter sequence, in a new tab.

Adds Tools > Build from Sequence...
"""

from miles import Molecule
from miles.builtin_plugins.sequence.biopolymer import build_sequence
from miles.builtin_plugins.sequence.dialog import SequenceDialog
from miles.plugin import App, PluginContext, plural

NAME = "Build from Sequence"
VERSION = "1.0"
DESCRIPTION = "Build a peptide, DNA or RNA from its sequence"


def build_from_sequence(app: App) -> None:
    """Ask for a sequence, build it in the background and open the structure in a new tab."""
    dialog = SequenceDialog(app)
    if not dialog.run():
        return
    title, sequence, kind = dialog.title, dialog.sequence_text, dialog.chosen_kind()

    def built(molecule: Molecule) -> None:
        app.new(molecule, title)

    message = f"Building {plural(len(sequence), 'residue')}"
    app.run_in_background(build_sequence, sequence, kind, message=message, done=built)


def register(context: PluginContext) -> None:
    """Add the command after Tools > Fragments and SMILES..."""
    context.add_command("Tools > Build from Sequence...", build_from_sequence,
                        before=("Vibrations", "protein tools"))
