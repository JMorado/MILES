"""Prepare Protein: complete side chains, add hydrogens for a pH and remove water, as one undoable step.

Adds Tools > Prepare Protein...
"""

from miles import Molecule, has_polymer
from miles.builtin_plugins.prepare_protein.dialog import PrepareDialog
from miles.builtin_plugins.prepare_protein.preparation import PreparationReport, prepare
from miles.plugin import App, PluginContext, Severity, plural

NAME = "Prepare Protein"
VERSION = "1.0"
DESCRIPTION = "Complete side chains, add hydrogens for a pH and remove water"


def prepare_protein(app: App) -> None:
    """Ask what to do and prepare the protein or nucleic acid of the current tab."""
    doc = app.current
    if doc is None:
        return
    if not has_polymer(doc.molecule):
        app.message("Prepare Protein works on proteins and nucleic acids")
        return
    dialog = PrepareDialog(app, doc.molecule)
    if not dialog.run():
        return
    complete, hydrogens = dialog.complete.isChecked(), dialog.hydrogens.isChecked()
    ph, waters = dialog.ph.value(), dialog.water.isChecked()
    revision = doc.revision

    def job(mol: Molecule) -> tuple[Molecule, PreparationReport]:
        return mol, prepare(mol, complete=complete, hydrogens=hydrogens, ph=ph, waters=waters)

    def done(result: tuple[Molecule, PreparationReport]) -> None:
        if doc not in app.documents:
            return
        if doc.revision != revision:
            app.message("The structure changed while it was being prepared; prepare it again",
                        severity=Severity.WARNING)
            return
        mol, report = result
        doc.set_molecule(mol)
        done_parts = []
        if complete:
            done_parts.append(plural(report.side_chains_completed, "side chain") + " completed")
        if hydrogens:
            done_parts.append(f"hydrogens for pH {ph:.1f}")
        if waters:
            done_parts.append(plural(report.waters_removed, "water molecule") + " removed")
        app.message("Prepared: " + ", ".join(done_parts), severity=Severity.WARNING)

    app.run_in_background(job, doc.molecule.copy(), message="Preparing the structure", done=done)


def register(context: PluginContext) -> None:
    """Add the command before Tools > Show Binding Pocket."""
    context.add_command("Tools > Prepare Protein...", prepare_protein, before="tool choice")
