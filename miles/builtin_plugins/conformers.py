"""Generate Conformers: low-energy conformers of a small molecule, as frames of the timeline.

Adds Tools > Generate Conformers.
"""

import numpy as np

import miles
from miles.plugin import App, PluginContext, Severity, plural

NAME = "Generate Conformers"
VERSION = "1.0"
DESCRIPTION = "Search for low-energy conformers of a small molecule"

ATOM_LIMIT = 200   # the conformer search is meant for small molecules


def generate_conformers(app: App) -> None:
    """Search for conformers in the background and show them as frames, as one undoable step."""
    doc = app.current
    if doc is None:
        return
    mol = doc.molecule
    if any(a.resname for a in mol.atoms) or len(mol) > ATOM_LIMIT:
        app.message(f"The conformer search is for small molecules (up to {ATOM_LIMIT} atoms)",
                    severity=Severity.WARNING)
        return
    if len(mol) < 3:
        return
    revision = doc.revision

    def done(found: list[tuple[float, np.ndarray]]) -> None:
        if doc not in app.documents:
            return
        if doc.revision != revision:
            app.message("The structure changed during the conformer search; search again",
                        severity=Severity.WARNING)
            return
        doc.set_frames(np.array([xyz for _, xyz in found]), energies=[energy for energy, _ in found])
        app.message(f"Found {plural(len(found), 'conformer')}; step through them with the timeline")

    app.run_in_background(miles.conformers, mol.copy(), message="Searching for conformers", done=done,
                          failed=lambda text: app.message(text, severity=Severity.WARNING))


def register(context: PluginContext) -> None:
    """Add the command after Tools > Build from Sequence..."""
    context.add_command("Tools > Generate Conformers", generate_conformers,
                        before=("Vibrations", "protein tools"))
