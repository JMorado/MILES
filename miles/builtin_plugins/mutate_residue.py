"""Mutate Residue: replace the side chain of an amino acid from its right-click menu.

Right-clicking an atom of an amino acid offers Mutate <residue> with the
other amino acids; the new side chain is built in place, as one undoable
step, and selected.
"""

import functools

from miles import AMINO_ACIDS, amino_acid_codes
from miles.builtin_plugins.shared.mutation import mutate_residue
from miles.plugin import App, DocumentView, MenuItem, MenuTarget, PluginContext, Severity

NAME = "Mutate Residue"
VERSION = "1.0"
DESCRIPTION = "Replace the side chain of an amino acid from the right-click menu"


def mutate(app: App, doc: DocumentView, atom: int, code: str) -> None:
    """Replace the side chain of the residue containing an atom, and select the new atoms."""
    try:
        new = doc.edit(mutate_residue, atom, code)
    except ValueError as exc:
        app.message(str(exc), severity=Severity.WARNING)
        return
    doc.select(new)


def register(context: PluginContext) -> None:
    """Add the Mutate submenu to the right-click menu of amino acid atoms."""
    app = context.app

    def items(doc: DocumentView, target: MenuTarget) -> list[MenuItem]:
        mol = doc.molecule
        i = target.atom
        if i is None or i >= len(mol) or mol.atoms[i].resname not in AMINO_ACIDS:
            return []
        atom = mol.atoms[i]
        choices = [MenuItem(code, functools.partial(mutate, app, doc, i, code))
                   for code in amino_acid_codes() if code != atom.resname]
        return [MenuItem(f"Mutate {atom.resname} {atom.resseq}", items=choices)]

    context.add_context_items(items, "atom")
