"""Binding Pocket: the residues around a ligand, with its hydrogen bonds, salt bridges and pi stacking.

Adds Tools > Show Binding Pocket (for the selected ligand, or the first
one) and the same item in the right-click menu of a ligand atom. The
pocket is drawn as tubes and the ligand as balls; dashed lines with their
lengths show the interactions until the structure changes or Clear
Interactions in the right-click menu removes them.
"""

import numpy as np

from miles import AMINO_ACIDS, Representation, has_polymer, is_polymer, is_water
from miles.builtin_plugins.binding_pocket.interactions import (
    Interaction,
    InteractionKind,
    find_interactions,
    ligand_atoms,
    pocket_atoms,
)
from miles.builtin_plugins.shared.hbonds import HBOND_ANGLE, HBOND_MAX
from miles.plugin import (
    ANGSTROM,
    Annotations,
    App,
    AtomDisplay,
    Change,
    DocumentView,
    DrawContext,
    Label,
    Line,
    MenuItem,
    MenuTarget,
    PluginContext,
    Severity,
    plural,
)

NAME = "Binding Pocket"
VERSION = "1.0"
DESCRIPTION = "Show the pocket around a ligand and its interactions"

HYDROGEN_BONDS = "miles.builtin_plugins.hydrogen_bonds"
"""The plugin whose distance and angle settings decide what a hydrogen bond is."""

COLORS = {
    InteractionKind.HYDROGEN_BOND: "#3A9BDC",
    InteractionKind.SALT_BRIDGE: "#E0703A",
    InteractionKind.PI_STACKING: "#33A86A",
}
"""Colour of each kind of interaction."""


def register(context: PluginContext) -> None:
    """Add the command, the right-click items and the interaction lines."""
    app = context.app

    def found(doc: DocumentView) -> list[Interaction]:
        """Return the interactions shown in a tab (none when the list is empty)."""
        return context.document_cache(doc).setdefault("interactions", [])

    def show_pocket_of(doc: DocumentView, atom: int) -> None:
        """Show the binding pocket of the ligand containing an atom, with its interactions."""
        mol = doc.molecule
        ligand = ligand_atoms(mol, atom)
        pocket = pocket_atoms(mol, ligand, radius.value)
        if not pocket:
            app.message(f"No residues within {radius.value:g} angstrom of this molecule",
                        severity=Severity.WARNING)
            return
        distance = context.setting_value(HYDROGEN_BONDS, "max_distance", HBOND_MAX)
        angle = context.setting_value(HYDROGEN_BONDS, "min_angle", HBOND_ANGLE)
        interactions = find_interactions(mol, ligand, pocket, max_distance=distance, min_angle=angle)
        found(doc)[:] = interactions

        def apply(display: AtomDisplay) -> AtomDisplay:
            hidden = display.hidden.copy()
            hidden[pocket + ligand] = False
            return (display.with_rep(Representation.TUBE, pocket).with_rep(Representation.BALL, ligand)
                    .with_hidden(hidden))

        doc.edit_display(apply)
        doc.fit(ligand + pocket)
        residues = len({(mol.atoms[i].chain, mol.atoms[i].resseq) for i in pocket})
        counts = [plural(sum(x.kind == kind for x in interactions), str(kind)) for kind in
                  dict.fromkeys(x.kind for x in interactions)]
        app.message(f"Pocket of {plural(residues, 'residue')}" + "".join(f", {c}" for c in counts))
        lines.update(doc)

    def show_pocket(app: App) -> None:
        doc = app.current
        if doc is None or not len(doc.molecule):
            return
        # a metal ion bound to the protein has a pocket too; without a polymer there is no pocket at all
        ligands = doc.find("ligand or ion") if has_polymer(doc.molecule) else []
        chosen = [i for i in doc.selection if i in set(ligands)] or ligands
        if not chosen:
            app.message("This structure has no ligand", severity=Severity.WARNING)
            return
        show_pocket_of(doc, chosen[0])

    def clear(doc: DocumentView) -> None:
        if found(doc):
            found(doc).clear()
            lines.update(doc)

    def draw(doc: DocumentView, drawing: DrawContext) -> Annotations | None:
        interactions = found(doc)
        mol = doc.molecule
        if not interactions:
            return None
        result = Annotations()
        for interaction in interactions:
            if max(interaction.first + interaction.second) >= len(mol):
                continue
            a = mol.coords[list(interaction.first)].mean(axis=0)
            b = mol.coords[list(interaction.second)].mean(axis=0)
            color = COLORS.get(interaction.kind, "#888888")
            result.lines.append(Line(a, b, color, 1.8, (4.0, 2.0)))
            length = float(np.linalg.norm(a - b))
            result.labels.append(Label(f"{length:.1f} {ANGSTROM}", (a + b) / 2, color, size=9.5,
                                       offset=(0.0, -9.0), halo=0.9))
        return result

    def atom_items(doc: DocumentView, target: MenuTarget) -> list[MenuItem]:
        mol = doc.molecule
        i = target.atom
        if i is None or i >= len(mol):
            return []
        atom = mol.atoms[i]
        if atom.resname in AMINO_ACIDS or not has_polymer(mol) or is_polymer(atom) or is_water(atom):
            return []
        return [MenuItem("Show Binding Pocket", lambda: show_pocket_of(doc, i))]

    def selection_items(doc: DocumentView, target: MenuTarget) -> list[MenuItem]:
        if not found(doc):
            return []
        return [MenuItem("Clear Interactions", lambda: clear(doc))]

    radius = context.add_setting("radius", "Pocket radius", 5.0, minimum=3.0, maximum=10.0, step=0.5,
                                 unit="angstrom")
    lines = context.add_annotations(draw)
    context.on_document_changed(lambda doc, change: clear(doc) if change == Change.STRUCTURE else None)
    context.add_context_items(atom_items, "atom")
    context.add_context_items(selection_items, "selection")
    context.add_command("Tools > Show Binding Pocket", show_pocket, before="tool choice")
