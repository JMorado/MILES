"""Molecular Surface: the solvent-excluded surface of the selection or the whole structure.

Adds the View > Surface submenu: Show Surface, Hide Surface, the colouring
(element, hydrophobicity, charge, electrostatic potential or one colour)
and See-Through. The electrostatic potential takes the partial charges of
the charge model chosen in View > Charges (the charges of ionisable groups
for proteins and nucleic acids). The surface is worked out on a worker
thread, so large structures do not stop the window; it follows edits and
trajectory frames, and exported pictures and movie frames get the exact
surface of what they show.
"""

from dataclasses import dataclass

import numpy as np

from miles import ChargeModel, Molecule, has_polymer
from miles.plugin import (
    SEPARATOR,
    App,
    Change,
    DocumentView,
    DrawContext,
    MeshItem,
    PluginContext,
    SceneItems,
    Severity,
    Theme,
)
from miles.plugin.render import Mesh, SurfaceColoring, SurfaceRequest, ionisable_charges, surface_mesh

NAME = "Molecular Surface"
VERSION = "1.0"
DESCRIPTION = "Show the molecular surface, coloured by element, hydrophobicity, charge or potential"

COLORINGS = {SurfaceColoring.ELEMENT: "Element Colors", SurfaceColoring.HYDROPHOBICITY: "Hydrophobicity",
             SurfaceColoring.CHARGE: "Charge", SurfaceColoring.POTENTIAL: "Electrostatic Potential",
             SurfaceColoring.UNIFORM: "One Color"}
"""The colourings with their titles, in the order of the menu."""


@dataclass(frozen=True)
class Shown:
    """The surface a tab shows: the atoms it wraps, by atom id (None: everything drawn)."""

    uids: tuple[int, ...] | None


def surface_charges(doc: DocumentView) -> tuple[np.ndarray | None, str]:
    """Charges for the electrostatic potential on a tab's surface, and where they come from.

    Proteins and nucleic acids use the charges of their ionisable groups,
    which need no hydrogens and no protonation model, unless the file gave
    charges and they are chosen; other structures use the chosen charge
    model.

    Returns
    -------
    tuple of (numpy.ndarray or None, str)
        The charges, or None with the reason there are none.
    """
    from_file = doc.charge_model == "file" and doc.file_charges_apply()
    if has_polymer(doc.molecule) and not from_file:
        return ionisable_charges(doc.molecule), "charges of ionisable groups"
    charges = doc.charges()
    if charges is None:
        return None, doc.charge_error()
    return charges, f"{ChargeModel(doc.charge_model).label} charges"


def compute(mol: Molecule, request: SurfaceRequest, theme: Theme, charges: np.ndarray | None,
            display: object) -> Mesh:
    """Work a surface out (on a worker thread: everything it gets is a copy)."""
    return surface_mesh(mol, request, theme, charges, display)   # type: ignore[arg-type]


def register(context: PluginContext) -> None:
    """Add the Surface submenu and draw the surface of each tab."""
    app = context.app

    def shown(doc: DocumentView) -> Shown | None:
        found: Shown | None = context.document_cache(doc).get("shown")
        return found

    def request(doc: DocumentView) -> SurfaceRequest | None:
        """Return what the tab's surface wraps, with the atoms as they are now, and how it is coloured."""
        wanted = shown(doc)
        if wanted is None:
            return None
        atoms = None
        if wanted.uids is not None:
            atoms = tuple(np.flatnonzero(np.isin(doc.molecule.uids, wanted.uids)).tolist())
            if not atoms:   # every atom it wrapped was deleted
                return None
        return SurfaceRequest(atoms, SurfaceColoring(coloring.value(doc)), see_through.value(doc))

    def tell_range(doc: DocumentView, mesh: Mesh | None) -> None:
        """Say what the colours of a potential surface that was asked for mean, once it is there."""
        source = context.document_cache(doc).pop("explain", None)
        if source and mesh is not None and not mesh.empty:
            r = mesh.value_range
            app.message(f"Electrostatic potential from {source}{SEPARATOR}"
                        f"red -{r:.1f} to blue +{r:.1f} kcal/(mol e)", severity=Severity.WARNING)

    def scene(doc: DocumentView, drawing: DrawContext) -> SceneItems | None:
        wanted = request(doc)
        if wanted is None:
            return None
        if drawing.dragging and not drawing.export:   # atoms are being dragged: keep what is there
            kept = surfaces.latest(doc)
            return SceneItems(meshes=[MeshItem(kept, wanted.alpha)]) if kept is not None else None
        potential = wanted.coloring == SurfaceColoring.POTENTIAL
        charges, _ = surface_charges(doc) if potential else (None, "")
        key = (doc.revision, wanted.atoms, wanted.coloring, drawing.theme.name,
               (doc.charge_model, doc.file_charges_apply()) if potential else None)
        mesh = surfaces.get(doc, key, compute, doc.molecule.copy(), wanted, drawing.theme,
                            None if charges is None else np.array(charges), doc.atom_display,
                            now=drawing.export)
        if not drawing.export and surfaces.has(doc, key):
            tell_range(doc, mesh)
        return SceneItems(meshes=[MeshItem(mesh, wanted.alpha)]) if mesh is not None else None

    layer = context.add_scene_items(scene)
    surfaces = context.background_results(layer, "Computing the surface")

    def explain(doc: DocumentView) -> None:
        """Say where a potential surface's charges come from, and what its colours mean once it arrives."""
        wanted = request(doc)
        if wanted is None or wanted.coloring != SurfaceColoring.POTENTIAL:
            return
        charges, source = surface_charges(doc)
        if charges is None:
            app.message(f"No electrostatic potential: {source}", severity=Severity.WARNING)
            return
        context.document_cache(doc)["explain"] = source
        layer.update(doc)

    def current() -> DocumentView | None:
        doc = app.current
        return doc if doc is not None and len(doc.molecule) else None

    def show_surface(app: App) -> None:
        doc = current()
        if doc is None:
            return
        uids = tuple(doc.molecule.uids[sorted(doc.selection)].tolist()) if doc.selection else None
        context.document_cache(doc)["shown"] = Shown(uids)
        layer.update(doc)
        explain(doc)

    def hide_surface(app: App) -> None:
        doc = current()
        if doc is not None:
            context.document_cache(doc).pop("shown", None)
            layer.update(doc)

    def colored(doc: DocumentView, value: str) -> None:
        layer.update(doc)
        explain(doc)

    def charges_changed(doc: DocumentView, change: str) -> None:
        """Draw a potential surface again for a new charge model, and say what its colours mean now."""
        if change == Change.PROPERTIES and shown(doc) is not None:
            wanted = request(doc)
            if wanted is not None and wanted.coloring == SurfaceColoring.POTENTIAL:
                layer.update(doc)
                explain(doc)

    context.add_menu("View > Surface", after="Depth Fade")
    context.add_command("View > Surface > Show Surface", show_surface)
    context.add_command("View > Surface > Hide Surface", hide_surface)
    context.add_separator("View > Surface")
    options = {str(k): v for k, v in COLORINGS.items()}
    coloring = context.add_choice("View > Surface", "surface_coloring", options, colored,
                                  default=SurfaceColoring.ELEMENT)
    context.add_separator("View > Surface")
    see_through = context.add_toggle("View > Surface > See-Through", "surface_see_through",
                                     lambda doc, on: layer.update(doc))
    context.on_document_changed(charges_changed)
