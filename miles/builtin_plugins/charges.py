"""Charges and Dipole: atoms coloured by partial charge, charge labels and the dipole moment.

Adds the View > Charges submenu: Color by Charge, Charge Labels and Dipole
Moment, toggles of each tab, and the charge model (Gasteiger, MMFF94 or
the charges of the file). The dipole is an arrow through the centre of
mass, with its size next to the tip and in the status line.
"""

import numpy as np

from miles import ChargeModel, dipole_moment
from miles.plugin import Annotations, DocumentView, DrawContext, Label, PluginContext, SceneItems, Severity
from miles.plugin.render import LabelMode, arrow_rows, charge_colors

NAME = "Charges and Dipole"
VERSION = "1.0"
DESCRIPTION = "Colour atoms by partial charge, label the charges and draw the dipole moment"

DIPOLE_SCALE = 1.0          # angstrom of arrow per debye (default of the setting)
DIPOLE_MIN = 0.01           # debye; smaller moments are not drawn
DIPOLE_RADIUS = 0.1         # angstrom, shaft of the dipole arrow
DIPOLE_MARGIN = 0.9         # angstrom the arrow reaches beyond the outermost atom, at least


def charge_text(charge: float) -> str:
    """Charge with two decimals and its sign; values that round to zero get no sign."""
    rounded = round(charge, 2)
    return "0.00" if rounded == 0 else f"{rounded:+.2f}"


def register(context: PluginContext) -> None:
    """Add the Charges submenu, the colours, labels and dipole arrow, and the dipole in the status line."""
    app = context.app

    def remember_geometry(doc: DocumentView) -> None:
        """Keep the coordinates the dipole reported by the file belongs to."""
        if doc.info.dipole is not None:
            context.document_cache(doc)["geometry"] = doc.molecule.coords.copy()

    for doc in app.documents:
        remember_geometry(doc)
    context.on_document_opened(remember_geometry)

    def dipole(doc: DocumentView) -> np.ndarray | None:
        """Dipole moment in debye, or None without charges.

        With charges from the file, the dipole the calculation reported is
        used while the geometry is the one it was computed for; otherwise
        the dipole of the partial charges.
        """
        reported, coords = doc.info.dipole, doc.molecule.coords
        geometry = context.document_cache(doc).get("geometry")
        same = (geometry is not None and geometry.shape == coords.shape
                and bool(np.allclose(geometry, coords, atol=1e-4)))
        if doc.charge_model == ChargeModel.FILE and reported is not None and same:
            return np.asarray(reported, dtype=float)
        charges = doc.charges()
        return None if charges is None else dipole_moment(doc.molecule, charges)

    def dipole_arrow(doc: DocumentView) -> tuple[np.ndarray, np.ndarray, float] | None:
        """Tail, tip and size in debye of the dipole arrow.

        The arrow runs through the centre of mass from the negative towards
        the positive end. Both ends reach :data:`DIPOLE_MARGIN` beyond the
        outermost atom; a moment large enough for the arrow's scale (the
        Dipole arrow length setting, angstrom per debye) to make it longer
        than that sets its length.
        """
        mol = doc.molecule
        if not show_dipole.value(doc) or not len(mol):
            return None
        moment = dipole(doc)
        size = float(np.linalg.norm(moment)) if moment is not None else 0.0
        if moment is None or size < DIPOLE_MIN:
            return None
        centre = mol.center_of_mass()
        reach = float(np.linalg.norm(mol.coords - centre, axis=1).max()) + DIPOLE_MARGIN
        half = moment / size * max(size * dipole_scale.value / 2, reach)
        return centre - half, centre + half, size

    def scene(doc: DocumentView, drawing: DrawContext) -> SceneItems | None:
        colors = None
        if color_by_charge.value(doc):
            charges = doc.charges()
            if charges is not None and len(charges) == len(doc.molecule):
                colors = charge_colors(charges)
        arrow = dipole_arrow(doc)
        if colors is None and arrow is None:
            return None
        items = SceneItems(atom_colors=colors)
        if arrow is not None:
            items.arrows = arrow_rows([arrow[0]], [arrow[1]], DIPOLE_RADIUS, drawing.theme.text)
        return items

    def annotations(doc: DocumentView, drawing: DrawContext) -> Annotations | None:
        labels = []
        mol = doc.molecule
        if charge_labels.value(doc):
            charges = doc.charges()
            if charges is not None and len(charges) == len(mol):
                below = 13.0 if doc.display.labels != LabelMode.NONE else 0.0
                labels = [Label(charge_text(float(charges[i])), int(i), size=9.5, offset=(0.0, below))
                          for i in np.nonzero(doc.shown_atoms())[0]]
        arrow = dipole_arrow(doc)
        if arrow is not None:
            labels.append(Label(f"{arrow[2]:.2f} D", arrow[1], away_from=arrow[0]))
        return Annotations(labels=labels) if labels else None

    def status(doc: DocumentView) -> str | None:
        moment = dipole(doc) if show_dipole.value(doc) and len(doc.molecule) else None
        return f"dipole {float(np.linalg.norm(moment)):.2f} D" if moment is not None else None

    drawn = context.add_scene_items(scene)
    written = context.add_annotations(annotations)
    reported = context.add_status(status)

    def rescaled(value: float) -> None:
        drawn.update()
        written.update()

    dipole_scale = context.add_setting("dipole_scale", "Dipole arrow length", DIPOLE_SCALE, minimum=0.1,
                                       maximum=5.0, step=0.1, unit="angstrom per debye", changed=rescaled)

    def refresh(doc: DocumentView) -> None:
        drawn.update(doc)
        written.update(doc)
        reported.update(doc)

    def report_error(doc: DocumentView) -> None:
        """Say why nothing is shown when the chosen model has no charges."""
        error = doc.charge_error()
        if len(doc.molecule) and error:
            app.message(error, severity=Severity.WARNING)

    def switched(doc: DocumentView, on: bool) -> None:
        refresh(doc)
        if on:
            report_error(doc)

    context.add_menu("View > Charges", after=("Surface", "Depth Fade"))
    color_by_charge = context.add_toggle("View > Charges > Color by Charge", "charge_colors", switched)
    charge_labels = context.add_toggle("View > Charges > Charge Labels", "charge_labels", switched)
    show_dipole = context.add_toggle("View > Charges > Dipole Moment", "dipole", switched)
    context.add_separator("View > Charges")

    def model_chosen(doc: DocumentView, model: str) -> None:
        refresh(doc)
        if color_by_charge.value(doc) or charge_labels.value(doc) or show_dipole.value(doc):
            report_error(doc)

    # the model belongs to the document (it is undoable), so the choice reads and writes it there
    models = {str(m): m.label for m in ChargeModel}
    def available(doc: DocumentView, model: str) -> bool:
        """From File only when the file gave charges."""
        return model != ChargeModel.FILE or doc.info.charges is not None

    def write(doc: DocumentView, model: str) -> None:
        doc.set_charge_model(model)

    context.add_choice("View > Charges", "charge_model", models, model_chosen, enabled=available,
                       read=lambda doc: doc.charge_model, write=write)
