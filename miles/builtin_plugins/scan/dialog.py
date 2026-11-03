"""Relaxed scan of a distance, angle or dihedral: range, number of steps and the run."""

from PySide6.QtWidgets import QDoubleSpinBox, QGridLayout, QSpinBox

from miles import Molecule
from miles.builtin_plugins.scan.relaxed_scan import Coordinate, coordinate_kind, measure
from miles.plugin import ANGSTROM, DEGREE, App, PluginDialog, atom_name
from miles.plugin.widgets import label

# range, decimals, single step and default (offset of start, offset of stop, steps) per coordinate
SETTINGS = {
    Coordinate.DISTANCE: ((0.5, 20.0), 3, 0.05, (0.0, 1.0, 10)),
    Coordinate.ANGLE: ((1.0, 179.0), 1, 5.0, (-20.0, 20.0, 8)),
    Coordinate.DIHEDRAL: ((-720.0, 720.0), 1, 10.0, (0.0, 360.0, 36)),
}


def value_text(kind: Coordinate, value: float) -> str:
    """Format a coordinate value with its unit: angstrom for distances, degrees otherwise."""
    if kind == Coordinate.DISTANCE:
        return f"{value:.3f} {ANGSTROM}"
    return f"{value:.1f}{DEGREE}"


class ScanDialog(PluginDialog):
    """Start, stop and steps of a scan; the scan itself runs after the dialog, in the background.

    Parameters
    ----------
    app : App
        The application.
    mol : Molecule
        Structure to scan (only measured).
    atoms : list of int
        Two, three or four atoms, in order.
    """

    def __init__(self, app: App, mol: Molecule, atoms: list[int]) -> None:
        super().__init__(app, "Scan Coordinate")
        self.setStyleSheet(self.styleSheet() + "QAbstractSpinBox { min-width: 90px; }")
        self.atoms = list(atoms)
        self.kind = coordinate_kind(atoms)
        current = measure(mol.coords, atoms)
        (low, high), decimals, single, (before, after, steps) = SETTINGS[self.kind]
        names = "-".join(atom_name(mol, i) for i in atoms)
        now = value_text(self.kind, current)
        self.body.addWidget(label(f"{self.kind.value.capitalize()} {names}, now {now}"))
        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(8)
        suffix = f" {ANGSTROM}" if self.kind == Coordinate.DISTANCE else DEGREE
        self.start = self._field(low, high, decimals, single, suffix, current + before)
        self.stop = self._field(low, high, decimals, single, suffix, current + after)
        self.steps = QSpinBox()
        self.steps.setRange(1, 360)
        self.steps.setValue(steps)
        rows = [("From", self.start), ("To", self.stop), ("Steps", self.steps)]
        for row, (text, widget) in enumerate(rows):
            grid.addWidget(label(text), row, 0)
            grid.addWidget(widget, row, 1)
        self.body.addLayout(grid)
        note = label("Each point sets the coordinate, then relaxes the rest of the structure with "
                     "MMFF94 (UFF when MMFF94 has no parameters). The scan runs in the background.", "Muted")
        note.setWordWrap(True)
        note.setMinimumWidth(360)
        self.body.addWidget(note)
        self.buttons("Scan")

    @staticmethod
    def _field(low: float, high: float, decimals: int, single: float, suffix: str,
               value: float) -> QDoubleSpinBox:
        field = QDoubleSpinBox()
        field.setRange(low, high)
        field.setDecimals(decimals)
        field.setSingleStep(single)
        field.setSuffix(suffix)
        field.setValue(min(max(value, low), high))
        return field
