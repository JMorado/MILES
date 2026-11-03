"""Preparing a protein: missing side chains, hydrogens for a pH, water."""

from PySide6.QtWidgets import QCheckBox, QDoubleSpinBox, QGridLayout

from miles import Molecule
from miles.builtin_plugins.prepare_protein.preparation import incomplete_residues, water_count
from miles.plugin import App, PluginDialog, plural
from miles.plugin.widgets import label


class PrepareDialog(PluginDialog):
    """Options for :func:`miles.builtin_plugins.prepare_protein.preparation.prepare`."""

    def __init__(self, app: App, mol: Molecule) -> None:
        super().__init__(app, "Prepare Protein")
        missing = len(incomplete_residues(mol))
        waters = water_count(mol)
        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(8)
        self.complete = QCheckBox("Complete missing side chains")
        self.complete.setChecked(missing > 0)
        self.complete.setEnabled(missing > 0)
        grid.addWidget(self.complete, 0, 0)
        grid.addWidget(label(f"{plural(missing, 'residue')} incomplete", "Muted"), 0, 1)
        self.hydrogens = QCheckBox("Add hydrogens for pH")
        self.hydrogens.setChecked(True)
        grid.addWidget(self.hydrogens, 1, 0)
        self.ph = QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setSingleStep(0.1)
        self.ph.setDecimals(1)
        self.ph.setValue(7.4)
        grid.addWidget(self.ph, 1, 1)
        self.water = QCheckBox("Remove water")
        self.water.setEnabled(waters > 0)
        grid.addWidget(self.water, 2, 0)
        grid.addWidget(label(plural(waters, "water molecule"), "Muted"), 2, 1)
        self.body.addLayout(grid)
        self.buttons("Prepare")
