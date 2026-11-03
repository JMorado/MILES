"""Superposing the structure of one tab onto another."""

from PySide6.QtWidgets import QCheckBox, QComboBox, QLabel

from miles import Molecule
from miles.builtin_plugins.superpose.matching import Superposition, apply, superpose
from miles.plugin import ANGSTROM, App, Change, DocumentView, PluginDialog
from miles.plugin.widgets import theme_palette


class SuperposeDialog(PluginDialog):
    """Choose the structure to superpose the current one onto."""

    def __init__(self, app: App, mobile: DocumentView, targets: list[DocumentView]) -> None:
        super().__init__(app, "Superpose")
        self.mobile = mobile
        self.targets = targets
        self.body.addWidget(QLabel(f"Move {mobile.name} onto:"))
        self.target = QComboBox()
        self.target.addItems([doc.name for doc in targets])
        theme_palette(self.target, app.theme)   # readable when the app and system themes differ
        self.body.addWidget(self.target)
        self.copy = QCheckBox("Also add a copy to that tab, to compare them")
        self.copy.setChecked(True)
        self.body.addWidget(self.copy)
        self.buttons("Superpose")

    @property
    def chosen(self) -> DocumentView:
        """The target document."""
        return self.targets[self.target.currentIndex()]


def superpose_documents(mobile: DocumentView, target: DocumentView, add_copy: bool) -> Superposition:
    """Move one document's structure onto another's, as an undoable edit.

    Parameters
    ----------
    mobile : DocumentView
        Document whose molecule moves.
    target : DocumentView
        Document that stays in place.
    add_copy : bool
        Also append a copy of the moved structure to the target document.

    Returns
    -------
    Superposition
        The transformation and the RMSD.

    Raises
    ------
    ValueError
        If the two structures cannot be matched.
    """
    result = superpose(mobile.molecule, target.molecule)
    mobile.edit(apply, result, kind=Change.COORDS)
    if add_copy:
        target.edit(Molecule.extend, mobile.molecule.copy())
    return result


def describe(result: Superposition) -> str:
    """One line summarising a superposition."""
    return f"RMSD {result.rmsd:.2f} {ANGSTROM} over {len(result.pairs)} atoms ({result.method})"
