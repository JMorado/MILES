"""Dialogs of plugins in the colours of the interface (:class:`PluginDialog`)."""

from PySide6.QtWidgets import QDialog, QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from miles.api import App
from miles.ui.controls import dialog_stylesheet

MARGINS = (18, 16, 18, 16)
SPACING = 10


class PluginDialog(QDialog):
    """A dialog in the colours of the interface, with the usual margins and button row.

    Put the contents in :attr:`body` (a vertical layout), add the buttons
    with :meth:`buttons` and show it with :meth:`run`, which deletes the
    dialog once it is closed, so dialogs never pile up in the window.

    Parameters
    ----------
    app : App
        The application, for its window and theme.
    title : str
        Title of the dialog window.

    Examples
    --------
    >>> dialog = PluginDialog(app, "Count Atoms")
    >>> field = QSpinBox()
    >>> dialog.body.addWidget(field)
    >>> dialog.buttons("Count")
    >>> if dialog.run():
    ...     print(field.value())
    """

    def __init__(self, app: App, title: str) -> None:
        super().__init__(app.widget)
        self.setWindowTitle(title)
        self.setStyleSheet(dialog_stylesheet(app.theme))
        self.body = QVBoxLayout(self)
        self.body.setContentsMargins(*MARGINS)
        self.body.setSpacing(SPACING)

    def buttons(self, primary: str = "OK", cancel: str | None = "Cancel",
                extra: QWidget | None = None) -> QPushButton:
        """Add the row of buttons at the bottom: Cancel and the default button, which accepts.

        Parameters
        ----------
        primary : str, optional
            Text of the default button.
        cancel : str or None, optional
            Text of the button that closes the dialog without accepting, or
            None for none.
        extra : QWidget, optional
            Shown at the left of the row, e.g. a check box.

        Returns
        -------
        QPushButton
            The default button, e.g. to grey it out while the input is not valid.
        """
        row = QHBoxLayout()
        if extra is not None:
            row.addWidget(extra)
        row.addStretch(1)
        if cancel is not None:
            rejecting = QPushButton(cancel)
            rejecting.clicked.connect(self.reject)
            row.addWidget(rejecting)
        accepting = QPushButton(primary)
        accepting.setObjectName("Primary")
        accepting.setDefault(True)
        accepting.clicked.connect(self.accept)
        row.addWidget(accepting)
        self.body.addLayout(row)
        return accepting

    def run(self) -> bool:
        """Show the dialog until it is closed; returns whether it was accepted.

        The dialog is deleted afterwards (once control returns to the window),
        so read its fields right after this call.
        """
        accepted = self.exec() == QDialog.DialogCode.Accepted
        self.setParent(None)   # out of the window now, deleted once the caller lets go of it
        self.deleteLater()
        return accepted
