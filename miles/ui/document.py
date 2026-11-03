"""The window's documents: :class:`miles.core.document.Document` with Qt signals.

The editing, undo and trajectory logic lives in :mod:`miles.core.document`,
which needs no Qt; this class announces its changes as signals for the
widgets.
"""

from PySide6.QtCore import QObject, Signal

# Change, ChangeEvent, TOPOLOGY and UNDO_LIMIT are documented in miles.core.document;
# code in the window imports them from here
from miles.core.document import TOPOLOGY as TOPOLOGY
from miles.core.document import UNDO_LIMIT as UNDO_LIMIT
from miles.core.document import Change as Change
from miles.core.document import ChangeEvent as ChangeEvent
from miles.core.document import Document as CoreDocument
from miles.core.jobinfo import JobInfo
from miles.core.molecule import Molecule

__all__ = ["Document"]


class Document(QObject, CoreDocument):
    """A document whose changes are also Qt signals.

    Signals
    -------
    changed(str)
        A :class:`Change` value after every change.
    changed_detail(ChangeEvent)
        The same change with the atoms and properties it concerns.
    renamed()
        The name or the modified flag changed.
    """

    changed = Signal(str)
    changed_detail = Signal(object)
    renamed = Signal()

    def __init__(self, mol: Molecule | None = None, info: JobInfo | None = None,
                 path: str | None = None) -> None:
        QObject.__init__(self)
        CoreDocument.__init__(self, mol, info, path)

    def _announce(self, event: ChangeEvent) -> None:
        CoreDocument._announce(self, event)
        self.changed.emit(str(event.kind))
        self.changed_detail.emit(event)

    def _renamed(self) -> None:
        self.renamed.emit()
