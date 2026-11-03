"""Building peptides, DNA and RNA from a typed or pasted sequence."""

from typing import override

from PySide6.QtWidgets import QPlainTextEdit

from miles.builtin_plugins.sequence.biopolymer import SequenceKind, clean_sequence
from miles.plugin import App, PluginDialog
from miles.plugin.widgets import Segmented, label

KINDS = [(SequenceKind.HELIX, "Protein Helix"), (SequenceKind.STRAND, "Protein Strand"),
         (SequenceKind.DNA, "DNA"), (SequenceKind.DNA_SINGLE, "Single DNA"), (SequenceKind.RNA, "RNA")]
EXAMPLES = {SequenceKind.HELIX: "AEAAAKEAAAKA", SequenceKind.STRAND: "TWTWTE",
            SequenceKind.DNA: "CGCGAATTCGCG", SequenceKind.DNA_SINGLE: "ACGTACGT",
            SequenceKind.RNA: "GGACUUCGGUCC"}


class SequenceDialog(PluginDialog):
    """Sequence field and polymer type; after :meth:`run` accepts, :attr:`sequence_text` is checked.

    The structure itself is built after the dialog, in the background.
    """

    def __init__(self, app: App) -> None:
        super().__init__(app, "Build from Sequence")
        self.setStyleSheet(self.styleSheet()
                           + "QPlainTextEdit { padding: 6px; font-family: Menlo, Consolas, monospace; }")
        self.sequence_text = ""
        self.title = ""
        self.kind = Segmented([name for _, name in KINDS])
        self.kind.changed.connect(self._kind_changed)
        self.body.addWidget(self.kind)
        self.sequence = QPlainTextEdit()
        self.sequence.setFixedSize(520, 120)
        self.body.addWidget(self.sequence)
        self.message = label("One-letter codes; spaces, numbers and FASTA headers are ignored.", "Muted")
        self.message.setWordWrap(True)
        self.body.addWidget(self.message)
        self.buttons("Build")
        self._kind_changed(0)

    def chosen_kind(self) -> SequenceKind:
        """Return the kind of polymer chosen."""
        return KINDS[self.kind.index()][0]

    def _kind_changed(self, index: int) -> None:
        self.sequence.setPlaceholderText(f"e.g. {EXAMPLES[self.chosen_kind()]}")

    @override
    def accept(self) -> None:
        """Close only with a sequence that can be built; otherwise say what is wrong."""
        try:
            sequence = clean_sequence(self.sequence.toPlainText(), self.chosen_kind())
        except ValueError as exc:
            self.message.setText(str(exc))
            return
        short = sequence if len(sequence) <= 12 else sequence[:10] + "..."
        self.sequence_text = sequence
        self.title = f"{KINDS[self.kind.index()][1]} {short}"
        super().accept()
