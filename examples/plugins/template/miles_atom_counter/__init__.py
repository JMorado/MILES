"""Atom Counter: how many atoms of each element, a starting point for a MILES plugin.

Adds Tools > Count Atoms, which says how many atoms of each element the
selection has (or the whole structure, when nothing is selected), and
View > Element Counts, which keeps that count in the status line of a tab.
"""

from collections import Counter

from miles.plugin import SEPARATOR, App, DocumentView, PluginContext

NAME = "Atom Counter"
VERSION = "0.1.0"
DESCRIPTION = "Count the atoms of each element in the structure or the selection"
REQUIRES_API = 1   # version of the plugin interface this plugin was written for


def element_counts(doc: DocumentView) -> str:
    """Count the atoms of each element in the selection, or in the whole structure.

    Returns
    -------
    str
        Elements in alphabetical order with their counts, e.g.
        ``"C 2  .  H 6  .  O 1"`` (with MILES's separator dot), or an empty
        string for an empty structure.
    """
    symbols = doc.molecule.symbols
    chosen = doc.selection or range(len(symbols))
    counts = Counter(symbols[i] for i in chosen)
    return SEPARATOR.join(f"{symbol} {n}" for symbol, n in sorted(counts.items()))


def count_atoms(app: App) -> None:
    """Show the counts of the current tab at the top of the window."""
    doc = app.current
    if doc is None or not len(doc.molecule):
        app.message("There are no atoms to count")
        return
    app.message(element_counts(doc))


def register(context: PluginContext) -> None:
    """Add the command, the toggle and the text in the status line."""

    def status_text(doc: DocumentView) -> str | None:
        return element_counts(doc) if toggle.value(doc) else None

    context.add_command("Tools > Count Atoms", count_atoms, enabled=lambda app: app.current is not None)
    status = context.add_status(status_text)
    toggle = context.add_toggle("View > Element Counts", "element_counts", lambda doc, on: status.update(doc))
    context.on_selection_changed(status.update)
