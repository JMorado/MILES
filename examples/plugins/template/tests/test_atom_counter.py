"""Tests of the Atom Counter plugin in a MILES window that runs without a screen."""

import miles_atom_counter

from miles.plugin import SEPARATOR, App
from miles.testing import document_from_smiles, load_plugin


def test_counts_the_whole_structure_and_the_selection(miles_app: App) -> None:
    """Every atom is counted, or only the selected ones."""
    load_plugin(miles_app, miles_atom_counter)
    doc = document_from_smiles(miles_app, "CCO", "ethanol")
    assert miles_atom_counter.element_counts(doc) == SEPARATOR.join(["C 2", "H 6", "O 1"])
    doc.select([0, 1])
    assert miles_atom_counter.element_counts(doc) == "C 2"


def test_the_command_and_the_toggle(miles_app: App) -> None:
    """The command runs and the toggle switches without errors of the plugin."""
    info = load_plugin(miles_app, miles_atom_counter)
    doc = document_from_smiles(miles_app, "O", "water")
    miles_app.run("Tools > Count Atoms")
    doc.style(element_counts=True)   # the toggle, as View > Element Counts sets it
    assert miles_app.commands().count("View > Element Counts") == 1
    assert not info.error
