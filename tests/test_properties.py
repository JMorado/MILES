"""Molecular properties, 2D drawings and their commands."""

import time
import xml.etree.ElementTree as ET

import pytest
from PySide6.QtCore import QSettings
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication, QFileDialog

from miles.builtin_plugins.drawing.depiction import depiction_svg
from miles.builtin_plugins.drawing.files import SVG_MIME, svg_image, write_drawing
from miles.builtin_plugins.properties.descriptors import molecular_properties
from miles.builtin_plugins.properties.panel import property_note, property_rows, property_text
from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
from miles.builtin_plugins.shared.molecules import hydrogens_missing, is_small_molecule
from miles.core.chem import ChemError, from_smiles
from miles.ui.mainwindow import MainWindow

SVG = "{http://www.w3.org/2000/svg}"
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"


def mol_of(smiles):
    return from_smiles(smiles)[0]


@pytest.mark.parametrize(("smiles", "formula", "weight", "exact", "inchikey"), [
    (CAFFEINE, "C8H10N4O2", 194.19, 194.08038, "RYYVLZVUVIJVGH-UHFFFAOYSA-N"),
    (ASPIRIN, "C9H8O4", 180.16, 180.04226, "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"),
    ("CCO", "C2H6O", 46.07, 46.04186, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"),
])
def test_composition_and_identifiers(smiles, formula, weight, exact, inchikey):
    props = molecular_properties(mol_of(smiles))
    assert props.formula == formula
    assert props.weight == pytest.approx(weight, abs=0.01)
    assert props.exact_mass == pytest.approx(exact, abs=1e-4)
    assert props.inchikey == inchikey
    assert props.inchi is not None and props.inchi.startswith("InChI=1S/" + formula)
    assert props.charge == 0
    assert props.small_molecule and not props.hydrogens_missing


def test_descriptors_match_reference_values():
    aspirin = molecular_properties(mol_of(ASPIRIN))
    assert aspirin.tpsa == pytest.approx(63.6, abs=0.05)          # Ertl
    assert aspirin.logp == pytest.approx(1.31, abs=0.01)          # Wildman-Crippen
    assert (aspirin.donors, aspirin.acceptors, aspirin.rotatable_bonds) == (1, 4, 3)
    caffeine = molecular_properties(mol_of(CAFFEINE))
    assert caffeine.tpsa == pytest.approx(61.82, abs=0.05)
    assert caffeine.logp == pytest.approx(-1.03, abs=0.01)
    assert caffeine.molar_refractivity == pytest.approx(52.04, abs=0.01)
    assert (caffeine.donors, caffeine.acceptors, caffeine.rotatable_bonds) == (0, 6, 0)
    assert caffeine.lipinski_violations == []
    ethanol = molecular_properties(mol_of("CCO"))
    assert ethanol.tpsa == pytest.approx(20.23, abs=0.05)


def test_stereo_and_charges_reach_the_identifiers():
    alanine = molecular_properties(mol_of("C[C@@H](C(=O)O)N"))   # L-alanine
    assert alanine.inchikey == "QNAYBMKLOCPYGJ-REOHCLBHSA-N"
    assert "@" in (alanine.smiles or "")
    zwitterion = molecular_properties(mol_of("[NH3+]CC(=O)[O-]"))
    assert zwitterion.charge == 0
    assert "+" in (zwitterion.smiles or "") and "-" in (zwitterion.smiles or "")
    acetate = molecular_properties(mol_of("CC(=O)[O-]"))
    assert acetate.charge == -1


def test_rule_of_five_violations():
    long_chain = molecular_properties(mol_of("C" * 30))          # very lipophilic
    assert long_chain.lipinski_violations == ["logP"]


def test_polymers_and_missing_hydrogens_get_composition_only():
    protein = build_sequence("ACDEFGHIKLMNPQRSTVWY", SequenceKind.HELIX)
    assert not is_small_molecule(protein)
    props = molecular_properties(protein)
    assert props.formula.startswith("C") and props.weight > 2000
    assert props.charge is not None and props.smiles is None and props.logp is None
    assert "small molecules" in property_note(props)
    caffeine = mol_of(CAFFEINE)
    heavy = caffeine.subset([i for i, a in enumerate(caffeine.atoms) if a.symbol != "H"])
    assert hydrogens_missing(heavy) and not hydrogens_missing(caffeine)
    bare = molecular_properties(heavy)
    assert bare.hydrogens_missing and bare.charge is None and bare.smiles is None
    assert "Hydrogens" in property_note(bare)


def test_rows_and_text():
    props = molecular_properties(mol_of(ASPIRIN))
    sections = property_rows(props)
    names = [row.name for section in sections for row in section]
    assert names[:2] == ["Formula", "Molecular weight"]
    assert {"SMILES", "InChI", "InChIKey", "logP", "TPSA", "Rule of five"} <= set(names)
    text = property_text(props)
    assert "Formula: C9H8O4" in text and "InChIKey: BSYNRYMUTXBXSQ-UHFFFAOYSA-N" in text
    assert property_note(props) == ""


def test_depiction_is_standalone_svg_and_leaves_the_molecule_alone():
    mol = mol_of("CC(C)C[C@H](N)C(=O)N[C@@H](Cc1ccccc1)C(=O)O")
    before = mol.coords.copy()
    svg = depiction_svg(mol, "dipeptide")
    root = ET.fromstring(svg)
    assert root.tag == f"{SVG}svg" and root.find(f"{SVG}svg") is None   # not nested
    assert float(root.get("width")) > 50 and root.get("viewBox")
    assert root.find(f"{SVG}title").text == "dipeptide"
    texts = "".join(t.text or "" for t in root.iter(f"{SVG}text"))
    assert "N" in texts and "O" in texts
    assert "H3C" not in texts and "CH3" not in texts   # skeletal chain ends
    assert len(list(root.iter(f"{SVG}polygon"))) + len(list(root.iter(f"{SVG}line"))) > 20
    assert (mol.coords == before).all() and len(mol) == len(before)


def test_a_single_carbon_keeps_its_label():
    root = ET.fromstring(depiction_svg(mol_of("C")))
    assert "CH" in "".join(t.text or "" for t in root.iter(f"{SVG}text"))


def test_depiction_refuses_large_structures():
    protein = build_sequence("ACDEFGHIKLMNPQRSTVWY" * 3, SequenceKind.HELIX)
    with pytest.raises(ChemError, match="select a ligand"):
        depiction_svg(protein)


def test_drawing_files(app, tmp_path):
    svg = depiction_svg(mol_of(CAFFEINE))
    write_drawing(svg, str(tmp_path / "caffeine.svg"))
    assert ET.parse(tmp_path / "caffeine.svg").getroot().tag == f"{SVG}svg"
    write_drawing(svg, str(tmp_path / "caffeine.png"))
    image = QImage(str(tmp_path / "caffeine.png"))
    assert image.width() > 300 and image.height() > 200
    assert not svg_image(svg).isNull()


@pytest.fixture
def window(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    win.resize(1000, 700)
    win.show()
    app.processEvents()
    yield win
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def wait_until(condition, timeout=5.0):
    end = time.time() + timeout
    while not condition():
        assert time.time() < end, "timed out"
        QApplication.processEvents()
        time.sleep(0.01)


def test_properties_panel_follows_the_selection_and_edits(window):
    window.open_structure(mol_of(ASPIRIN), window.current_doc().info)
    canvas = window.current_canvas()
    panel = window.stage.floating_panels()["properties"]
    window.api.run("Properties")
    assert panel.isVisible() and panel.props is not None
    assert panel.props.inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    panel.copy_all.click()
    assert "TPSA: 63.6" in QApplication.clipboard().text()

    oxygens = [i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "O"]
    canvas.set_selection(oxygens)
    wait_until(lambda: panel.props is not None and panel.props.formula == "HO4")   # with the OH hydrogen
    assert panel.subtitle.text().startswith("Selection") and "with 1 bonded hydrogen" in panel.subtitle.text()

    canvas.set_selection([])
    wait_until(lambda: panel.props is not None and panel.props.formula == "C9H8O4")
    hydrogen = next(i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "H")
    canvas.doc.edit(lambda mol: mol.remove_atoms([hydrogen]))
    wait_until(lambda: panel.props is not None and panel.props.formula == "C9H7O4")

    first = canvas
    window.new_document()   # another tab: the panel follows it, and hides for an empty one
    assert not panel.isVisible() and panel.wanted
    window.api.run("Properties")   # nothing to show here: Cmd I says so and the panel stays wanted
    assert not panel.isVisible() and panel.wanted
    assert window.current_canvas().toast_text == panel.unavailable
    window.tabbar.setCurrentIndex(window.canvases.index(first))
    assert panel.isVisible()
    window.api.run("Properties")   # Cmd I again closes it
    assert not panel.isVisible() and not panel.wanted


def test_drawing_commands(window, tmp_path, monkeypatch):
    window.open_structure(mol_of(CAFFEINE), window.current_doc().info)
    target = tmp_path / "drawing.svg"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args, **kwargs: (str(target), ""))
    window.api.run("File > Export 2D Drawing")
    assert ET.parse(target).getroot().tag == f"{SVG}svg"
    window.api.run("Edit > Copy 2D Drawing")
    mime = QApplication.clipboard().mimeData()
    assert mime.hasImage() and mime.hasFormat(SVG_MIME)
