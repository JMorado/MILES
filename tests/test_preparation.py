import time

import numpy as np
import pytest

from miles.builtin_plugins.prepare_protein.preparation import (
    add_hydrogens_at_ph,
    complete_side_chains,
    incomplete_residues,
    prepare,
    remove_water,
)
from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
from miles.core.document import Document
from miles.core.molecule import Atom, Representation
from miles.core.neighbors import pairs_within
from miles.files.formats import parse_text

WATER = (
    "HETATM  100  O   HOH A 101      20.000   0.000   0.000  1.00  0.00           O\n"
    "HETATM  101  O   HOH A 102      24.000   0.000   0.000  1.00  0.00           O\n"
)


def experimental(sequence, kind=SequenceKind.STRAND):
    """A built peptide without hydrogens, as it would come from a PDB file."""
    mol = build_sequence(sequence, kind)
    mol.remove_atoms([i for i, a in enumerate(mol.atoms) if a.symbol == "H"])
    return mol


def hydrogens_on(mol, resseq, *names):
    count = 0
    for i, a in enumerate(mol.atoms):
        if a.resseq == resseq and a.name in names:
            count += sum(mol.atoms[j].symbol == "H" for j in mol.neighbors(i))
    return count


def truncated():
    """GKDAEGR with the ends of the Lys and Glu side chains missing."""
    mol = experimental("GKDAEGR")
    missing = {("LYS", "CE"), ("LYS", "NZ"), ("GLU", "OE1"), ("GLU", "OE2")}
    mol.remove_atoms([i for i, a in enumerate(mol.atoms) if (a.resname, a.name) in missing])
    return mol


def test_incomplete_side_chains_are_found_and_rebuilt():
    mol = truncated()
    found = incomplete_residues(mol)
    assert [(mol.atoms[i].resname, mol.atoms[i].name) for i in found] == [("LYS", "CA"), ("GLU", "CA")]
    assert complete_side_chains(mol) == 2
    assert incomplete_residues(mol) == []
    names = {a.name for a in mol.atoms if a.resname == "LYS"}
    assert {"CE", "NZ"} <= names


@pytest.mark.parametrize(("ph", "n_term", "lys", "asp", "glu", "c_term"), [
    (7.4, 3, 3, 0, 0, 0),
    (1.0, 3, 3, 1, 1, 1),
    (13.0, 2, 2, 0, 0, 0),
])
def test_protonation_depends_on_ph(ph, n_term, lys, asp, glu, c_term):
    mol = truncated()
    complete_side_chains(mol)
    added = add_hydrogens_at_ph(mol, ph)
    assert added == sum(a.symbol == "H" for a in mol.atoms)
    assert hydrogens_on(mol, 1, "N") == n_term
    assert hydrogens_on(mol, 2, "NZ") == lys
    assert hydrogens_on(mol, 3, "OD1", "OD2") == asp
    assert hydrogens_on(mol, 5, "OE1", "OE2") == glu
    assert hydrogens_on(mol, 7, "O", "OXT") == c_term


def test_hydrogens_are_placed_sensibly_and_named():
    mol = truncated()
    doc = Document(mol)
    doc.set_representation(Representation.TUBE, [5])
    uid = mol.uids[5]
    doc.edit(prepare)
    assert mol.uids[5] == uid and doc.display.rep[5] == Representation.TUBE   # the same atom, style kept
    hydrogens = [i for i, a in enumerate(mol.atoms) if a.symbol == "H"]
    lengths = [np.linalg.norm(mol.coords[h] - mol.coords[mol.neighbors(h)[0]]) for h in hydrogens]
    assert min(lengths) > 0.9 and max(lengths) < 1.15
    assert len(pairs_within(mol.coords, 0.8)) == 0
    lysine = sorted(a.name for a in mol.atoms if a.resname == "LYS" and a.symbol == "H")
    assert lysine == ["H", "HA", "HB2", "HB3", "HD2", "HD3", "HE2", "HE3", "HG2", "HG3", "HZ1", "HZ2", "HZ3"]
    assert all(a.resname == "LYS" for a in mol.atoms if a.name == "HZ1")
    assert sorted(a.name for a in mol.atoms if a.resseq == 1 and a.name.startswith("H")) == [
        "H1", "H2", "H3", "HA2", "HA3"]


def test_existing_hydrogens_are_replaced():
    mol = build_sequence("GKG", SequenceKind.STRAND)
    before = sum(a.symbol == "H" for a in mol.atoms)
    add_hydrogens_at_ph(mol)
    after = sum(a.symbol == "H" for a in mol.atoms)
    assert after == before + 1   # NH3+ at the N-terminus and on Lys, COO- at the C-terminus
    assert hydrogens_on(mol, 2, "NZ") == 3


def test_water_is_removed():
    mol = parse_text(_pdb(experimental("GAG")) + WATER)
    assert remove_water(mol) == 2
    assert not any(a.resname == "HOH" for a in mol.atoms)
    report = prepare(parse_text(_pdb(experimental("GAG")) + WATER), waters=True, hydrogens=False)
    assert report.waters_removed == 2 and report.hydrogens_added == 0


def test_preparing_a_large_protein_is_fast():
    mol = experimental("ACDEFGHIKLMNPQRSTVWY" * 10, SequenceKind.HELIX)
    start = time.perf_counter()
    added = add_hydrogens_at_ph(mol)
    assert added > 1500
    assert time.perf_counter() - start < 5.0


def _pdb(mol):
    lines = []
    for k, (a, (x, y, z)) in enumerate(zip(mol.atoms, mol.coords, strict=True), start=1):
        name = a.name if len(a.name) == 4 else f" {a.name:<3}"
        lines.append(f"ATOM  {k:5d} {name} {a.resname:>3} {a.chain:1}{a.resseq:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.symbol:>2}")
    return "\n".join(lines) + "\n"


def test_every_ionisable_group_gets_the_right_hydrogens():
    mol = experimental("ADEHKRCYG")
    add_hydrogens_at_ph(mol, 7.4)
    assert hydrogens_on(mol, 1, "N") == 3          # NH3+ at the N-terminus
    assert hydrogens_on(mol, 2, "OD1", "OD2") == 0
    assert hydrogens_on(mol, 3, "OE1", "OE2") == 0
    assert hydrogens_on(mol, 4, "ND1", "NE2") == 1   # one neutral tautomer
    assert hydrogens_on(mol, 5, "NZ") == 3
    assert (hydrogens_on(mol, 6, "NE"), hydrogens_on(mol, 6, "NH1"), hydrogens_on(mol, 6, "NH2")) == (1, 2, 2)
    assert hydrogens_on(mol, 7, "SG") == 1
    assert hydrogens_on(mol, 8, "OH") == 1
    assert hydrogens_on(mol, 9, "O", "OXT") == 0
    backbone = [hydrogens_on(mol, k, "N") for k in range(2, 10)]
    assert backbone == [1] * 8
    assert all(hydrogens_on(mol, k, "CA") == (2 if k == 9 else 1) for k in range(1, 10))


def test_histidine_is_doubly_protonated_below_its_pka():
    mol = experimental("AHA")
    add_hydrogens_at_ph(mol, 5.0)
    assert hydrogens_on(mol, 2, "ND1") == 1 and hydrogens_on(mol, 2, "NE2") == 1


def test_variant_names_fix_the_state():
    mol = experimental("AHDA")
    for a in mol.atoms:
        a.resname = {"HIS": "HID", "ASP": "ASH"}.get(a.resname, a.resname)
    add_hydrogens_at_ph(mol, 7.4)
    assert hydrogens_on(mol, 2, "ND1") == 1 and hydrogens_on(mol, 2, "NE2") == 0
    assert hydrogens_on(mol, 3, "OD2") == 1


def test_metal_bound_histidine_and_disulfide_cysteine_carry_no_hydrogen():
    mol = experimental("AHCCA")
    ne2 = next(i for i, a in enumerate(mol.atoms) if a.resseq == 2 and a.name == "NE2")
    zinc = mol.add_atom(Atom("Zn", resname="ZN", resseq=101, chain="A", hetero=True),
                        mol.coords[ne2] + np.array([2.0, 0.0, 0.0]))
    mol.add_bond(ne2, zinc, 1.0)
    sg = [i for i, a in enumerate(mol.atoms) if a.name == "SG"]
    mol.add_bond(sg[0], sg[1], 1.0)
    add_hydrogens_at_ph(mol, 7.4)
    assert hydrogens_on(mol, 2, "NE2") == 0 and hydrogens_on(mol, 2, "ND1") == 1
    assert hydrogens_on(mol, 3, "SG") == 0 and hydrogens_on(mol, 4, "SG") == 0
