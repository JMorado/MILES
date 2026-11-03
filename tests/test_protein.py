import numpy as np

from miles.core.document import Document
from miles.core.molecule import Atom, Molecule
from miles.core.perception import assign_bond_orders, infer_atom_types, perceive_bonds
from miles.core.protein import default_representations, estimate_secondary_structure
from miles.render.cartoon import cartoon_mesh


def finish(m: Molecule, has_bonds: bool) -> None:
    if not has_bonds:
        perceive_bonds(m)
        assign_bond_orders(m)
    infer_atom_types(m, only_missing=False)


def ideal_helix(n: int = 14) -> Molecule:
    """Backbone CA and O atoms on an ideal alpha helix (3.6 residues per turn, 1.5 A rise)."""
    atoms, coords = [], []
    for k in range(n):
        t = np.radians(100.0 * k)
        ca = np.array([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * k])
        atoms.append(Atom("C", name="CA", resname="ALA", resseq=k + 1, chain="A"))
        coords.append(ca)
        atoms.append(Atom("O", name="O", resname="ALA", resseq=k + 1, chain="A"))
        coords.append(ca + np.array([0.0, 0.0, 1.2]))
    return Molecule.from_atoms(atoms, coords)


def test_helix_is_detected_from_ca_geometry():
    m = ideal_helix()
    estimate_secondary_structure(m)
    assert sum(a.ss == "H" for a in m.atoms if a.name == "CA") >= 10


def test_default_representation_and_mesh():
    m = ideal_helix()
    m.add_atom(Atom("O", name="O", resname="HOH", resseq=200, chain="A", hetero=True), [20.0, 0.0, 0.0])
    finish(m, False)
    reps = default_representations(m)
    assert reps is not None and reps[-1] == "hidden"
    mesh = cartoon_mesh(m, reps == "cartoon")
    assert not mesh.empty
    assert mesh.indices.max() < len(mesh.vertices)
    assert np.allclose(np.linalg.norm(mesh.normals, axis=1), 1.0, atol=1e-4)


def test_peptide_bonds_keep_the_carbonyl():
    from miles.core import perception
    from miles.core.chem import from_smiles
    mol, _ = from_smiles("CC(=O)NC")            # N-methylacetamide
    mol.set_bonds({})
    for a in mol.atoms:
        a.geom = a.nsites = a.valence = None
    perception.perceive_bonds(mol)
    perception.assign_bond_orders(mol)
    perception.infer_atom_types(mol, only_missing=False)
    c = next(i for i, a in enumerate(mol.atoms)
             if a.symbol == "C" and any(mol.atoms[n].symbol == "O" for n in mol.neighbors(i)))
    o = next(n for n in mol.neighbors(c) if mol.atoms[n].symbol == "O")
    n = next(k for k in mol.neighbors(c) if mol.atoms[k].symbol == "N")
    assert mol.bond_order(c, o) == 2.0
    assert mol.bond_order(c, n) == 1.0
    assert mol.atoms[n].geom == "trigonal"


def test_gentle_clean_keeps_experimental_structure():
    from miles.core.forcefield import clean_geometry
    m = ideal_helix(20)
    finish(m, False)
    before = m.coords.copy()
    assert clean_geometry(m) == "gentle"
    assert np.sqrt(((m.coords - before) ** 2).sum(axis=1).mean()) < 0.5


def test_clean_selection_leaves_surroundings_alone():
    from miles.core.forcefield import clean_selection
    from miles.core.fragments import ring_fragment
    m = ideal_helix(12)
    finish(m, False)
    ring = ring_fragment("benzene").mol
    offset = m.extend(ring)
    m.set_positions(range(offset, len(m)), m.coords[offset:] + np.array([6.0, 0.0, 8.0]))
    sel = list(range(offset, len(m)))
    m.set_positions(sel, m.coords[sel] + np.random.default_rng(2).normal(0, 0.2, (len(sel), 3)))
    before = m.coords.copy()
    clean_selection(m, sel)
    assert np.allclose(m.coords[:offset], before[:offset])
    cc = [np.linalg.norm(m.coords[i] - m.coords[j]) for (i, j) in m.bonds
          if i >= offset and m.atoms[i].symbol == m.atoms[j].symbol == "C"]
    assert np.allclose(cc, 1.40, atol=0.02)


def test_new_hydrogens_follow_their_residue():
    from miles.core import building, perception
    from miles.core.fragments import ring_fragment
    m = ideal_helix(8)
    finish(m, False)
    doc = Document(m)
    doc.display = doc.display.with_rep(default_representations(m))
    ring = ring_fragment("benzene").mol
    ring.remove_atoms([i for i, a in enumerate(ring.atoms) if a.symbol == "H"])
    for a in ring.atoms:
        a.resname, a.resseq, a.chain, a.hetero = "BNZ", 200, "A", True
        a.geom = a.nsites = a.valence = None
    offset = doc.edit(Molecule.extend, ring)
    perception.infer_atom_types(m)
    lig = list(range(offset, len(m)))
    start = len(m)
    assert doc.edit(building.add_hydrogens, lig) == 6
    new = m.atoms[start:]
    assert all(h.resname == "BNZ" and h.hetero for h in new)
    assert not doc.display.rep[start:].any()   # the ligand has no representation of its own
    doc.edit(building.add_hydrogens)
    on_ca = [k for k in range(start + 6, len(m)) if m.atoms[k].resname == "ALA"]
    assert on_ca and all(doc.display.rep[k] == "cartoon" for k in on_ca)   # hydrogens drawn like their atom
