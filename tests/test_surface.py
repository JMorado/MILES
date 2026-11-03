import time

import numpy as np
import pytest

from miles.core.chem import from_smiles
from miles.core.elements import element
from miles.core.fragments import amino_acid_template
from miles.core.molecule import Atom, Molecule
from miles.core.surface import SurfaceColoring, SurfaceMesh, atom_property, molecular_surface


def area(mesh: SurfaceMesh) -> float:
    v = mesh.vertices[mesh.faces].astype(float)
    return float(0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1).sum())


def edge_counts(mesh: SurfaceMesh) -> np.ndarray:
    f = mesh.faces
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    return np.unique(edges, axis=0, return_counts=True)[1]


def components(mesh: SurfaceMesh) -> int:
    label = np.arange(len(mesh.vertices))
    f = mesh.faces.astype(np.int64)
    while True:
        low = np.minimum(np.minimum(label[f[:, 0]], label[f[:, 1]]), label[f[:, 2]])
        new = label.copy()
        for k in range(3):
            np.minimum.at(new, f[:, k], low)
        new = new[new]
        if np.array_equal(new, label):
            return len(np.unique(label))
        label = new


def outward_fraction(mesh: SurfaceMesh, mol: Molecule) -> float:
    center = mol.coords[mesh.atoms]
    return float((np.einsum("ij,ij->i", mesh.normals, mesh.vertices - center) > 0).mean())


def single_atom(symbol: str = "C") -> Molecule:
    return Molecule.from_atoms([Atom(symbol)], np.zeros((1, 3)))


def test_single_atom_is_a_closed_sphere():
    mol = single_atom()
    mesh = molecular_surface(mol)
    r = element("C").vdw_radius
    radii = np.linalg.norm(mesh.vertices, axis=1)
    assert np.allclose(radii, r, atol=0.06)
    assert area(mesh) == pytest.approx(4 * np.pi * r * r, rel=0.05)
    assert (edge_counts(mesh) == 2).all()
    assert components(mesh) == 1
    assert outward_fraction(mesh, mol) == 1.0
    assert np.allclose(np.linalg.norm(mesh.normals, axis=1), 1.0, atol=1e-4)
    assert mesh.faces.dtype == np.uint32 and mesh.vertices.dtype == np.float32


def test_faces_wind_counter_clockwise_from_outside():
    mesh = molecular_surface(single_atom())
    v = mesh.vertices[mesh.faces].astype(float)
    face_normals = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    assert (np.einsum("ij,ij->i", face_normals, v.mean(axis=1)) > 0).all()


def test_distant_atoms_give_separate_surfaces():
    mol = Molecule.from_atoms([Atom("O"), Atom("O")], np.array([[0.0, 0, 0], [25.0, 0, 0]]))
    mesh = molecular_surface(mol)
    assert components(mesh) == 2
    assert set(mesh.atoms.tolist()) == {0, 1}
    assert (edge_counts(mesh) == 2).all()


def test_close_atoms_merge_into_one_surface():
    mol = Molecule.from_atoms([Atom("C"), Atom("C")], np.array([[0.0, 0, 0], [1.54, 0, 0]]))
    assert components(molecular_surface(mol)) == 1


def test_subset_of_atoms():
    mol = Molecule.from_atoms([Atom("C"), Atom("N")], np.array([[0.0, 0, 0], [20.0, 0, 0]]))
    mesh = molecular_surface(mol, atoms=[1])
    assert set(mesh.atoms.tolist()) == {1}
    assert molecular_surface(mol, atoms=[]).empty


def test_benzene_is_watertight_and_similar_without_hydrogens():
    benzene = from_smiles("c1ccccc1")[0]
    start = time.perf_counter()
    mesh = molecular_surface(benzene)
    assert time.perf_counter() - start < 0.5
    assert (edge_counts(mesh) == 2).all()
    assert components(mesh) == 1
    assert outward_fraction(mesh, benzene) > 0.95
    carbons = [i for i, a in enumerate(benzene.atoms) if a.symbol == "C"]
    without = molecular_surface(benzene, atoms=carbons)
    assert area(without) == pytest.approx(area(mesh), rel=0.2)


def test_peptide_surface():
    peptide = from_smiles("N" + "C(C)C(=O)N" * 29 + "CC(=O)O")[0]
    start = time.perf_counter()
    mesh = molecular_surface(peptide)
    assert time.perf_counter() - start < 2.0
    assert components(mesh) == 1
    assert (edge_counts(mesh) != 2).mean() < 0.001     # rare ambiguous corners of Surface Nets
    assert outward_fraction(mesh, peptide) > 0.95


def test_large_cloud_is_fast():
    n = 5000
    rng = np.random.default_rng(0)
    side = int(np.ceil(n ** (1 / 3)))
    lattice = np.stack(np.meshgrid(*[np.arange(side)] * 3, indexing="ij"), axis=-1).reshape(-1, 3)[:n] * 1.6
    mol = Molecule.from_atoms([Atom("C") for _ in range(n)], lattice + rng.uniform(-0.1, 0.1, lattice.shape))
    start = time.perf_counter()
    mesh = molecular_surface(mol)
    elapsed = time.perf_counter() - start
    print(f"surface of {n} atoms: {elapsed:.2f} s, {len(mesh.faces)} triangles")
    assert elapsed < 5.0
    assert not mesh.empty


def residues() -> Molecule:
    """Lys, Asp and Leu side by side, numbered 1-3 in chain A."""
    mol = Molecule()
    for n, code in enumerate(("LYS", "ASP", "LEU"), start=1):
        piece = amino_acid_template(code).mol.copy()
        piece.coords = piece.coords + np.array([12.0 * n, 0.0, 0.0])
        for a in piece.atoms:
            a.resname, a.resseq, a.chain = code, n, "A"
        mol.extend(piece)
    return mol


def values_of(mol: Molecule, values: np.ndarray, resname: str, name: str) -> float:
    return float(next(values[i] for i, a in enumerate(mol.atoms) if a.resname == resname and a.name == name))


def test_hydrophobicity_signs():
    mol = residues()
    values = atom_property(mol, SurfaceColoring.HYDROPHOBICITY)
    assert values_of(mol, values, "LEU", "CD1") > 0
    assert values_of(mol, values, "LYS", "NZ") < 0
    assert values_of(mol, values, "ASP", "CG") < 0
    assert np.all(np.abs(values) <= 1)


def test_charge_signs():
    mol = residues()
    values = atom_property(mol, SurfaceColoring.CHARGE)
    assert values_of(mol, values, "LYS", "NZ") > 0
    assert values_of(mol, values, "ASP", "OD1") < 0
    assert values_of(mol, values, "ASP", "OD2") < 0
    assert values_of(mol, values, "LEU", "CD1") == 0
    hydrogens_on_nz = [i for i, a in enumerate(mol.atoms) if a.resname == "LYS" and a.name == "HZ"]
    assert all(values[i] > 0 for i in hydrogens_on_nz)


def test_charges_of_small_molecule_groups():
    acetate = from_smiles("CC(=O)[O-]")[0]
    values = atom_property(acetate, SurfaceColoring.CHARGE)
    oxygens = [i for i, a in enumerate(acetate.atoms) if a.symbol == "O"]
    assert all(values[i] < 0 for i in oxygens)
    ammonium = from_smiles("C[NH3+]")[0]
    values = atom_property(ammonium, SurfaceColoring.CHARGE)
    assert values[next(i for i, a in enumerate(ammonium.atoms) if a.symbol == "N")] > 0


def test_element_and_uniform_colouring_are_zero():
    mol = residues()
    for coloring in (SurfaceColoring.ELEMENT, SurfaceColoring.UNIFORM):
        assert not atom_property(mol, coloring).any()
