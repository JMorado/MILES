import numpy as np

from miles.builtin_plugins.binding_pocket.interactions import (
    InteractionKind,
    aromatic_rings,
    find_interactions,
    ligand_atoms,
    pocket_atoms,
)
from miles.core.chem import from_smiles
from miles.core.molecule import Atom, Molecule
from miles.files.formats import parse_text


def pdb_line(serial: int, name: str, resname: str, resseq: int, xyz: tuple[float, float, float],
             element: str, hetero: bool = False) -> str:
    record = "HETATM" if hetero else "ATOM"
    label = f" {name:<3}" if len(name) < 4 else name
    x, y, z = xyz
    return (f"{record:<6}{serial:>5} {label}{'':1}{resname:>3} A{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}")


LYSINE_AND_ACETATE = "\n".join([
    pdb_line(1, "CB", "LYS", 1, (-5.0, 0.0, 0.0), "C"),
    pdb_line(2, "CG", "LYS", 1, (-3.8, 0.8, 0.0), "C"),
    pdb_line(3, "CD", "LYS", 1, (-2.5, 0.0, 0.0), "C"),
    pdb_line(4, "CE", "LYS", 1, (-1.3, 0.8, 0.0), "C"),
    pdb_line(5, "NZ", "LYS", 1, (0.0, 0.0, 0.0), "N"),
    pdb_line(6, "C1", "ACT", 50, (3.2, 0.0, 0.0), "C", hetero=True),
    pdb_line(7, "O1", "ACT", 50, (2.6, 1.08, 0.0), "O", hetero=True),
    pdb_line(8, "O2", "ACT", 50, (2.6, -1.08, 0.0), "O", hetero=True),
    pdb_line(9, "C2", "ACT", 50, (4.7, 0.0, 0.0), "C", hetero=True),
    pdb_line(10, "O", "HOH", 90, (0.0, 3.0, 0.0), "O", hetero=True),
    "END",
]) + "\n"


def lysine_and_acetate() -> tuple[Molecule, list[int]]:
    mol = parse_text(LYSINE_AND_ACETATE)
    ligand = ligand_atoms(mol, next(i for i, a in enumerate(mol.atoms) if a.resname == "ACT"))
    return mol, ligand


def test_ligand_and_pocket_are_whole_residues_without_water():
    mol, ligand = lysine_and_acetate()
    assert sorted(mol.atoms[i].name for i in ligand) == ["C1", "C2", "O1", "O2"]
    pocket = pocket_atoms(mol, ligand)
    assert {mol.atoms[i].resname for i in pocket} == {"LYS"}
    assert len(pocket) == 5  # the whole lysine, though only some atoms are within 5 angstrom


def test_salt_bridge_and_hydrogen_bonds():
    mol, ligand = lysine_and_acetate()
    found = find_interactions(mol, ligand, pocket_atoms(mol, ligand))
    kinds = [f.kind for f in found]
    assert InteractionKind.SALT_BRIDGE in kinds
    bridge = next(f for f in found if f.kind == InteractionKind.SALT_BRIDGE)
    assert sorted(mol.atoms[i].name for i in bridge.first) == ["O1", "O2"]
    assert [mol.atoms[i].name for i in bridge.second] == ["NZ"]
    assert abs(bridge.distance - 2.6) < 0.01
    hbonds = [f for f in found if f.kind == InteractionKind.HYDROGEN_BOND]
    assert len(hbonds) == 2  # NZ to both oxygens, accepted on distance (no hydrogens in the file)
    assert len(found) == len({(f.kind, f.first, f.second) for f in found})


def water_pair(towards: bool) -> Molecule:
    hydrogen = (0.96, 0.0, 0.0) if towards else (-0.96, 0.0, 0.0)
    atoms = [Atom("O"), Atom("H"), Atom("H"), Atom("O"), Atom("H"), Atom("H")]
    coords = np.array([(0, 0, 0), hydrogen, (-0.24, 0.93, 0), (2.9, 0, 0), (3.2, 0.9, 0), (3.2, -0.9, 0)],
                      dtype=float)
    mol = Molecule.from_atoms(atoms, coords)
    for i, j in ((0, 1), (0, 2), (3, 4), (3, 5)):
        mol.add_bond(i, j)
    return mol


def test_hydrogen_bond_needs_a_hydrogen_pointing_at_the_acceptor():
    good = water_pair(towards=True)
    found = find_interactions(good, [0, 1, 2], [3, 4, 5])
    assert [(f.kind, f.first, f.second) for f in found] == [(InteractionKind.HYDROGEN_BOND, (0,), (3,))]
    assert find_interactions(water_pair(towards=False), [0, 1, 2], [3, 4, 5]) == []


def test_stacked_benzene_rings():
    first = from_smiles("c1ccccc1")[0]
    carbons = [i for i, a in enumerate(first.atoms) if a.symbol == "C"]
    centered = first.coords - first.coords[carbons].mean(axis=0)
    normal = np.linalg.svd(centered[carbons])[2][2]
    first.coords = centered
    second = first.copy()
    second.coords = second.coords + 3.7 * normal
    mol = first.copy()
    offset = mol.extend(second)
    ligand = list(range(offset))
    pocket = list(range(offset, len(mol)))
    assert len(aromatic_rings(mol, ligand)) == 1
    found = find_interactions(mol, ligand, pocket)
    assert [f.kind for f in found] == [InteractionKind.PI_STACKING]
    assert abs(found[0].distance - 3.7) < 1e-6
    assert len(found[0].first) == len(found[0].second) == 6


def test_unrelated_atoms_do_not_interact():
    mol, ligand = lysine_and_acetate()
    mol.set_positions(ligand, mol.coords[ligand] + 20.0)
    assert pocket_atoms(mol, ligand) == []
    assert find_interactions(mol, ligand, []) == []
