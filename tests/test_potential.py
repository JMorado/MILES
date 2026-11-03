"""Electrostatic potential of partial charges on molecular surfaces."""

import time

import numpy as np
import pytest

from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
from miles.core import potential
from miles.core.charges import ChargeModel, partial_charges
from miles.core.chem import from_smiles
from miles.core.molecule import Molecule
from miles.core.potential import (
    COULOMB,
    PROBE_OFFSET,
    Dielectric,
    default_dielectric,
    electrostatic_potential,
    surface_potential,
)
from miles.core.surface import SurfaceColoring, SurfaceMesh, ionisable_charges, molecular_surface
from miles.render.surface import CHARGE_COLORS, SurfaceRequest, surface_mesh
from miles.render.theme import LIGHT


def naive(points, coords, charges, dielectric):
    """Reference sum, one point and one charge at a time."""
    result = []
    for p in points:
        total = 0.0
        for c, q in zip(coords, charges, strict=True):
            d = max(float(np.linalg.norm(p - c)), potential.MIN_DISTANCE)
            total += q / (4 * d * d if dielectric == Dielectric.DISTANCE else d)
        result.append(COULOMB * total)
    return np.array(result)


def sphere_surface(n, radius=20.0):
    """A sphere of ``n`` vertices with outward normals (no faces are needed to evaluate)."""
    rng = np.random.default_rng(3)
    normals = rng.normal(size=(n, 3))
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    return SurfaceMesh(vertices=(radius * normals).astype(np.float32), normals=normals.astype(np.float32),
                       faces=np.zeros((1, 3), np.uint32), atoms=np.zeros(n, int))


def molecule_of(coords, symbol="C"):
    mol = Molecule()
    for xyz in coords:
        mol.add_atom(symbol, xyz)
    return mol


@pytest.mark.parametrize("dielectric", list(Dielectric))
def test_sum_matches_a_direct_calculation(dielectric):
    rng = np.random.default_rng(1)
    coords = rng.uniform(-5, 5, (30, 3))
    charges = rng.normal(size=30) * 0.4
    points = rng.uniform(-8, 8, (50, 3))
    expected = naive(points, coords, charges, dielectric)
    assert np.allclose(electrostatic_potential(points, coords, charges, dielectric), expected,
                       rtol=1e-4, atol=1e-3)


def test_single_charge_and_close_points():
    value = electrostatic_potential([[0.0, 0.0, 2.0]], [[0.0, 0.0, 0.0]], [1.0])
    assert value[0] == pytest.approx(COULOMB / 2.0, rel=1e-5)
    screened = electrostatic_potential([[0.0, 0.0, 2.0]], [[0.0, 0.0, 0.0]], [1.0], Dielectric.DISTANCE)
    assert screened[0] == pytest.approx(COULOMB / 16.0, rel=1e-5)
    on_top = electrostatic_potential([[0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]], [1.0])
    assert np.isfinite(on_top).all()
    assert not electrostatic_potential(np.zeros((3, 3)), np.zeros((2, 3)), [0.0, 0.0]).any()


def test_surface_potential_is_taken_outside_the_surface():
    mol, _ = from_smiles("CC(=O)[O-]")
    charges = partial_charges(mol, ChargeModel.GASTEIGER)
    surface = molecular_surface(mol)
    values = surface_potential(mol, surface, charges)
    points = surface.vertices + PROBE_OFFSET * surface.normals
    assert np.allclose(values, naive(points[:40], mol.coords, charges, Dielectric.VACUUM).tolist()
                       + values[40:].tolist(), rtol=1e-4, atol=1e-3)


@pytest.mark.parametrize(("smiles", "symbol", "sign"), [("CC(=O)[O-]", "O", -1), ("C[NH3+]", "N", 1)])
def test_sign_near_charged_groups(smiles, symbol, sign):
    mol, _ = from_smiles(smiles)
    charges = partial_charges(mol, ChargeModel.GASTEIGER)
    surface = molecular_surface(mol)
    values = surface_potential(mol, surface, charges)
    centres = mol.coords[[i for i, s in enumerate(mol.symbols) if s == symbol]]
    near = np.min(np.linalg.norm(surface.vertices[:, None, :] - centres[None], axis=2), axis=1) < 2.2
    assert near.any()
    assert sign * values[near].mean() > 5.0   # kcal/(mol e)


def test_interpolation_follows_the_exact_potential(monkeypatch):
    rng = np.random.default_rng(5)
    directions = rng.normal(size=(400, 3))
    radii = 15 * rng.random((400, 1)) ** (1 / 3)
    coords = directions / np.linalg.norm(directions, axis=1, keepdims=True) * radii
    mol = molecule_of(coords)                      # charges inside a ball, at least 3.4 A from the points
    charges = rng.normal(size=400) * 0.3
    surface = sphere_surface(20000, radius=17.0)
    exact = surface_potential(mol, surface, charges, Dielectric.DISTANCE)
    monkeypatch.setattr(potential, "EXACT_LIMIT", 400 * 3000)   # at most 3000 exact points
    approximate = surface_potential(mol, surface, charges, Dielectric.DISTANCE)
    error = np.abs(approximate - exact)
    assert np.corrcoef(approximate, exact)[0, 1] > 0.98
    assert np.percentile(error, 95) < 0.15 * np.percentile(np.abs(exact), 95)


def test_large_surface_is_fast():
    rng = np.random.default_rng(7)
    mol = molecule_of(rng.uniform(-25, 25, (9000, 3)))
    charges = rng.normal(size=9000) * 0.3
    surface = sphere_surface(120_000, radius=30.0)
    start = time.perf_counter()
    values = surface_potential(mol, surface, charges, Dielectric.DISTANCE)
    elapsed = time.perf_counter() - start
    assert values.shape == (120_000,) and np.isfinite(values).all()
    assert elapsed < 3.0   # about 1 s on a laptop; generous for slow machines


def test_dielectric_and_charges_for_polymers():
    peptide = build_sequence("KAAAAD", SequenceKind.HELIX)
    small, _ = from_smiles("CCO")
    assert default_dielectric(peptide) == Dielectric.DISTANCE
    assert default_dielectric(small) == Dielectric.VACUUM
    charges = ionisable_charges(peptide)
    names = [(a.resname, a.name) for a in peptide.atoms]
    assert charges[names.index(("LYS", "NZ"))] == pytest.approx(1.0)
    assert charges[[i for i, n in enumerate(names) if n[0] == "ASP"]].sum() < 0
    ion = molecule_of([[0.0, 0.0, 0.0]], "Mg")
    assert ionisable_charges(ion)[0] == 2.0


def test_surface_mesh_colours_by_potential():
    mol, _ = from_smiles("C[NH3+]")
    charges = partial_charges(mol, ChargeModel.GASTEIGER)
    request = SurfaceRequest(None, SurfaceColoring.POTENTIAL)
    mesh = surface_mesh(mol, request, LIGHT, charges)
    assert mesh.value_range > 0
    blue = np.abs(mesh.colors - CHARGE_COLORS[2]).sum(axis=1) < 0.05
    assert blue.any()                                    # a cation is positive all round
    plain = surface_mesh(mol, request, LIGHT)            # without charges: one neutral colour
    assert plain.value_range == 0 and np.allclose(plain.colors, plain.colors[0])
