import time

import numpy as np
import pytest

from miles.builtin_plugins.grids.combine import Combination, combine_volumes, same_points
from miles.core.volume import (
    DENSITY_LEVEL,
    ORBITAL_LEVEL,
    Volume,
    isosurface,
    mapped_isosurface,
    signed_isosurfaces,
    symmetric_range,
)
from miles.files import session
from miles.files.common import JobInfo
from miles.files.formats import read_file

BOHR = 0.529177210903


def write_cube(path, values, origin_bohr, step_bohr, atoms, comment="test grid"):
    """Write a Gaussian cube file with an orthogonal grid (all lengths in bohr)."""
    n = values.shape
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"test job\n{comment}\n")   # the second line says what the grid holds
        f.write(f"{len(atoms):5d}" + "".join(f"{x:12.6f}" for x in origin_bohr) + "\n")
        for axis in range(3):
            v = [0.0, 0.0, 0.0]
            v[axis] = step_bohr
            f.write(f"{n[axis]:5d}" + "".join(f"{x:12.6f}" for x in v) + "\n")
        for z, xyz in atoms:
            f.write(f"{z:5d}{float(z):12.6f}" + "".join(f"{x:12.6f}" for x in xyz) + "\n")
        for i in range(n[0]):
            for j in range(n[1]):
                row = values[i, j]
                for k in range(0, n[2], 6):
                    f.write("".join(f"{x:13.5E}" for x in row[k:k + 6]) + "\n")


def edge_list(faces):
    """Every triangle edge as a sorted vertex pair."""
    return np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)


def grid(n, step):
    """Coordinates (angstrom) of an n^3 grid centred on the origin."""
    origin = -(n - 1) / 2 * step
    axis = origin + step * np.arange(n)
    return np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1), origin


def blob_volume(n=40, step=0.2, width=1.0):
    points, origin = grid(n, step)
    values = np.exp(-(points ** 2).sum(axis=-1) / width ** 2).astype(np.float32)
    return Volume("blob", np.full(3, origin), np.eye(3) * step, values)


def test_cube_file_is_read_with_its_grid(tmp_path):
    n, step = 16, 0.4
    points, origin = grid(n, step)
    values = np.exp(-(points ** 2).sum(axis=-1)) * np.sign(points[..., 2] + 1e-9)
    path = tmp_path / "orbital.cube"
    atoms = [(8, [0.0, 0.0, 0.0]), (1, [0.0, 1.43, -1.11]), (1, [0.0, -1.43, -1.11])]
    write_cube(path, values, [origin / BOHR] * 3, step / BOHR, atoms, comment="orbital 5")
    mol, info = read_file(str(path))
    assert mol.symbols == ["O", "H", "H"]
    assert np.allclose(mol.coords[1], np.array([0.0, 1.43, -1.11]) * BOHR, atol=1e-4)
    (volume,) = info.volumes
    assert volume.name == "orbital 5"
    assert volume.values.shape == (n, n, n)
    assert np.allclose(volume.origin, origin, atol=1e-4)
    assert np.allclose(volume.axes, np.eye(3) * step, atol=1e-5)
    assert np.allclose(volume.values, values, atol=1e-4)
    assert volume.signed


def test_isosurface_of_a_blob_is_a_closed_sphere():
    volume = blob_volume()
    level = 0.5
    mesh = isosurface(volume, level)
    radius = np.sqrt(np.log(1 / level))
    distances = np.linalg.norm(mesh.vertices, axis=1)
    assert np.allclose(distances, radius, atol=0.02)
    _, counts = np.unique(edge_list(mesh.faces), axis=0, return_counts=True)
    assert (counts == 2).all()
    assert (np.einsum("ij,ij->i", mesh.normals, mesh.vertices) > 0).all()
    tri = mesh.vertices[mesh.faces]
    face_normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    assert (np.einsum("ij,ij->i", face_normals, tri.mean(axis=1)) > 0).mean() > 0.99
    assert (mesh.atoms == -1).all()


def test_lobe_cut_by_the_grid_edge_is_closed():
    points, origin = grid(20, 0.2)
    values = np.exp(-((points - [1.9, 0, 0]) ** 2).sum(axis=-1)).astype(np.float32)
    mesh = isosurface(Volume("edge", np.full(3, origin), np.eye(3) * 0.2, values), 0.3)
    _, counts = np.unique(edge_list(mesh.faces), axis=0, return_counts=True)
    assert (counts == 2).all()


def test_signed_lobes_of_an_antisymmetric_function():
    points, origin = grid(36, 0.2)
    z = points[..., 2]
    values = (z * np.exp(-(points ** 2).sum(axis=-1))).astype(np.float32)
    volume = Volume("p orbital", np.full(3, origin), np.eye(3) * 0.2, values)
    positive, negative = signed_isosurfaces(volume, 0.1)
    assert (positive.vertices[:, 2] > 0).all() and (negative.vertices[:, 2] < 0).all()
    assert abs(len(positive.vertices) - len(negative.vertices)) <= 0.02 * len(positive.vertices)
    for lobe, sign in ((positive, 1), (negative, -1)):
        centre = np.array([0.0, 0.0, sign * 0.7])
        outward = np.einsum("ij,ij->i", lobe.normals, lobe.vertices - centre)
        assert (outward > 0).mean() > 0.95


def test_non_orthogonal_grid():
    axes = np.array([[0.2, 0.0, 0.0], [0.1, 0.2, 0.0], [0.0, 0.05, 0.2]])
    n = 50
    index = np.stack(np.meshgrid(*[np.arange(n)] * 3, indexing="ij"), axis=-1).reshape(-1, 3)
    origin = -(index.mean(axis=0) @ axes)
    points = origin + index @ axes
    values = np.exp(-(points ** 2).sum(axis=-1)).reshape(n, n, n).astype(np.float32)
    mesh = isosurface(Volume("skewed", origin, axes, values), 0.5)
    assert np.allclose(np.linalg.norm(mesh.vertices, axis=1), np.sqrt(np.log(2)), atol=0.03)
    assert (np.einsum("ij,ij->i", mesh.normals, mesh.vertices) > 0).all()


def test_suggested_levels():
    points, origin = grid(10, 0.5)
    density = np.exp(-(points ** 2).sum(axis=-1)) * 50
    orbital = points[..., 0] * np.exp(-(points ** 2).sum(axis=-1))
    tiny = density * 1e-5
    assert Volume("density", origin, np.eye(3), density).suggested_level() == DENSITY_LEVEL
    assert Volume("orbital", origin, np.eye(3), orbital).suggested_level() == ORBITAL_LEVEL
    assert 0 < Volume("tiny", origin, np.eye(3), tiny).suggested_level() < tiny.max()
    assert not Volume("density", origin, np.eye(3), density).signed
    assert Volume("zero", origin, np.eye(3), np.zeros((4, 4, 4))).suggested_level() == 0


def test_isosurface_of_a_large_grid_is_fast():
    volume = blob_volume(n=100, step=0.1, width=2.0)
    start = time.perf_counter()
    mesh = isosurface(volume, 0.3)
    elapsed = time.perf_counter() - start
    print(f"100^3 grid: {elapsed:.2f} s, {len(mesh.faces)} triangles")
    assert not mesh.empty
    assert elapsed < 2.0


def test_session_keeps_volumes(tmp_path):
    volume = blob_volume(n=8, step=0.5)
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("O")
    path = tmp_path / "grid.miles"
    session.write(str(path), mol, JobInfo(title="water", volumes=[volume]))
    _, info = session.read(str(path))
    (restored,) = info.volumes
    assert restored.name == "blob"
    assert restored.values.shape == volume.values.shape
    assert np.allclose(restored.values, volume.values)
    assert np.allclose(restored.axes, volume.axes) and np.allclose(restored.origin, volume.origin)


def test_old_sessions_without_volumes_still_load(tmp_path):
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("C")
    data = session.to_dict(mol, JobInfo(title="old"))
    del data["volumes"], data["vibrations"]
    restored, info = session.from_dict(data)
    assert len(restored) == len(mol)
    assert info.volumes == [] and info.vibrations == []


@pytest.mark.parametrize("level", [0.2, -0.2])
def test_flat_grid_gives_no_surface(level):
    volume = Volume("flat", np.zeros(3), np.eye(3), np.zeros((5, 5, 5), np.float32))
    assert isosurface(volume, level).empty


def test_sample_is_exact_for_linear_fields():
    axes = np.array([[0.5, 0.0, 0.0], [0.1, 0.4, 0.0], [0.0, 0.0, 0.3]])
    origin = np.array([1.0, 2.0, 3.0])
    index = np.indices((5, 6, 7)).astype(float)

    def field(i, j, k):
        return 2 * i - j + 0.5 * k + 1

    volume = Volume("linear", origin, axes, field(*index).astype(np.float32))
    fractional = np.array([[1.3, 2.7, 4.1], [0.0, 0.0, 0.0], [3.99, 4.5, 5.5]])
    values = volume.sample(origin + fractional @ axes)
    assert np.allclose(values, [field(*f) for f in fractional], atol=1e-5)
    outside = volume.sample(origin + np.array([[10.0, 10.0, 10.0]]) @ axes)   # clamped to the far corner
    assert np.allclose(outside, field(4, 5, 6))


def test_mapped_isosurface_carries_the_other_grid():
    density = blob_volume()
    points, origin = grid(40, 0.2)
    potential = Volume("z", np.full(3, origin), np.eye(3) * 0.2, points[..., 2].astype(np.float32))
    surface, values = mapped_isosurface(density, 0.5, potential)
    assert not surface.empty and len(values) == len(surface.vertices)
    assert np.allclose(values, surface.vertices[:, 2], atol=1e-4)   # trilinear is exact for a linear field
    radius = np.sqrt(np.log(2.0))   # where exp(-r^2) = 0.5
    assert 0.8 * radius < symmetric_range(values) <= radius + 0.05


def test_symmetric_range_ignores_outliers():
    values = np.concatenate([np.linspace(-1.0, 1.0, 1000), [50.0]])
    assert symmetric_range(values) < 1.0
    assert symmetric_range([]) == 1.0
    assert symmetric_range(np.zeros(4)) == 1.0


def test_mapped_volume_mesh_is_coloured_by_sign():
    from miles.render.surface import MAPPED_COLORS, volume_mesh

    density = blob_volume()
    points, origin = grid(40, 0.2)
    potential = Volume("z", np.full(3, origin), np.eye(3) * 0.2, points[..., 2].astype(np.float32))
    mesh = volume_mesh(density, 0.5, potential, mapped_range=0.5)
    top = mesh.colors[mesh.vertices[:, 2] > 0.5]
    bottom = mesh.colors[mesh.vertices[:, 2] < -0.5]
    assert len(top) and len(bottom)
    assert np.allclose(top, MAPPED_COLORS[2], atol=1e-6)      # positive: blue
    assert np.allclose(bottom, MAPPED_COLORS[0], atol=1e-6)   # negative: red


def test_difference_of_a_grid_with_itself_is_zero_and_sum_doubles():
    a = blob_volume()
    zero = combine_volumes(a, a)
    assert np.abs(zero.values).max() == 0.0 and zero.signed and zero.name == "blob - blob"
    double = combine_volumes(a, a, Combination.SUM, "twice")
    assert np.allclose(double.values, 2 * a.values) and double.name == "twice" and not double.signed


def test_grids_on_other_points_are_resampled():
    a = blob_volume(n=30, step=0.2)
    shifted = Volume("shifted", a.origin + np.array([0.3, 0.0, 0.0]), a.axes.copy(), a.values.copy())
    assert not same_points(a, shifted)
    difference = combine_volumes(a, shifted)
    nx, ny, nz = a.values.shape
    points = a.origin + np.indices((nx, ny, nz)).reshape(3, -1).T @ a.axes
    expected = a.values.reshape(-1) - shifted.sample(points)
    assert np.allclose(difference.values.reshape(-1), expected, atol=1e-6)
    assert difference.values.max() > 0.1 and difference.values.min() < -0.1   # a shifted blob: two lobes


def test_session_keeps_the_sign_marker(tmp_path):
    from miles.core.chem import from_smiles

    a = blob_volume(n=8, step=0.5)
    mol, _ = from_smiles("O")
    path = tmp_path / "grids.miles"
    session.write(str(path), mol, JobInfo(volumes=[a, combine_volumes(a, a)]))
    _, info = session.read(str(path))
    assert [v.both_signs for v in info.volumes] == [None, True]
