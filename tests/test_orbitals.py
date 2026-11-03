import json
import math
import time
from pathlib import Path

import numpy as np
import pytest

from miles.core.orbitals import (
    DENSITY_NAME,
    SPHERICAL,
    SPIN_NAME,
    Orbital,
    OrbitalSet,
    Shell,
    Spin,
    contracted_shell,
    solid_harmonics,
)
from miles.core.units import BOHR
from miles.core.volume import DENSITY_LEVEL, Volume
from miles.files import session
from miles.files.formats import read_file
from miles.files.molden import read_molden_orbitals

DATA = Path(__file__).parent / "data"
H2 = DATA / "h2_sto3g.molden"

STO3G_EXPONENTS = [3.42525091, 0.62391373, 0.16885540]
STO3G_COEFFICIENTS = [0.15432897, 0.53532814, 0.44463454]


def integral(a: Volume, b: Volume | None = None) -> float:
    """Integral of a * b (or a squared) over the grid, in atomic units."""
    other = a if b is None else b
    cell = abs(np.linalg.det(a.axes)) / BOHR**3
    return float((a.values.astype(float) * other.values.astype(float)).sum() * cell)


def single_function(shell: str, index: int, flags: str = "", title: str = "") -> str:
    """Molden text for one atom with one primitive shell and one orbital on basis function ``index``."""
    return f"""[Molden Format]
[Title]
{title}
{flags}
[Atoms] Angs
C  1  6  0.0  0.0  0.0
[GTO]
1 0
{shell} 1 1.0
0.8 1.0

[MO]
Sym= A
Ene= -0.5
Spin= Alpha
Occup= 2.0
{index + 1} 1.0
"""


def test_h2_orbitals_are_read():
    orbitals = read_molden_orbitals(H2.read_text())
    assert orbitals is not None
    assert [o.energy for o in orbitals.orbitals] == [-0.5782, 0.6703]
    assert [o.occupation for o in orbitals.orbitals] == [2.0, 0.0]
    assert [o.symmetry for o in orbitals.orbitals] == ["1sg", "1su"]
    assert not orbitals.unrestricted
    assert orbitals.coefficients.shape == (2, 2)
    assert np.allclose(orbitals.atoms, [[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]])
    assert orbitals.homo() == 0
    assert [orbitals.label(k) for k in range(2)] == ["HOMO", "LUMO"]
    assert orbitals.name(1) == "MO 2 (LUMO)"


def test_h2_orbitals_are_orthonormal_on_a_grid():
    orbitals = read_molden_orbitals(H2.read_text())
    assert orbitals is not None
    bonding = orbitals.evaluate(0, spacing=0.1)
    antibonding = orbitals.evaluate(1, spacing=0.1)
    assert bonding.name == "MO 1 (HOMO)"
    assert np.allclose(bonding.origin, [-4.0, -4.0, -4.0 - 0.7 * BOHR])
    assert np.allclose(bonding.axes, np.eye(3) * 0.1)
    assert integral(bonding) == pytest.approx(1.0, abs=1e-3)
    assert integral(antibonding) == pytest.approx(1.0, abs=1e-3)
    assert integral(bonding, antibonding) == pytest.approx(0.0, abs=1e-3)
    assert bonding.values.min() > -1e-6
    assert antibonding.values.max() == pytest.approx(-antibonding.values.min(), rel=0.05)


def test_coefficients_of_bare_primitives_give_the_same_basis():
    norms = [(2 * a / math.pi) ** 0.75 for a in STO3G_EXPONENTS]
    pairs = zip(STO3G_EXPONENTS, STO3G_COEFFICIENTS, norms, strict=True)
    rows = "\n".join(f"{a} {c * n}" for a, c, n in pairs)
    text = H2.read_text()
    bare = text.split("[GTO]")[0] + f"[GTO]\n1 0\ns 3 1.0\n{rows}\n\n2 0\ns 3 1.0\n{rows}\n\n[MO]" \
        + text.split("[MO]")[1]
    standard, scaled = read_molden_orbitals(text), read_molden_orbitals(bare)
    assert standard is not None and scaled is not None
    for a, b in zip(standard.shells, scaled.shells, strict=True):
        assert np.allclose(a.coefficients, b.coefficients)


def test_shell_scale_factor_scales_the_exponents():
    text = single_function("s", 0).replace("s 1 1.0", "s 1 2.0")
    orbitals = read_molden_orbitals(text)
    assert orbitals is not None
    assert orbitals.shells[0].exponents == pytest.approx([3.2])


@pytest.mark.parametrize("angmom", [2, 3, 4, 5])
def test_solid_harmonics_are_orthonormal(angmom: int):
    def moment(powers: tuple[int, ...]) -> float:
        if any(p % 2 for p in powers):
            return 0.0
        return float(np.prod([math.prod(range(p - 1, 0, -2)) for p in powers]))

    functions = SPHERICAL[angmom]
    reference = moment((2 * angmom, 0, 0))
    gram = np.array([[sum(fa * fb * moment(tuple(p + q for p, q in zip(pa, pb, strict=True)))
                          for pa, fa in a.items() for pb, fb in b.items()) / reference
                      for b in functions] for a in functions])
    assert np.allclose(gram, np.eye(2 * angmom + 1))


@pytest.mark.parametrize(("shell", "flags", "index", "size"), [
    ("d", "", 3, 6),           # cartesian xy
    ("d", "", 0, 6),           # cartesian xx
    ("d", "[5D]", 4, 5),       # d(-2)
    ("f", "", 9, 10),          # cartesian xyz
    ("f", "[7F]", 5, 7),       # f(+3)
    ("g", "", 9, 15),          # cartesian xxyy
    ("g", "[9G]", 8, 9),       # g(-4)
])
def test_single_functions_are_normalised_on_a_grid(shell: str, flags: str, index: int, size: int):
    orbitals = read_molden_orbitals(single_function(shell, index, flags))
    assert orbitals is not None
    assert orbitals.coefficients.shape == (1, size)
    volume = orbitals.evaluate(0, spacing=0.08, margin=5.0)
    assert integral(volume) == pytest.approx(1.0, abs=2e-3)


def test_orca_files_flip_high_m_spherical_functions():
    plain = read_molden_orbitals(single_function("f", 5, "[7F]"))
    title = "Molden file created by orca_2mkl for BaseName=x"
    orca = read_molden_orbitals(single_function("f", 5, "[7F]", title))
    assert plain is not None and orca is not None
    assert orca.coefficients[0, 5] == -plain.coefficients[0, 5]
    low = read_molden_orbitals(single_function("f", 2, "[7F]", "orca_2mkl"))
    assert low is not None and low.coefficients[0, 2] == 1.0


def test_files_without_a_basis_have_no_orbitals():
    assert read_molden_orbitals("[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n") is None
    assert read_molden_orbitals(H2.read_text().split("[MO]")[0]) is None


def test_inconsistent_coefficients_are_rejected():
    with pytest.raises(ValueError, match="basis function 3"):
        read_molden_orbitals(H2.read_text() + "   3    0.1\n")


def test_labels_count_from_the_frontier_orbitals():
    energies = [-1.0, -0.5, -0.4, 0.1, 0.3]
    orbitals = OrbitalSet([Orbital(e, 2.0 if e < 0 else 0.0) for e in energies], [], np.zeros((5, 0)),
                          np.zeros((1, 3)))
    assert orbitals.homo() == 2
    assert [orbitals.label(k) for k in range(5)] == ["HOMO-2", "HOMO-1", "HOMO", "LUMO", "LUMO+1"]
    assert orbitals.label(5) == ""


def test_unrestricted_labels_name_the_spin():
    orbitals = OrbitalSet([Orbital(-0.6, 1.0, Spin.ALPHA), Orbital(-0.4, 1.0, Spin.ALPHA),
                           Orbital(0.2, 0.0, Spin.ALPHA), Orbital(-0.5, 1.0, Spin.BETA),
                           Orbital(0.1, 0.0, Spin.BETA), Orbital(0.4, 0.0, Spin.BETA)],
                          [], np.zeros((6, 0)), np.zeros((1, 3)))
    assert orbitals.unrestricted
    assert orbitals.homo() == 1
    assert [orbitals.label(k) for k in range(6)] == [
        "alpha HOMO-1", "alpha HOMO", "alpha LUMO", "beta HOMO", "beta LUMO", "beta LUMO+1"]
    assert orbitals.name(4) == "MO 5 (beta LUMO)"


def test_read_file_attaches_the_orbitals():
    mol, info = read_file(str(H2))
    assert mol.symbols == ["H", "H"]
    assert np.allclose(mol.coords, [[0.0, 0.0, -0.7 * BOHR], [0.0, 0.0, 0.7 * BOHR]], atol=1e-4)
    assert info.orbitals is not None
    assert len(info.orbitals.orbitals) == 2


def test_orbitals_survive_a_session_round_trip():
    mol, info = read_file(str(H2))
    data = json.loads(json.dumps(session.to_dict(mol, info)))
    _, restored = session.from_dict(data)
    assert restored.orbitals is not None and info.orbitals is not None
    assert restored.orbitals.orbitals == info.orbitals.orbitals
    assert np.allclose(restored.orbitals.coefficients, info.orbitals.coefficients, atol=1e-6)
    before, after = info.orbitals.evaluate(1), restored.orbitals.evaluate(1)
    assert after.name == before.name
    assert np.allclose(after.values, before.values, atol=1e-6)
    data.pop("orbitals")
    assert session.from_dict(data)[1].orbitals is None


def test_evaluation_is_fast_for_a_medium_molecule():
    rng = np.random.default_rng(3)
    atoms = rng.uniform(0.0, 7.0 / BOHR, size=(20, 3))
    shells: list[Shell] = []
    for center in atoms:
        shells += [
            contracted_shell(center, 0, np.array([3047.5, 457.4, 103.9, 29.2, 9.29, 3.16]),
                             np.array([0.0018, 0.014, 0.069, 0.23, 0.47, 0.36])),
            contracted_shell(center, 0, np.array([7.87, 1.88, 0.54]), np.array([-0.12, -0.17, 1.15])),
            contracted_shell(center, 0, np.array([0.17]), np.array([1.0])),
            contracted_shell(center, 1, np.array([7.87, 1.88, 0.54]), np.array([0.07, 0.32, 0.74])),
            contracted_shell(center, 1, np.array([0.17]), np.array([1.0])),
            contracted_shell(center, 2, np.array([0.8]), np.array([1.0]), spherical=True),
        ]
    n = sum(s.size for s in shells)
    orbitals = OrbitalSet([Orbital(-0.3, 2.0)], shells, rng.normal(size=(1, n)) / math.sqrt(n), atoms)
    start = time.perf_counter()
    volume = orbitals.evaluate(0)
    elapsed = time.perf_counter() - start
    assert min(volume.values.shape) >= 55
    assert elapsed < 1.0


def test_solid_harmonics_match_the_textbook():
    tables = solid_harmonics(4)
    assert tables[2][0] == pytest.approx({(0, 0, 2): 1.0, (2, 0, 0): -0.5, (0, 2, 0): -0.5})
    assert tables[2][4] == pytest.approx({(1, 1, 0): math.sqrt(3)})                      # d, m = -2
    f5 = math.sqrt(5 / 8)
    assert tables[3][5] == pytest.approx({(3, 0, 0): f5, (1, 2, 0): -3 * f5})             # f, m = +3
    assert tables[4][0] == pytest.approx({(0, 0, 4): 1.0, (4, 0, 0): 3 / 8, (0, 4, 0): 3 / 8,
                                          (2, 2, 0): 3 / 4, (2, 0, 2): -3.0, (0, 2, 2): -3.0})


def _one_shell_molden(kind: str, size: int, title: str) -> str:
    """One atom with one shell and one normalised orbital whose coefficients grow as 1, 2, 3, ..."""
    flags = {"f": "[7F]", "g": "[9G]", "h": ""}[kind]
    lines = ["[Molden Format]", "[Title]", title, "[Atoms] AU", "Cu 1 29 0.0 0.0 0.0", "[GTO]", "  1 0",
             f" {kind}    1 1.00", "   1.50000000    1.00000000", "", flags, "[MO]", " Sym= 1a", " Ene= -0.5",
             " Spin= Alpha", " Occup= 2.0"]
    norm = math.sqrt(sum((k + 1) ** 2 for k in range(size)))   # a spherical shell is orthonormal
    lines += [f" {k + 1} {(k + 1) / norm:.10f}" for k in range(size)]
    return "\n".join(line for line in lines if line != "") + "\n"


@pytest.mark.parametrize(("kind", "size"), [("f", 7), ("g", 9), ("h", 11)])
def test_orca_files_have_m_3_and_4_flipped(kind, size):
    plain = read_molden_orbitals(_one_shell_molden(kind, size, "written by hand"))
    orca = read_molden_orbitals(_one_shell_molden(kind, size, "Molden file created by orca_2mkl"))
    assert plain is not None and orca is not None and plain.shells[0].size == size
    signs = np.sign(orca.coefficients[0]) * np.sign(plain.coefficients[0])
    flipped = {5, 6, 7, 8} & set(range(size))   # m = +3, -3, +4, -4 in Molden order
    assert [k for k in range(size) if signs[k] < 0] == sorted(flipped)


EXAMPLES = Path(__file__).parent.parent / "examples" / "qm"


def grid_integral(volume, spacing):
    return float(volume.values.sum()) * (spacing / BOHR) ** 3


def test_electron_density_integrates_to_the_electrons():
    orbitals = read_molden_orbitals(H2.read_text())
    density = orbitals.density()
    assert density.name == DENSITY_NAME and not density.signed
    assert grid_integral(density, 0.25) == pytest.approx(2.0, abs=0.01)
    formaldehyde = read_molden_orbitals((EXAMPLES / "formaldehyde.molden").read_text())
    assert grid_integral(formaldehyde.density(spacing=0.15), 0.15) == pytest.approx(16.0, abs=0.5)
    with pytest.raises(ValueError, match="alpha and beta"):
        orbitals.spin_density()


def test_spin_density_of_an_open_shell_molecule_holds_one_electron():
    orbitals = read_molden_orbitals((EXAMPLES / "nitrogen_dioxide.molden").read_text())
    spin = orbitals.spin_density()
    assert spin.name == SPIN_NAME
    assert grid_integral(spin, 0.25) == pytest.approx(1.0, abs=0.05)   # one unpaired electron
    assert spin.values.min() < 0 < spin.values.max()
    assert spin.signed and spin.suggested_level() == pytest.approx(DENSITY_LEVEL)
