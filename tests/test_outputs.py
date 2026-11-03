"""Reading quantum chemistry outputs: what Open Babel misses or gets wrong, and Gaussian input.

Every input here is a few hand-written lines in the layout of the program
concerned, reproducing a case where reading went wrong.
"""

import math
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from miles.core.geometry import angle, dihedral, distance
from miles.core.orbitals import Shell, contracted_shell, overlap_matrix
from miles.files import gaussian, outputs
from miles.files.common import FileFormatError
from miles.files.formats import parse_text, read_file
from miles.files.molden import read_molden_orbitals
from miles.files.openbabel import HARTREE, input_format

WATER = [8, 1, 1]
GAMESS_HEADER = ("\n COORDINATES OF ALL ATOMS ARE (ANGS)\n"
                 "   ATOM   CHARGE       X              Y              Z\n"
                 " ------------------------------------------------------------\n")


def water_rows(shift: float, fmt: str) -> str:
    """The three atoms of water in the row layout of a program, moved by ``shift`` along z."""
    z = (0.12 + shift, -0.48 - shift, -0.48 - shift)
    y = (0.0, 0.76, -0.76)
    if fmt == "gaussian":
        return "".join(f"{k + 1:7d}{n:11d}{0:12d}{0.0:16.6f}{y[k]:12.6f}{z[k]:12.6f}\n"
                       for k, n in enumerate(WATER))
    if fmt == "gamess":
        return "".join(f" {s:11s}{float(n):4.1f}{0.0:15.10f}{y[k]:15.10f}{z[k]:15.10f}\n"
                       for k, (s, n) in enumerate(zip("OHH", WATER, strict=True)))
    return "".join(f"  {s}{0.0:14.6f}{y[k]:12.6f}{z[k]:12.6f}\n" for k, s in enumerate("OHH"))


# -- steps of a job --------------------------------------------------------------------------------

def test_files_that_report_empty_structures_at_the_end_are_read(tmp_path):
    # Open Babel's MOPAC and GAMESS input readers report empty structures
    # after the last one forever; reading must stop (it hung before)
    path = tmp_path / "water.mop"
    path.write_text(" PM7\nwater\n\n  O 0.0 1 0.0 1 0.12 1\n"
                    "  H 0.0 1 0.76 1 -0.48 1\n  H 0.0 1 -0.76 1 -0.48 1\n")
    code = f"from miles.files.formats import read_file; print(len(read_file({str(path)!r})[0]))"
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert done.returncode == 0 and done.stdout.split()[-1] == "3"


def test_firefly_optimisation_steps_are_read():
    # Firefly writes a bare "NSERCH=" line before each step; its "NSERCH= k ENERGY="
    # summaries and the geometry printed again at the end are not steps
    header = GAMESS_HEADER
    text = ""
    for step, shift in enumerate([0.02, 0.01, 0.0]):
        text += f"\n NSERCH=   {step}\n" + header + water_rows(shift, "gamess")
        text += f"\n          NSERCH=  {step}     ENERGY=    -76.40{step}\n"
    text += "\n      ***** EQUILIBRIUM GEOMETRY LOCATED *****\n" + header + water_rows(0.0, "gamess")
    steps = outputs.gamess_steps(text, WATER)
    assert len(steps) == 3
    assert steps[0][0, 2] == pytest.approx(0.14) and steps[-1][0, 2] == pytest.approx(0.12)
    assert outputs.gamess_steps(text, [8, 1, 1, 1]) == []   # blocks with other atoms are skipped


def test_gamess_us_optimisation_steps_are_read():
    header = GAMESS_HEADER
    text = ""
    for step, shift in enumerate([0.02, 0.0]):
        text += f"\n BEGINNING GEOMETRY SEARCH POINT NSERCH=   {step} ...\n"
        text += header + water_rows(shift, "gamess")
        text += f"\n NSERCH:   {step}  E=     -76.40{step}  GRAD. MAX=  0.001\n"
    text += "\n      ***** EQUILIBRIUM GEOMETRY LOCATED *****\n" + header + water_rows(0.0, "gamess")
    assert len(outputs.gamess_steps(text, WATER)) == 2


def orca_scan(steps: int) -> str:
    text = "                                 * O   R   C   A *\n"
    for step in range(steps):
        text += (f"\n         *               RELAXED SURFACE SCAN STEP   {step + 1}               *\n\n"
                 "---------------------------------\nCARTESIAN COORDINATES (ANGSTROEM)\n"
                 "---------------------------------\n" + water_rows(0.01 * (steps - 1 - step), "orca")
                 + f"\nFINAL SINGLE POINT ENERGY       -76.4{step}0000000\n")
    return text


def test_orca_scans_show_the_last_geometry_with_every_step(tmp_path):
    # Open Babel kept only the first geometry of an ORCA scan
    path = tmp_path / "scan.out"
    path.write_text(orca_scan(3))
    mol, info = read_file(str(path))
    assert info.frames.shape == (3, 3, 3)
    assert mol.coords[0, 2] == pytest.approx(0.12)
    assert info.energies == pytest.approx([(-76.40 - 0.01 * k) * HARTREE for k in range(3)])


def gaussian_steps(blocks: int, energies: int) -> str:
    text = " Entering Gaussian System\n"
    for k in range(blocks):
        text += ("                         Standard orientation:\n"
                 " ---------------------------------------------------------------------\n"
                 " Center     Atomic      Atomic             Coordinates (Angstroms)\n"
                 " Number     Number       Type             X           Y           Z\n"
                 " ---------------------------------------------------------------------\n"
                 + water_rows(0.01 * k, "gaussian")
                 + " ---------------------------------------------------------------------\n")
        if k < energies:
            text += f" SCF Done:  E(RHF) =  -76.0{k}00000000     A.U. after    9 cycles\n"
    return text


def test_gaussian_unrelaxed_scans_keep_every_point():
    steps = outputs.job_steps(gaussian_steps(3, 3), "g09", WATER)
    assert steps is not None
    coords, energies = steps
    assert coords.shape == (3, 3, 3)
    assert energies == pytest.approx([(-76.0 - 0.01 * k) * HARTREE for k in range(3)])
    # an optimisation prints its last geometry again without a new energy: Open Babel's steps are kept
    assert outputs.job_steps(gaussian_steps(3, 2), "g09", WATER) is None


def test_qchem_steps_charge_and_multiplicity(tmp_path):
    def rows(step: int) -> list[str]:
        return water_rows(0.01 * (1 - step), "orca").splitlines(True)

    text = " Welcome to Q-Chem\n$molecule\n +1 2\n O 0 0 0.12\n H 0 0.76 -0.48\n H 0 -0.76 -0.48\n$end\n"
    for k in range(2):
        text += ("             Standard Nuclear Orientation (Angstroms)\n"
                 "    I     Atom           X                Y                Z\n"
                 " ----------------------------------------------------------------\n"
                 + "".join(f"    {n + 1}  {row}" for n, row in enumerate(rows(k)))
                 + " ----------------------------------------------------------------\n"
                 f" Total energy in the final basis set =      -76.4{k}000000\n")
    path = tmp_path / "opt.out"
    path.write_text(text)
    _, info = read_file(str(path))
    assert info.frames.shape == (2, 3, 3) and len(info.energies) == 2
    assert (info.charge, info.multiplicity) == (1, 2)


@pytest.mark.parametrize(("fmt", "text", "expected"), [
    ("gamout", " CHARGE OF MOLECULE                           =    1\n"
               " SPIN MULTIPLICITY                            =    2\n", (1, 2)),
    ("gamout", " CHARGE OF MOLECULE                  =   -2\n"
               "     MULT  =       1     ICHARG=      -2\n", (-2, 1)),
    ("nwo", "          Charge           :     1\n          Spin multiplicity:     2\n", (1, 2)),
    ("orca", "Total Charge           Charge          ....    1\n"
             "Multiplicity           Mult            ....    2\n", (1, 2)),
    ("jout", "  net molecular charge:    1\n  multiplicity:            2\n", (1, 2)),
    ("mopout", " *  SEXTET     - SPIN STATE DEFINED AS A SEXTET\n *    CHARGE ON SYSTEM =    -3\n", (-3, 6)),
    ("mopout", " PM7 water\n", (None, None)),
    ("qcout", "$molecule\n -1 1\n Cl 0 0 0\n$end\n", (-1, 1)),
])
def test_charge_and_multiplicity_of_outputs(fmt, text, expected):
    assert outputs.charge_and_multiplicity(text, fmt) == expected


# -- ORCA normal modes -----------------------------------------------------------------------------

def orca_frequencies(symmetry: bool, intensities_in_km: bool) -> str:
    """Normal modes of a diatomic in ORCA's layout: five zero modes and the stretch along z."""
    lines = ["-----------------------", "VIBRATIONAL FREQUENCIES", "-----------------------", ""]
    lines += [f"   {k}:       0.00 cm**-1" for k in range(5)]
    lines += ["   5:    2100.50 cm**-1", "", "------------", "NORMAL MODES", "------------", ""]
    lines.append("            " + "".join(f"{k:11d}" for k in range(6)))
    if symmetry:
        labels = ["1-Ag", "1-Bu", "1-Au", "2-Bu", "2-Au", "2-Ag"]
        lines.append("            " + "".join(f"{label:>11s}" for label in labels))
    for row in range(6):
        values = [1.0 if k == row and k < 5 else 0.0 for k in range(5)]
        values.append(0.7071 if row == 2 else -0.7071 if row == 5 else 0.0)
        lines.append(f"      {row:2d}  " + "".join(f"{v:11.6f}" for v in values))
    lines += ["", "-----------", "IR SPECTRUM", "-----------", ""]
    if intensities_in_km:
        lines += [" Mode   freq       eps      Int      T**2         TX        TY        TZ",
                  "       cm**-1   L/(mol*cm) km/mol    a.u.",
                  "----------------------------------------------------------------------------",
                  "  5:   2100.50   0.012345   62.40  0.002900  ( 0.000000  0.000000  0.053852)"]
    else:
        lines += [" Mode    freq (cm**-1)   T**2         TX         TY         TZ",
                  "-------------------------------------------------------------------",
                  "  5:   2100.50   0.002900  (  0.000000   0.000000   0.053852)"]
    return "\n".join(lines) + "\n\n"


def test_orca_modes_take_the_intensity_in_km_per_mol():
    # Open Babel read the eps column of ORCA 5 and later as the intensity
    modes = outputs.orca_modes(orca_frequencies(symmetry=False, intensities_in_km=True), 2)
    assert modes is not None and len(modes) == 1
    assert modes[0].frequency == pytest.approx(2100.5) and modes[0].intensity == pytest.approx(62.4)
    assert modes[0].displacements[:, 2] == pytest.approx([0.7071, -0.7071])


def test_orca_modes_labelled_by_symmetry_are_read():
    # ORCA 6 adds a row of symmetry labels under the mode numbers, where Open Babel read no modes
    modes = outputs.orca_modes(orca_frequencies(symmetry=True, intensities_in_km=True), 2)
    assert modes is not None and len(modes) == 1 and modes[0].intensity == pytest.approx(62.4)


def test_orca_modes_without_intensities_in_km_per_mol_have_none():
    modes = outputs.orca_modes(orca_frequencies(symmetry=False, intensities_in_km=False), 2)
    assert modes is not None and modes[0].intensity is None


# -- formatted checkpoint files --------------------------------------------------------------------

def fchk(geometries: int = 1, modes: bool = False) -> str:
    coords = [0.0, 0.0, 0.2267671, 0.0, 1.4362330, -0.9070685, 0.0, -1.4362330, -0.9070685]

    def array(name: str, kind: str, values: list[float]) -> str:
        rows = "".join("".join(f"{v:16.8E}" for v in values[k:k + 5]) + "\n"
                       for k in range(0, len(values), 5))
        return f"{name:43s}{kind}   N={len(values):12d}\n" + rows

    text = "water\nSP        RHF                                                         STO-3G\n"
    for k in range(geometries):
        text += f"{'Number of atoms':43s}I{3:17d}\n{'Charge':43s}I{1:17d}\n{'Multiplicity':43s}I{2:17d}\n"
        text += f"{'Atomic numbers':43s}I   N={3:12d}\n{8:12d}{1:12d}{1:12d}\n"
        scale = 1 + 0.01 * (geometries - 1 - k)
        text += array("Current cartesian coordinates", "R", [c * scale for c in coords])
        text += f"{'Total Energy':43s}R{-74.96 - 0.01 * k:27.15E}\n"
    text += array("Mulliken Charges", "R", [-0.6, 0.3, 0.3])
    text += array("Dipole Moment", "R", [0.0, 0.0, 0.629265633])
    if modes:
        e2 = [1600.0, 3700.0, 3800.0] + [1.0] * 3 + [1.0] * 3 + [70.0, 2.0, 20.0]
        vectors = [0.0, 0.0, 0.07, 0.0, 0.4, -0.55, 0.0, -0.4, -0.55] * 3
        text += array("Vib-E2", "R", e2) + array("Vib-Modes", "R", vectors)
    return text


def test_checkpoint_charges_dipole_and_modes(tmp_path):
    path = tmp_path / "water.fchk"
    path.write_text(fchk(modes=True))
    _, info = read_file(str(path))
    assert (info.charge, info.multiplicity) == (1, 2)
    assert info.charges == pytest.approx([-0.6, 0.3, 0.3])
    # the file gives the dipole in atomic units; Open Babel passed it on as if in debye
    assert float(np.linalg.norm(info.dipole)) == pytest.approx(1.5994, abs=1e-3)
    assert [m.frequency for m in info.vibrations] == pytest.approx([1600.0, 3700.0, 3800.0])
    assert [m.intensity for m in info.vibrations] == pytest.approx([70.0, 2.0, 20.0])


def test_checkpoint_with_every_step_shows_the_last(tmp_path):
    # Q-Chem writes the records again for every step; Open Babel read only the first
    path = tmp_path / "opt.fchk"
    path.write_text(fchk(geometries=3))
    mol, info = read_file(str(path))
    assert info.frames.shape == (3, 3, 3) and len(info.energies) == 3
    assert mol.coords[0, 2] == pytest.approx(0.2267671 * 0.529177, abs=1e-5)


# -- detecting the program -------------------------------------------------------------------------

def test_gamess_uk_is_not_read_as_gamess(tmp_path):
    path = tmp_path / "job.out"
    path.write_text("          GAMESS-UK  version 8.0\n")
    assert input_format(str(path)) == "gukout"


GAMESS_UK = textwrap.dedent("""\
                             GAMESS-UK
                                            molecular geometry
             *******************************************************************************
             *     atom   atomic                coordinates                 number of      *
             *            charge       x             y              z       shells         *
             *******************************************************************************
             *                                                                             *
             *                                                                             *
             *    o         8.0   0.0000000      0.0000000      0.2267671       3          *
             *                                                                1s  2sp      *
             *                                                                             *
             *                                                                             *
             *    h         1.0   0.0000000      1.4362330     -0.9070685       1          *
             *                                                                1s           *
             *                                                                             *
             *                                                                             *
             *    h         1.0   0.0000000     -1.4362330     -0.9070685       1          *
             *                                                                1s           *
             *                                                                             *
             *******************************************************************************

                         cartesians to normal mode coordinates transformation matrix

                        1           2           3
          1   x   0.1   0.2
""")


def test_gamess_uk_frequency_outputs_do_not_abort(tmp_path):
    # Open Babel's GAMESS-UK reader aborts the whole program on the normal modes
    path = tmp_path / "freq.out"
    path.write_text(GAMESS_UK)
    code = f"from miles.files.formats import read_file; print(len(read_file({str(path)!r})[0]))"
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert done.returncode == 0 and done.stdout.split()[-1] == "3"


def test_outputs_of_programs_open_babel_cannot_read_say_so(tmp_path):
    path = tmp_path / "job.out"
    path.write_text("          Psi4: An Open-Source Ab Initio Electronic Structure Package\n")
    with pytest.raises(FileFormatError, match="output of Psi4"):
        read_file(str(path))


def test_formats_open_babel_only_writes_are_not_read_with_it(tmp_path):
    assert input_format(str(tmp_path / "job.gjf")) is None   # Gaussian input has its own reader
    assert input_format(str(tmp_path / "job.gzmat")) == "gzmat"


def test_files_with_no_structure_say_how_they_were_read(tmp_path):
    path = tmp_path / "job.inp"
    path.write_text("! B3LYP def2-SVP\n* xyz 0 1\nO 0 0 0\n*\n")   # ORCA input, read as GAMESS input
    with pytest.raises(FileFormatError, match="read as GAMESS Input"):
        read_file(str(path))


# -- Molden orbitals in other conventions ----------------------------------------------------------

def one_d_shell(coefficients: list[list[float]]) -> str:
    """A hydrogen atom with one cartesian d shell and the given orbitals (xx, yy, zz, xy, xz, yz)."""
    lines = ["[Molden Format]", "[Atoms] AU", "H 1 1 0.0 0.0 0.0", "[GTO]", "  1 0", " d    1 1.00",
             "   1.0000000000      1.0000000000", "", "[MO]"]
    for k, vector in enumerate(coefficients):
        lines += [" Sym= A", f" Ene= {0.5 + k}", " Spin= Alpha", " Occup= 0.0"]
        lines += [f" {i + 1} {v:.12f}" for i, v in enumerate(vector)]
    return "\n".join(lines) + "\n"


# Two normalised orbitals with each function normalised: d_xy, and d_x2-y2
# (<xx|yy> = 1/3, so xx - yy has norm sqrt(4/3))
ORBITALS = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], [math.sqrt(0.75), -math.sqrt(0.75), 0.0, 0.0, 0.0, 0.0]]
SCALES = {  # how each program's xx and xy functions relate to normalised ones
    "Molden": (1.0, 1.0), "CFOUR": (math.sqrt(3.0), 1.0), "Psi4 1.3": (1.0, 1.0 / math.sqrt(3.0)),
    "Turbomole": (math.sqrt(3.0), math.sqrt(3.0)),
}


@pytest.mark.parametrize("program", list(SCALES))
def test_cartesian_conventions_of_programs_give_normalised_orbitals(program):
    xx, xy = SCALES[program]
    factors = [xx, xx, xx, xy, xy, xy]
    written = [[c / f for c, f in zip(vector, factors, strict=True)] for vector in ORBITALS]
    found = read_molden_orbitals(one_d_shell(written))
    assert found is not None
    overlap = overlap_matrix(found.shells)
    norms = np.einsum("ij,jk,ik->i", found.coefficients, overlap, found.coefficients)
    assert norms == pytest.approx([1.0, 1.0], abs=1e-9)
    assert found.coefficients[0, 3] == pytest.approx(1.0)


def test_orbitals_in_an_unknown_convention_are_left_out_with_a_note(tmp_path):
    with pytest.raises(ValueError, match="not normalised"):
        read_molden_orbitals(one_d_shell([[2.0 * c for c in vector] for vector in ORBITALS]))
    path = tmp_path / "h.molden"
    path.write_text(one_d_shell([[2.0 * c for c in vector] for vector in ORBITALS]))
    mol, info = read_file(str(path))
    assert info.orbitals is None and len(mol) == 1
    assert info.notes and "orbitals were left out" in info.notes[0]


def test_overlap_matrix_matches_the_formula_for_two_s_functions():
    a, b, r = 0.8, 1.3, 1.7
    shells = [contracted_shell(np.zeros(3), 0, np.array([a]), np.array([1.0])),
              contracted_shell(np.array([0.0, 0.0, r]), 0, np.array([b]), np.array([1.0]))]
    s = overlap_matrix(shells)
    expected = (2.0 * math.sqrt(a * b) / (a + b)) ** 1.5 * math.exp(-a * b / (a + b) * r * r)
    assert s == pytest.approx(np.array([[1.0, expected], [expected, 1.0]]))


@pytest.mark.parametrize(("angmom", "spherical"),
                         [(1, False), (2, False), (2, True), (3, True), (4, False), (5, True)])
def test_overlap_matrix_has_unit_diagonal(angmom, spherical):
    shell = contracted_shell(np.array([0.1, -0.2, 0.3]), angmom, np.array([2.0, 0.5]), np.array([0.4, 0.7]),
                             spherical)
    other = Shell(np.array([0.5, 0.5, 0.5]), 0, np.array([1.0]), np.array([1.0 / (math.pi / 2.0) ** 0.75]))
    s = overlap_matrix([shell, other])
    assert np.diag(s) == pytest.approx(np.ones(shell.size + 1))
    assert s == pytest.approx(s.T)


# -- Gaussian input --------------------------------------------------------------------------------

def test_gaussian_input_with_cartesian_coordinates(tmp_path):
    path = tmp_path / "water.gjf"
    path.write_text("%chk=water.chk\n%mem=1GB\n#p hf/sto-3g opt   ! comment\n\nwater molecule\n\n0 1\n"
                    "O     0.000000   0.000000   0.117300\nH  -1  0.000000   0.757200  -0.469200\n"
                    "H,0.000000,-0.757200,-0.469200\n\n")
    mol, info = read_file(str(path))
    assert mol.formula() == "H2O" and len(mol.bonds) == 2
    assert (info.title, info.charge, info.multiplicity) == ("water molecule", 0, 1)


def test_gaussian_input_with_a_z_matrix_and_variables():
    text = ("#p b3lyp/6-31g(d) opt\n\nhydrogen peroxide\n\n0 1\nO1\nO2 O1 roo\nH3 O1 roh O2 aooh\n"
            "H4 O2 roh O1 aooh H3 dih\nX5 O1 1.0 O2 90.0 H3 0.0\n\n"
            "roo 1.45\nroh 0.97\naooh 100.0\ndih 120.0\n")
    _, _, _, symbols, coords = gaussian.parse(text)
    assert symbols == ["O", "O", "H", "H"]   # the dummy atom is left out
    assert distance(coords[0], coords[1]) == pytest.approx(1.45)
    assert angle(coords[2], coords[0], coords[1]) == pytest.approx(100.0)
    assert abs(dihedral(coords[3], coords[1], coords[0], coords[2])) == pytest.approx(120.0)


def test_gaussian_input_details():
    bohr = ("#p hf/sto-3g units=bohr\n\nh2\n\n0 1\nH 0 0 0\nH 0 0 1.4\n\n"
            "--Link1--\n#p hf\n\nnext\n\n0 1\nHe 0 0 0\n")
    _, _, _, symbols, coords = gaussian.parse(bohr)
    assert symbols == ["H", "H"] and coords[1, 2] == pytest.approx(1.4 * 0.529177, abs=1e-5)
    ghost = "#p hf\n\nghost\n\n0 1\nHe 0 0 0\nBq 0 0 2\nNe-Bq 0 0 4\n"
    assert gaussian.parse(ghost)[3] == ["He"]
    pasted = "#p hf/sto-3g\n\nwater\n\n0 1\nO 0 0 0.12\nH 0 0.76 -0.48\nH 0 -0.76 -0.48\n"
    assert parse_text(pasted).formula() == "H2O"
    with pytest.raises(FileFormatError, match="route"):
        gaussian.parse("***, water\ngeometry={O;H,O,0.96}\n")   # a Molpro input, also .com
    with pytest.raises(FileFormatError, match="not defined"):
        gaussian.parse("#p hf\n\nt\n\n0 1\nH\nH 1 r\n")


def test_molden_labels_in_capitals_are_elements(tmp_path):
    # Psi4 writes "CU"; Open Babel read it as a dummy atom, so the orbitals were dropped too
    path = tmp_path / "cuh.molden"
    path.write_text("[Molden Format]\n[Atoms] AU\nCU 1 29 0.0 0.0 0.0\nH 2 1 0.0 0.0 2.9\n[GTO]\n  1 0\n"
                    " s    1 1.00\n   1.0000000000      1.0000000000\n\n  2 0\n s    1 1.00\n"
                    "   1.0000000000      1.0000000000\n\n"
                    "[MO]\n Sym= A\n Ene= -0.5\n Spin= Alpha\n Occup= 2.0\n"
                    " 1 1.0\n 2 0.0\n")
    mol, info = read_file(str(path))
    assert mol.symbols == ["Cu", "H"] and len(mol.bonds) == 1
    assert info.orbitals is not None and not info.notes


def test_multiplicity_is_that_of_the_ground_state():
    # a TDDFT input asks for triplet excited states of a singlet with MULT=3
    gamess = (" INPUT CARD> $TDDFT NSTATE=5 MULT=3 $END\n CHARGE OF MOLECULE        =    0\n"
              " SPIN MULTIPLICITY         =    1\n")
    firefly = (" INPUT CARD> $TDDFT NSTATE=5 MULT=3 $END\n CHARGE OF MOLECULE        =    0\n"
               "     MULT  =       1     ICHARG=       0     MAXIT =      30\n")
    assert outputs.charge_and_multiplicity(gamess, "gamout") == (0, 1)
    assert outputs.charge_and_multiplicity(firefly, "gamout") == (0, 1)
    assert outputs.charge_and_multiplicity(" charge of molecule      1\n state multiplicity      2\n",
                                           "gukout") == (1, 2)


def test_modes_printed_twice_are_listed_once():
    from miles.core.vibrations import VibrationalMode
    from miles.files.openbabel import _last_copy

    shapes = np.eye(9).reshape(9, 3, 3)
    harmonic = [VibrationalMode(f, 1.0, shapes[k]) for k, f in enumerate((1600.0, 3700.0, 3800.0))]
    assert _last_copy(harmonic + harmonic) == harmonic
    assert _last_copy(harmonic) == harmonic
    degenerate = [VibrationalMode(667.0, 1.0, shapes[0]), VibrationalMode(667.0, 1.0, shapes[1])]
    assert len(_last_copy(degenerate)) == 2   # same frequency, different motion: both kept
