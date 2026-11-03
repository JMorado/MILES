"""Write the quantum chemistry example files in ``examples/qm``.

The files show what MILES can open from calculations: grids in cube
files, molecular orbitals and vibrations in Molden files. They are made
with PySCF, which MILES itself does not use, so the script runs in an
environment of its own:

    micromamba create -n miles-examples -c conda-forge python=3.12 pyscf geometric numpy
    micromamba run -n miles-examples python examples/make_qm_examples.py

Every molecule is optimised with B3LYP/def2-SVP (geomeTRIC), then:

* chloromethane: electron density and electrostatic potential on the same grid;
* benzene: the HOMO on a grid;
* formaldehyde: orbitals, harmonic frequencies and infrared intensities in one
  Molden file (intensities from finite differences of the dipole moment);
* nitrogen dioxide: orbitals of an unrestricted (open-shell) calculation.

It takes a few minutes and gives the same files every time it is run with
the same versions of PySCF and geomeTRIC.
"""

import contextlib
import os
import sys
import tempfile
from collections.abc import Callable, Iterator
from pathlib import Path

import numpy as np
from pyscf import dft, gto, lib
from pyscf.geomopt.geometric_solver import optimize
from pyscf.hessian import thermo
from pyscf.tools import cubegen, molden

OUT = Path(__file__).resolve().with_name("qm")
FUNCTIONAL = "b3lyp"
BASIS = "def2-svp"
METHOD = "B3LYP/def2-SVP"
GRID_POINTS = 48       # per axis; about 1.5 MB per cube file
DIPOLE_STEP = 0.005    # angstrom, for the dipole derivatives
IR_FACTOR = 42.2561    # km/mol per (debye/angstrom)^2/amu
BOHR = 0.529177210903  # angstrom

MOLECULES = {
    "chloromethane": ("""C 0.000 0.000 0.000; Cl 0.000 0.000 1.780; H 1.027 0.000 -0.363
                         H -0.513 0.889 -0.363; H -0.513 -0.889 -0.363""", 0),
    "benzene": ("; ".join(f"C {1.39 * np.cos(k * np.pi / 3):.4f} {1.39 * np.sin(k * np.pi / 3):.4f} 0.0"
                          for k in range(6)) + "; " +
                "; ".join(f"H {2.47 * np.cos(k * np.pi / 3):.4f} {2.47 * np.sin(k * np.pi / 3):.4f} 0.0"
                          for k in range(6)), 0),
    "formaldehyde": ("C 0.000 0.000 0.000; O 0.000 0.000 1.205; H 0.000 0.940 -0.587; "
                     "H 0.000 -0.940 -0.587", 0),
    "nitrogen_dioxide": ("N 0.000 0.000 0.000; O 0.000 1.100 0.470; O 0.000 -1.100 0.470", 1),
}


def molecule(name: str) -> gto.Mole:
    """Build the starting geometry of one of the :data:`MOLECULES` (angstrom).

    Parameters
    ----------
    name : str
        Key of :data:`MOLECULES`.

    Returns
    -------
    pyscf.gto.Mole
        The molecule with the basis set; ``spin`` is the number of unpaired electrons.
    """
    atoms, unpaired = MOLECULES[name]
    return gto.M(atom=atoms, basis=BASIS, spin=unpaired, unit="Angstrom", verbose=0)


def solve(mol: gto.Mole) -> dft.rks.RKS:
    """Converged Kohn-Sham calculation, restricted for closed shells and unrestricted otherwise.

    Parameters
    ----------
    mol : pyscf.gto.Mole
        The molecule.

    Returns
    -------
    pyscf.dft.rks.RKS or pyscf.dft.uks.UKS
        The converged calculation.
    """
    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.xc = FUNCTIONAL
    mf.conv_tol = 1e-10
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {mol.atom}")
    return mf


@contextlib.contextmanager
def quiet() -> Iterator[None]:
    """Send everything written to standard output and error, also by libraries, to the null device."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    with open(os.devnull, "w") as null:
        os.dup2(null.fileno(), 1)
        os.dup2(null.fileno(), 2)
        try:
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            for fd in saved:
                os.close(fd)


def optimised(mol: gto.Mole) -> gto.Mole:
    """Geometry optimised with geomeTRIC; its log is discarded and its files go to a temporary folder.

    Parameters
    ----------
    mol : pyscf.gto.Mole
        Starting geometry.

    Returns
    -------
    pyscf.gto.Mole
        The optimised molecule.
    """
    with tempfile.TemporaryDirectory() as scratch, quiet():
        previous, lib.param.TMPDIR = lib.param.TMPDIR, scratch
        cwd = os.getcwd()
        os.chdir(scratch)
        try:
            mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
            mf.xc = FUNCTIONAL
            mf.verbose = 0
            result = optimize(mf, maxsteps=100)
        finally:
            os.chdir(cwd)
            lib.param.TMPDIR = previous
    result.verbose = 0
    return result


def retitle(path: Path, title: str, description: str) -> None:
    """Replace the two comment lines PySCF writes at the top of a cube file.

    The second line says what the grid holds; MILES names the grid
    after it. PySCF puts a date there, which would make every run differ.

    Parameters
    ----------
    path : pathlib.Path
        Cube file to change in place.
    title : str
        First line, the job.
    description : str
        Second line, the quantity on the grid.
    """
    lines = path.read_text().splitlines(keepends=True)
    lines[0], lines[1] = title + "\n", description + "\n"
    path.write_text("".join(lines))


def write_cube(path: Path, writer: Callable[..., np.ndarray], mol: gto.Mole, data: np.ndarray,
               title: str, description: str) -> None:
    """Write one PySCF cube file with fixed comment lines.

    Parameters
    ----------
    path : pathlib.Path
        Output file.
    writer : callable
        ``cubegen.density``, ``cubegen.mep`` or ``cubegen.orbital``.
    mol : pyscf.gto.Mole
        The molecule.
    data : numpy.ndarray
        Density matrix, or orbital coefficients for ``cubegen.orbital``.
    title, description : str
        Comment lines, see :func:`retitle`.
    """
    writer(mol, str(path), data, nx=GRID_POINTS, ny=GRID_POINTS, nz=GRID_POINTS)
    retitle(path, title, description)


def dipole(mol: gto.Mole, coords: np.ndarray) -> np.ndarray:
    """Dipole moment in debye of the molecule with its atoms at ``coords`` (angstrom)."""
    moved = mol.set_geom_(coords, unit="Angstrom", inplace=False)
    return np.asarray(solve(moved).dip_moment(unit="Debye", verbose=0))


def infrared(mol: gto.Mole, modes: np.ndarray) -> np.ndarray:
    """Infrared intensities from central differences of the dipole moment.

    Parameters
    ----------
    mol : pyscf.gto.Mole
        Molecule at its optimised geometry.
    modes : numpy.ndarray
        Normal modes as Cartesian displacements per unit mass-weighted
        coordinate (amu**-1/2), shape ``(n_modes, n_atoms, 3)``.

    Returns
    -------
    numpy.ndarray
        Intensities in km/mol, one per mode.
    """
    coords = mol.atom_coords(unit="Angstrom")
    derivative = np.zeros((mol.natm, 3, 3))   # d mu / d x, debye per angstrom
    for atom in range(mol.natm):
        for axis in range(3):
            forward, backward = coords.copy(), coords.copy()
            forward[atom, axis] += DIPOLE_STEP
            backward[atom, axis] -= DIPOLE_STEP
            derivative[atom, axis] = (dipole(mol, forward) - dipole(mol, backward)) / (2 * DIPOLE_STEP)
    change = np.einsum("iac,kia->kc", derivative, modes)
    return IR_FACTOR * (change ** 2).sum(axis=1)


def vibration_sections(mol: gto.Mole, frequencies: np.ndarray, intensities: np.ndarray,
                       modes: np.ndarray) -> str:
    """Molden ``[FREQ]``, ``[INT]``, ``[FR-COORD]`` and ``[FR-NORM-COORD]`` sections.

    ``[INT]`` comes straight after ``[FREQ]``: Open Babel reads intensities
    only when they come before the coordinates.

    Parameters
    ----------
    mol : pyscf.gto.Mole
        Molecule at the geometry of the frequency calculation.
    frequencies : numpy.ndarray
        Wavenumbers in cm**-1.
    intensities : numpy.ndarray
        Infrared intensities in km/mol.
    modes : numpy.ndarray
        Displacements, shape ``(n_modes, n_atoms, 3)``.

    Returns
    -------
    str
        The sections, ending with a newline.
    """
    lines = ["[FREQ]"] + [f"{f:12.4f}" for f in frequencies]
    lines += ["[INT]"] + [f"{i:12.4f}" for i in intensities]
    lines += ["[FR-COORD]"]
    lines += [f"{mol.atom_pure_symbol(k):2s} {x:14.8f} {y:14.8f} {z:14.8f}"
              for k, (x, y, z) in enumerate(mol.atom_coords(unit="Bohr"))]
    lines += ["[FR-NORM-COORD]"]
    for n, mode in enumerate(modes, start=1):
        lines.append(f"vibration {n}")
        lines += [f"{x:14.8f} {y:14.8f} {z:14.8f}" for x, y, z in mode]
    return "\n".join(lines) + "\n"


def chloromethane() -> None:
    """Density and electrostatic potential cubes on the same grid."""
    mf = solve(optimised(molecule("chloromethane")))
    dm = mf.make_rdm1()
    title = f"Chloromethane, {METHOD}"
    write_cube(OUT / "chloromethane_density.cube", cubegen.density, mf.mol, dm, title,
               "Electron density (e/bohr^3)")
    write_cube(OUT / "chloromethane_potential.cube", cubegen.mep, mf.mol, dm, title,
               "Electrostatic potential (hartree/e)")


def benzene() -> None:
    """Write the highest occupied orbital of benzene on a grid."""
    mf = solve(optimised(molecule("benzene")))
    homo = mf.mol.nelectron // 2 - 1
    write_cube(OUT / "benzene_homo.cube", cubegen.orbital, mf.mol, mf.mo_coeff[:, homo],
               f"Benzene, {METHOD}", f"HOMO (MO {homo + 1}), one of a degenerate pair")


def formaldehyde() -> None:
    """Orbitals, frequencies and infrared intensities in one Molden file."""
    mf = solve(optimised(molecule("formaldehyde")))
    analysis = thermo.harmonic_analysis(mf.mol, mf.Hessian().kernel())
    wavenumbers = np.asarray(analysis["freq_wavenumber"])
    frequencies = np.where(np.iscomplex(wavenumbers), -np.abs(wavenumbers.imag), wavenumbers.real)
    modes = np.asarray(analysis["norm_mode"]).real
    path = OUT / "formaldehyde.molden"
    molden.from_scf(mf, str(path))
    with open(path, "a") as f:
        f.write(vibration_sections(mf.mol, frequencies, infrared(mf.mol, modes), modes))


def nitrogen_dioxide() -> None:
    """Orbitals of an open-shell molecule, alpha and beta listed separately."""
    mf = solve(optimised(molecule("nitrogen_dioxide")))
    molden.from_scf(mf, str(OUT / "nitrogen_dioxide.molden"))


def main() -> None:
    """Make every example file."""
    lib.num_threads(1)   # parallel sums round differently from run to run
    OUT.mkdir(exist_ok=True)
    for make in (chloromethane, benzene, formaldehyde, nitrogen_dioxide):
        make()
        print(f"{make.__name__}: done")


if __name__ == "__main__":
    main()
