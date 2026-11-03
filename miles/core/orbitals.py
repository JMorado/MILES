r"""Molecular orbitals: basis sets, coefficients and orbitals on a grid.

A quantum chemistry calculation describes its result as a basis set of
contracted Gaussian functions centred on the atoms, and each molecular
orbital as a vector of coefficients over that basis
(:mod:`miles.files.molden` reads them from Molden files). To draw an
orbital it is evaluated on a regular grid,

.. math::

   \psi_i(\mathbf r) = \sum_\mu C_{i\mu}\, \chi_\mu(\mathbf r),

and the grid is returned as a :class:`~miles.core.volume.Volume`.

Cartesian functions are each normalised (so ``xy`` carries a factor
``sqrt(3)`` relative to ``xx``), and spherical ones are the real solid
harmonics in the order ``m = 0, +1, -1, +2, -2, ...``, for d to h shells
(h shells are always spherical). Contracted functions are normalised
exactly. Lengths are in bohr internally; grids are returned in angstrom,
with values in bohr**-3/2.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
from typing import ClassVar

import numpy as np

from miles.core.units import BOHR
from miles.core.volume import Volume

SHELL_TYPES = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5}

CARTESIAN = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1),
        (1, 1, 1)],
    4: [(4, 0, 0), (0, 4, 0), (0, 0, 4), (3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3),
        (2, 2, 0), (2, 0, 2), (0, 2, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)],
}
"""Powers (a, b, c) of x, y, z for each cartesian function, in Molden order."""

Monomials = dict[tuple[int, int, int], float]


def _times(poly: Monomials, axis: int, factor: float = 1.0) -> Monomials:
    """Multiply a polynomial by x, y or z (axis 0, 1 or 2) and a factor."""
    out: Monomials = {}
    for powers, c in poly.items():
        raised = (powers[0] + (axis == 0), powers[1] + (axis == 1), powers[2] + (axis == 2))
        out[raised] = out.get(raised, 0.0) + factor * c
    return out


def _add(*polys: Monomials) -> Monomials:
    out: Monomials = {}
    for poly in polys:
        for powers, c in poly.items():
            out[powers] = out.get(powers, 0.0) + c
    return {powers: c for powers, c in out.items() if abs(c) > 1e-14}


def solid_harmonics(lmax: int) -> dict[int, list[Monomials]]:
    """Real regular solid harmonics as polynomials in x, y and z.

    The harmonics follow Racah's normalisation and are built with the
    standard recurrences (Helgaker, Jorgensen and Olsen, *Molecular
    Electronic-Structure Theory*, section 6.4). With a radial part
    normalised for ``x**l``, each function is normalised.

    Parameters
    ----------
    lmax : int
        Highest angular momentum.

    Returns
    -------
    dict of int to list of dict
        For each ``l`` from 2 to ``lmax``, the ``2l + 1`` functions in Molden
        order (m = 0, +1, -1, +2, -2, ...), each as a mapping from the powers
        ``(a, b, c)`` of ``x**a y**b z**c`` to its coefficient.
    """
    harmonics: dict[tuple[int, int], Monomials] = {(0, 0): {(0, 0, 0): 1.0}}
    for j in range(lmax):   # build angular momentum j + 1 from j and j - 1
        scale = math.sqrt((2.0 if j == 0 else 1.0) * (2 * j + 1) / (2 * j + 2))
        top, bottom = harmonics[(j, j)], harmonics[(j, -j)]
        harmonics[(j + 1, j + 1)] = _add(_times(top, 0, scale), _times(bottom, 1, -scale) if j else {})
        harmonics[(j + 1, -j - 1)] = _add(_times(top, 1, scale), _times(bottom, 0, scale) if j else {})
        for m in range(-j, j + 1):
            denominator = math.sqrt((j + m + 1) * (j - m + 1))
            rising = _times(harmonics[(j, m)], 2, (2 * j + 1) / denominator)
            previous = harmonics.get((j - 1, m), {})
            r_squared = _add(*(_times(_times(previous, d), d) for d in range(3)))
            lowering = math.sqrt((j + m) * (j - m)) / denominator
            harmonics[(j + 1, m)] = _add(rising, {k: -lowering * c for k, c in r_squared.items()})
    tables = {}
    for angmom in range(2, lmax + 1):
        order = [0, *(sign * m for m in range(1, angmom + 1) for sign in (1, -1))]
        tables[angmom] = [harmonics[(angmom, m)] for m in order]
    return tables


SPHERICAL = solid_harmonics(5)
"""Real solid harmonics (m = 0, +1, -1, ...) for d to h shells, for a radial part normalised for x**l."""

CUTOFF = 1e-8
"""Each primitive is evaluated only where it exceeds this fraction of its peak."""

OCCUPIED = 1e-6
"""Orbitals with more electrons than this count as occupied."""

DENSITY_NAME = "Electron density"
"""Name of the grid made by :meth:`OrbitalSet.density`."""

SPIN_NAME = "Spin density"
"""Name of the grid made by :meth:`OrbitalSet.spin_density`."""

SCREEN = 1e-14
"""Pairs of shells whose most diffuse primitives overlap less than this are left out of the overlap matrix."""


class Spin(StrEnum):
    """Spin of an orbital; NONE for restricted calculations."""

    ALPHA = "alpha"
    BETA = "beta"
    NONE = ""


def double_factorial(n: int) -> int:
    """Return ``n!! = n (n - 2) (n - 4) ...``, and 1 for ``n <= 0``."""
    return math.prod(range(n, 0, -2)) if n > 0 else 1


def primitive_norm(exponents: np.ndarray, angmom: int) -> np.ndarray:
    """Normalisation constants of ``x**l exp(-alpha r**2)``."""
    return ((2 * exponents / math.pi) ** 0.75 * (4 * exponents) ** (angmom / 2)
            / math.sqrt(double_factorial(2 * angmom - 1)))


def self_overlap(exponents: np.ndarray, weights: np.ndarray, angmom: int) -> float:
    """Squared norm of ``x**l sum_i w_i exp(-alpha_i r**2)``."""
    total = exponents[:, None] + exponents[None, :]
    overlap = (math.pi / total) ** 1.5 * double_factorial(2 * angmom - 1) / (2 * total) ** angmom
    return float(weights @ overlap @ weights)


def cartesian_factor(a: int, b: int, c: int) -> float:
    """Factor that normalises ``x**a y**b z**c`` when the radial part is normalised for ``x**l``."""
    parts = double_factorial(2 * a - 1) * double_factorial(2 * b - 1) * double_factorial(2 * c - 1)
    return math.sqrt(double_factorial(2 * (a + b + c) - 1) / parts)


@dataclass
class Shell:
    """A contracted shell of Gaussian functions on one atom.

    Attributes
    ----------
    center : numpy.ndarray
        Position in bohr, shape ``(3,)``.
    angmom : int
        Angular momentum: 0 for s, 1 for p, and so on.
    exponents : numpy.ndarray
        Primitive exponents in bohr**-2.
    coefficients : numpy.ndarray
        Weights of the bare primitives ``exp(-alpha r**2)``, scaled so that the
        ``x**l`` function of the shell is normalised.
    spherical : bool
        Real solid harmonics (``2l + 1`` functions) instead of cartesian ones.
    """

    center: np.ndarray
    angmom: int
    exponents: np.ndarray
    coefficients: np.ndarray
    spherical: bool = False

    @property
    def size(self) -> int:
        """Number of basis functions in the shell."""
        if self.spherical and self.angmom >= 2:
            return 2 * self.angmom + 1
        return len(CARTESIAN[self.angmom])

    def functions(self) -> list[dict[tuple[int, int, int], float]]:
        """Angular part of each basis function as weighted monomials ``x**a y**b z**c``."""
        if self.spherical and self.angmom >= 2:
            return SPHERICAL[self.angmom]
        return [{powers: cartesian_factor(*powers)} for powers in CARTESIAN[self.angmom]]


def contracted_shell(center: np.ndarray, angmom: int, exponents: np.ndarray, coefficients: np.ndarray,
                     spherical: bool = False, normalized_primitives: bool = True) -> Shell:
    """Build a shell whose basis functions are normalised.

    Parameters
    ----------
    center : numpy.ndarray
        Position in bohr.
    angmom : int
        Angular momentum (0 to 4; 5 for spherical functions).
    exponents : numpy.ndarray
        Primitive exponents in bohr**-2.
    coefficients : numpy.ndarray
        Contraction coefficients; their overall scale does not matter.
    spherical : bool, optional
        Real solid harmonics instead of cartesian functions (d and higher).
    normalized_primitives : bool, optional
        Whether the coefficients multiply normalised primitives, as the Molden
        format specifies, rather than bare ``x**l exp(-alpha r**2)``.

    Returns
    -------
    Shell
        The normalised shell.

    Raises
    ------
    ValueError
        If the angular momentum is not supported or the contraction is zero.
    """
    if angmom not in CARTESIAN and not (spherical and angmom in SPHERICAL):
        raise ValueError(f"Angular momentum {angmom} is not supported")
    exponents = np.asarray(exponents, dtype=float)
    weights = np.asarray(coefficients, dtype=float)
    if normalized_primitives:
        weights = weights * primitive_norm(exponents, angmom)
    norm = self_overlap(exponents, weights, angmom)
    if norm <= 0.0:
        raise ValueError("Contraction with zero norm")
    return Shell(np.asarray(center, dtype=float), angmom, exponents, weights / math.sqrt(norm), spherical)


@dataclass
class Orbital:
    """One molecular orbital.

    Attributes
    ----------
    energy : float
        Orbital energy in hartree.
    occupation : float
        Number of electrons in it.
    spin : Spin
        Alpha or beta for unrestricted calculations, otherwise NONE.
    symmetry : str
        Symmetry label from the file (may be empty).
    """

    energy: float
    occupation: float
    spin: Spin = Spin.NONE
    symmetry: str = ""


@dataclass
class OrbitalSet:
    """Molecular orbitals with the basis set needed to evaluate them.

    The class constants :attr:`OCCUPIED`, :attr:`DENSITY_NAME` and
    :attr:`SPIN_NAME` repeat the module's, for code that only has the class.

    Attributes
    ----------
    orbitals : list of Orbital
        Energies, occupations, spins and symmetry labels, in file order.
    shells : list of Shell
        The basis set, in the order of the coefficients.
    coefficients : numpy.ndarray
        Shape ``(n_orbitals, n_basis_functions)``.
    atoms : numpy.ndarray
        Atom positions in bohr, shape ``(n_atoms, 3)``; they set the grid extent.
    """

    OCCUPIED: ClassVar[float] = OCCUPIED
    DENSITY_NAME: ClassVar[str] = DENSITY_NAME
    SPIN_NAME: ClassVar[str] = SPIN_NAME

    orbitals: list[Orbital]
    shells: list[Shell]
    coefficients: np.ndarray
    atoms: np.ndarray

    @property
    def unrestricted(self) -> bool:
        """Whether alpha and beta orbitals are listed separately."""
        return any(o.spin == Spin.BETA for o in self.orbitals)

    def _channel(self, spin: Spin) -> list[int]:
        """Orbitals of one spin, or all of them for a restricted calculation, by index."""
        if not self.unrestricted:
            return list(range(len(self.orbitals)))
        return [k for k, o in enumerate(self.orbitals) if (o.spin or Spin.ALPHA) == spin]

    def homo(self) -> int:
        """Index of the highest occupied orbital (alpha for unrestricted calculations).

        Returns
        -------
        int
            Index into :attr:`orbitals`; 0 when no orbital is occupied.
        """
        occupied = [k for k in self._channel(Spin.ALPHA) if self.orbitals[k].occupation > OCCUPIED]
        return max(occupied, key=lambda k: self.orbitals[k].energy) if occupied else 0

    def label(self, index: int) -> str:
        """Name of an orbital relative to the frontier orbitals.

        Parameters
        ----------
        index : int
            Index into :attr:`orbitals`.

        Returns
        -------
        str
            ``"HOMO"``, ``"HOMO-2"``, ``"LUMO+1"`` and so on, with ``"alpha "``
            or ``"beta "`` in front for unrestricted calculations; empty for an
            index out of range.
        """
        if not 0 <= index < len(self.orbitals):
            return ""
        orbital = self.orbitals[index]
        channel = self._channel(orbital.spin or Spin.ALPHA)
        if orbital.occupation > OCCUPIED:
            above = [k for k in channel if self.orbitals[k].occupation > OCCUPIED
                     and (self.orbitals[k].energy, k) > (orbital.energy, index)]
            text = f"HOMO-{len(above)}" if above else "HOMO"
        else:
            below = [k for k in channel if self.orbitals[k].occupation <= OCCUPIED
                     and (self.orbitals[k].energy, k) < (orbital.energy, index)]
            text = f"LUMO+{len(below)}" if below else "LUMO"
        return f"{orbital.spin or Spin.ALPHA} {text}" if self.unrestricted else text

    def name(self, index: int) -> str:
        """Name of the grid of an orbital.

        Parameters
        ----------
        index : int
            Index into :attr:`orbitals`.

        Returns
        -------
        str
            ``"MO 21 (HOMO)"``, counting from 1 in file order, or ``"MO 21"``
            when the orbital has no label.
        """
        label = self.label(index)
        return f"MO {index + 1} ({label})" if label else f"MO {index + 1}"

    def evaluate(self, index: int, spacing: float = 0.25, margin: float = 4.0) -> Volume:
        """Values of one orbital on a grid around the molecule.

        Parameters
        ----------
        index : int
            Index into :attr:`orbitals`.
        spacing : float, optional
            Grid spacing in angstrom.
        margin : float, optional
            Distance in angstrom from the outermost atoms to the grid edge.

        Returns
        -------
        Volume
            Grid in angstrom named by :meth:`name`, e.g. ``"MO 21 (HOMO)"``.

        Raises
        ------
        IndexError
            If there is no orbital ``index``.
        """
        if not 0 <= index < len(self.orbitals):
            raise IndexError(f"No orbital {index}")
        low, axes = self._grid(spacing, margin)
        values = self._values(index, axes)
        return Volume(self.name(index), low * BOHR, np.eye(3) * spacing, values.astype(np.float32))

    def density(self, spacing: float = 0.25, margin: float = 4.0) -> Volume:
        """Electron density on a grid around the molecule.

        The density is the sum of the squares of the occupied orbitals, each
        weighted by its occupation, in electrons per bohr**3. Grids are those
        of :meth:`evaluate`; tight core functions are only sampled at the grid
        points, so integrating the grid gives the number of electrons only
        roughly, while the valence density that isosurfaces show is accurate.

        Parameters
        ----------
        spacing : float, optional
            Grid spacing in angstrom.
        margin : float, optional
            Distance in angstrom from the outermost atoms to the grid edge.

        Returns
        -------
        Volume
            Grid named :data:`DENSITY_NAME`.
        """
        low, axes = self._grid(spacing, margin)
        values = self._occupied_sum(axes, lambda orbital: orbital.occupation)
        return Volume(DENSITY_NAME, low * BOHR, np.eye(3) * spacing, values.astype(np.float32))

    def spin_density(self, spacing: float = 0.25, margin: float = 4.0) -> Volume:
        """Spin density (alpha minus beta electron density) of an unrestricted calculation.

        Parameters
        ----------
        spacing : float, optional
            Grid spacing in angstrom.
        margin : float, optional
            Distance in angstrom from the outermost atoms to the grid edge.

        Returns
        -------
        Volume
            Grid named :data:`SPIN_NAME`, in electrons per bohr**3; it
            integrates to the number of unpaired electrons.

        Raises
        ------
        ValueError
            If the orbitals are not of an unrestricted calculation.
        """
        if not self.unrestricted:
            raise ValueError("A spin density needs separate alpha and beta orbitals")
        low, axes = self._grid(spacing, margin)
        values = self._occupied_sum(
            axes, lambda orbital: orbital.occupation * (-1.0 if orbital.spin == Spin.BETA else 1.0))
        return Volume(SPIN_NAME, low * BOHR, np.eye(3) * spacing, values.astype(np.float32), both_signs=True)

    def _grid(self, spacing: float, margin: float) -> tuple[np.ndarray, list[np.ndarray]]:
        """Corner (bohr) and axis coordinates (bohr) of a grid around the atoms."""
        step = spacing / BOHR
        low = self.atoms.min(axis=0) - margin / BOHR
        counts = np.ceil((self.atoms.max(axis=0) + margin / BOHR - low) / step).astype(int) + 1
        return low, [low[d] + step * np.arange(counts[d]) for d in range(3)]

    def _values(self, index: int, axes: list[np.ndarray]) -> np.ndarray:
        """Values of orbital ``index`` at the points of a grid."""
        values = np.zeros(tuple(len(a) for a in axes))
        offset = 0
        for shell in self.shells:
            weights = self.coefficients[index, offset:offset + shell.size]
            offset += shell.size
            if np.abs(weights).max() > 1e-10:
                _add_shell(values, axes, shell, weights)
        return values

    def _occupied_sum(self, axes: list[np.ndarray], weight: Callable[[Orbital], float]) -> np.ndarray:
        """Sum of ``weight(orbital) * psi**2`` over the occupied orbitals."""
        total = np.zeros(tuple(len(a) for a in axes))
        for k, orbital in enumerate(self.orbitals):
            if orbital.occupation > OCCUPIED:
                total += weight(orbital) * self._values(k, axes) ** 2
        return total


def _span(coords: np.ndarray, reach: float) -> slice:
    """Slice of the sorted values within ``reach`` of zero."""
    return slice(int(np.searchsorted(coords, -reach)), int(np.searchsorted(coords, reach, side="right")))


def _add_shell(values: np.ndarray, axes: list[np.ndarray], shell: Shell, weights: np.ndarray) -> None:
    """Add ``sum_m weights[m] * chi_m`` of one shell to the grid where the shell is not negligible.

    Gaussians factorise along x, y and z, so every primitive is an outer
    product of three short vectors over a box sized to its own exponent.
    """
    reach = math.sqrt(-math.log(CUTOFF) / shell.exponents.min()) + 1.0
    window = tuple(_span(axes[d] - shell.center[d], reach) for d in range(3))
    x, y, z = (axes[d][window[d]] - shell.center[d] for d in range(3))
    if not (x.size and y.size and z.size):
        return
    radial = np.zeros((x.size, y.size, z.size))
    for alpha, c in zip(shell.exponents, shell.coefficients, strict=True):
        r = math.sqrt(-math.log(CUTOFF) / alpha) + 1.0
        sx, sy, sz = _span(x, r), _span(y, r), _span(z, r)
        gx, gy, gz = np.exp(-alpha * x[sx] ** 2), np.exp(-alpha * y[sy] ** 2), np.exp(-alpha * z[sz] ** 2)
        radial[sx, sy, sz] += c * gx[:, None, None] * (gy[:, None] * gz[None, :])[None, :, :]
    monomials: dict[tuple[int, int, int], float] = {}
    for weight, function in zip(weights, shell.functions(), strict=True):
        for powers, factor in function.items():
            monomials[powers] = monomials.get(powers, 0.0) + float(weight) * factor
    angular = np.zeros_like(radial)
    for (a, b, c), factor in monomials.items():
        if factor != 0.0:
            angular += factor * (x**a)[:, None, None] * ((y**b)[:, None] * (z**c)[None, :])[None, :, :]
    values[window] += angular * radial


# -- overlap of the basis functions ----------------------------------------------------------------

def _monomials(angmom: int) -> list[tuple[int, int, int]]:
    return [(a, b, angmom - a - b) for a in range(angmom, -1, -1) for b in range(angmom - a, -1, -1)]


def _angular_matrix(shell: Shell) -> np.ndarray:
    """Coefficients of the shell's functions (rows) over the monomials of its degree (columns)."""
    columns = {powers: k for k, powers in enumerate(_monomials(shell.angmom))}
    matrix = np.zeros((shell.size, len(columns)))
    for row, function in enumerate(shell.functions()):
        for powers, factor in function.items():
            matrix[row, columns[powers]] += factor
    return matrix


def _overlaps_1d(la: int, lb: int, alpha: np.ndarray, beta: np.ndarray,
                 a: np.ndarray, b: np.ndarray) -> list[list[np.ndarray]]:
    """One-dimensional overlaps of ``x**i exp(-alpha x**2)`` at ``a`` with ``x**j exp(-beta x**2)`` at ``b``.

    ``i`` runs up to ``la`` and ``j`` up to ``lb``.

    Obara-Saika recurrence, for many pairs of primitives at once (one per
    element of the arrays).
    """
    p = alpha + beta
    centre = (alpha * a + beta * b) / p
    pa, pb, half = centre - a, centre - b, 0.5 / p
    s = [[np.zeros(0)] * (lb + 1) for _ in range(la + 1)]
    s[0][0] = np.sqrt(math.pi / p) * np.exp(-alpha * beta / p * (a - b) ** 2)
    for j in range(1, lb + 1):
        s[0][j] = pb * s[0][j - 1] + (j - 1) * half * s[0][j - 2] if j > 1 else pb * s[0][0]
    for i in range(1, la + 1):
        for j in range(lb + 1):
            value = pa * s[i - 1][j]
            if i > 1:
                value = value + (i - 1) * half * s[i - 2][j]
            if j > 0:
                value = value + j * half * s[i - 1][j - 1]
            s[i][j] = value
    return s


def overlap_matrix(shells: list[Shell]) -> np.ndarray:
    """Overlap of every pair of basis functions, computed exactly.

    Parameters
    ----------
    shells : list of Shell
        The basis set.

    Returns
    -------
    numpy.ndarray
        Symmetric matrix ``S[mu, nu] = <chi_mu | chi_nu>``, with ones on the
        diagonal for the normalised functions of :func:`contracted_shell`.
    """
    sizes = [shell.size for shell in shells]
    offsets = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
    overlap = np.zeros((offsets[-1], offsets[-1]))
    kinds: dict[tuple[int, bool], list[int]] = {}
    for k, shell in enumerate(shells):
        kinds.setdefault((shell.angmom, shell.spherical), []).append(k)
    primitives = {kind: _primitives(shells, members) for kind, members in kinds.items()}
    keys = sorted(kinds)
    for ka, kind_a in enumerate(keys):
        for kind_b in keys[ka:]:
            block = _overlap_block(primitives[kind_a], primitives[kind_b], kind_a[0], kind_b[0],
                                   kind_a == kind_b)
            if block is None:
                continue
            pairs, cartesian = block
            first = offsets[np.array(kinds[kind_a])[pairs[:, 0]]]
            second = offsets[np.array(kinds[kind_b])[pairs[:, 1]]]
            angular_a = _angular_matrix(shells[kinds[kind_a][0]])
            angular_b = _angular_matrix(shells[kinds[kind_b][0]])
            blocks = np.einsum("am,pmn,bn->pab", angular_a, cartesian, angular_b)
            rows = first[:, None, None] + np.arange(angular_a.shape[0])[None, :, None]
            columns = second[:, None, None] + np.arange(angular_b.shape[0])[None, None, :]
            overlap[rows, columns] = blocks
            overlap[columns, rows] = blocks
    return overlap


def _primitives(shells: list[Shell], members: list[int]) -> tuple[np.ndarray, ...]:
    """Exponent, weight, centre and shell (position in ``members``) of every primitive of some shells."""
    counts = [len(shells[k].exponents) for k in members]
    exponents = np.concatenate([shells[k].exponents for k in members])
    weights = np.concatenate([shells[k].coefficients for k in members])
    centres = np.repeat(np.array([shells[k].center for k in members]), counts, axis=0)
    owners = np.repeat(np.arange(len(members)), counts)
    return exponents, weights, centres, owners


def _overlap_block(prims_a: tuple[np.ndarray, ...], prims_b: tuple[np.ndarray, ...], la: int, lb: int,
                   same: bool) -> tuple[np.ndarray, np.ndarray] | None:
    """Overlaps of the cartesian monomials of every pair of shells of two kinds.

    Returns the pairs of shells (positions in their kinds), shape
    ``(n_pairs, 2)``, and their overlaps over the monomials of degree ``la``
    and ``lb``, shape ``(n_pairs, n_monomials_a, n_monomials_b)``; None when
    no pair overlaps.
    """
    exp_a, weight_a, centre_a, owner_a = prims_a
    exp_b, weight_b, centre_b, owner_b = prims_b
    grid = np.meshgrid(np.arange(exp_a.size), np.arange(exp_b.size), indexing="ij")
    ia, ib = grid[0].ravel(), grid[1].ravel()
    alpha, beta = exp_a[ia], exp_b[ib]
    distance = np.sum((centre_a[ia] - centre_b[ib]) ** 2, axis=1)
    keep = alpha * beta / (alpha + beta) * distance < -math.log(SCREEN)
    if same:
        keep &= owner_a[ia] <= owner_b[ib]
    if not keep.any():
        return None
    ia, ib, alpha, beta = ia[keep], ib[keep], alpha[keep], beta[keep]
    pair_ids = owner_a[ia] * (int(owner_b.max()) + 1) + owner_b[ib]
    unique, owner = np.unique(pair_ids, return_inverse=True)
    pairs = np.stack([unique // (int(owner_b.max()) + 1), unique % (int(owner_b.max()) + 1)], axis=1)
    weight = weight_a[ia] * weight_b[ib]
    axes = [_overlaps_1d(la, lb, alpha, beta, centre_a[ia, d], centre_b[ib, d]) for d in range(3)]
    mono_a, mono_b = _monomials(la), _monomials(lb)
    cartesian = np.zeros((len(unique), len(mono_a), len(mono_b)))
    for p, (xa, ya, za) in enumerate(mono_a):
        for q, (xb, yb, zb) in enumerate(mono_b):
            values = weight * axes[0][xa][xb] * axes[1][ya][yb] * axes[2][za][zb]
            cartesian[:, p, q] = np.bincount(owner, weights=values, minlength=len(unique))
    return pairs, cartesian
