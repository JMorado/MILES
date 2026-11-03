"""The molecule model: a table of atom columns, a table of bonds and pinned measurements.

Atoms are the rows of a column table. Every per-atom quantity (atomic
number, name, residue, coordinates, and any property a reader or plugin
adds) is one NumPy array, so whole-structure work is array arithmetic
instead of a loop over Python objects. Bonds are a second table of index
pairs and orders.

Storage is shared between a molecule and its copies until one of them
writes: :meth:`Molecule.copy` costs the same for ten atoms as for a million,
and an edit copies only the columns it changes. That makes undo snapshots
cheap. Arrays handed out by a molecule are read-only views; change a
molecule through its methods (or by assigning a whole column, as in
``mol.coords = new``), never by writing into an array it returned.

:attr:`Molecule.atoms` gives per-atom access for readable code
(``mol.atoms[i].symbol``); bulk work should read columns with
:meth:`Molecule.column` and :attr:`Molecule.coords`.
"""

import math
import sys
import threading
from collections import deque
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Any, Self, SupportsIndex, overload

import numpy as np
from numpy.typing import ArrayLike, DTypeLike

from miles.core.elements import ELEMENTS, MASSES, SYMBOLS

Bond = tuple[int, int]


class Representation(StrEnum):
    """How an atom is drawn. ``DEFAULT`` means "use the view's style".

    Which atom is drawn how is view state, kept by the document
    (:class:`miles.core.display.AtomDisplay`), not by the molecule.
    """

    DEFAULT = ""
    BALL = "ball"
    TUBE = "tube"
    WIRE = "wire"
    SPACEFILL = "spacefill"
    CARTOON = "cartoon"
    HIDDEN = "hidden"


@dataclass(slots=True)
class Atom:
    """A free-standing atom record, used to add atoms and to build molecules.

    Inside a molecule an atom is a row of columns, read and changed through
    :class:`AtomView`; this record carries the same fields between molecules.

    Attributes
    ----------
    symbol : str
        Element symbol (``"X"`` for a dummy atom).
    geom : str or None
        Ideal arrangement of the atom's bonds, one of ``geometry.TEMPLATES``.
    nsites : int or None
        Number of bonds the atom should carry.
    valence : float or None
        Bond-order sum the atom should reach.
    inferred : bool
        True when the three fields above were guessed from a file rather than
        set by the builder; Add Hydrogens then falls back to standard valences.
    formal_charge : int or None
        Formal charge as given by a file or SMILES (``-1`` for a carboxylate
        oxygen); None when nothing gave it, in which case it is worked out
        from the bonds when needed (see :func:`miles.core.charges.formal_charges`).
    isotope : int
        Mass number of a specific isotope (13 for carbon-13); 0 for the
        natural mixture.
    name, resname, resseq, chain, hetero : str, str, int, str, bool
        Residue data from PDB-like files (atom name, residue name and number,
        chain identifier, HETATM flag).
    icode : str
        Insertion code of the residue (``"A"`` in residue ``52A``), ``""`` for none.
    altloc : str
        Alternate location indicator of the atom (``"A"``, ``"B"``), ``""`` for none.
    ss : str
        Secondary structure of the residue: ``"H"`` helix, ``"E"`` strand or ``""``.
    """

    symbol: str
    geom: str | None = None
    nsites: int | None = None
    valence: float | None = None
    inferred: bool = False
    formal_charge: int | None = None
    isotope: int = 0
    name: str = ""
    altloc: str = ""
    resname: str = ""
    resseq: int = 0
    icode: str = ""
    chain: str = ""
    hetero: bool = False
    ss: str = ""


def hydrogen_atom() -> Atom:
    """Return a hydrogen with one single bond.

    Returns
    -------
    Atom
        Fresh hydrogen record.
    """
    return Atom("H", "linear", 1, 1)


# -- the columns ----------------------------------------------------------------------------

_SYMBOLS = SYMBOLS
_NUMBERS = {e.symbol: e.number for e in ELEMENTS.values()}
_MASSES = MASSES


@dataclass(frozen=True, slots=True)
class ColumnSpec:
    """Type and default value of an atom column."""

    dtype: np.dtype
    default: Any
    width: int = 1   # values per atom: 3 for coordinates


def _spec(dtype: DTypeLike, default: Any, width: int = 1) -> ColumnSpec:
    return ColumnSpec(np.dtype(dtype), default, width)


UNSET_CHARGE = -128
"""Value of the ``formal_charge`` column for atoms whose formal charge nothing gave."""

ATOM_COLUMNS: dict[str, ColumnSpec] = {
    "number": _spec(np.uint8, 0),           # atomic number; 0 is a dummy atom
    "geom": _spec(object, None),
    "nsites": _spec(np.int16, -1),          # -1: not set
    "valence": _spec(np.float64, math.nan),  # NaN: not set
    "inferred": _spec(np.bool_, False),
    "formal_charge": _spec(np.int8, UNSET_CHARGE),
    "isotope": _spec(np.uint16, 0),         # mass number; 0: the natural mixture
    "name": _spec(object, ""),
    "altloc": _spec(object, ""),
    "resname": _spec(object, ""),
    "resseq": _spec(np.int64, 0),
    "icode": _spec(object, ""),
    "chain": _spec(object, ""),
    "hetero": _spec(np.bool_, False),
    "ss": _spec(object, ""),
    "uid": _spec(np.int64, 0),
    "coords": _spec(np.float64, 0.0, 3),
}
"""The columns every molecule has, with their types and defaults."""

RECORD_FIELDS = ("geom", "nsites", "valence", "inferred", "formal_charge", "isotope", "name", "altloc",
                 "resname", "resseq", "icode", "chain", "hetero", "ss")
"""Fields of :class:`Atom` stored in a column of the same name (``symbol`` is stored as ``number``)."""

RESIDUE_KEY = ("chain", "resseq", "icode", "resname")
"""The columns that together tell one residue from another (see :meth:`Molecule.residue_key`)."""

_STRINGS = {"name", "altloc", "resname", "icode", "chain", "ss"}

_uid_lock = threading.Lock()
_next_uid = 1


def _new_token() -> int:
    """Return a number no molecule of this process has used for its topology."""
    return int(_new_uids(1)[0])


def _new_uids(count: int) -> np.ndarray:
    """Return identifiers that no other atom of this process has had."""
    global _next_uid
    with _uid_lock:
        start = _next_uid
        _next_uid += count
    return np.arange(start, start + count, dtype=np.int64)


def symbol_number(symbol: str) -> int:
    """Atomic number of an element symbol, 0 for anything that is not an element."""
    return _NUMBERS.get(symbol, 0)


def _to_column(name: str, value: Any) -> Any:
    """Turn a field value of :class:`Atom` into the value its column stores."""
    if name == "nsites":
        return -1 if value is None else int(value)
    if name == "valence":
        return math.nan if value is None else float(value)
    if name == "formal_charge":
        return UNSET_CHARGE if value is None else int(value)
    if name in _STRINGS:
        return sys.intern(str(value))
    return value


_SYMBOL_LIST: list[str] = _SYMBOLS.tolist()
_FIELD_VALUES: dict[str, Any] = {
    "number": _SYMBOL_LIST.__getitem__,
    "nsites": lambda v: None if v < 0 else int(v),
    "valence": lambda v: None if v != v else float(v),
    "formal_charge": lambda v: None if v == UNSET_CHARGE else int(v),
    "inferred": bool, "hetero": bool, "resseq": int, "isotope": int,
}
"""How a stored value of a column becomes the field value of :class:`Atom` (inverse of :func:`_to_column`)."""


def _from_column(name: str, value: Any) -> Any:
    """Turn a stored value of any column, properties included, into a plain Python value."""
    convert = _FIELD_VALUES.get(name)
    if convert is not None:
        return convert(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


class _Table:
    """Named columns of equal length whose storage is shared with copies until written.

    Each column lives in a buffer that may be longer than the table, so
    appending rows is usually a write into spare room rather than a copy. A
    table may append to or write into a buffer only when it owns it and no
    copy can see its rows; otherwise it copies the buffer first.
    """

    __slots__ = ("_owned", "_shared", "_store", "length")

    def __init__(self, columns: dict[str, np.ndarray], length: int) -> None:
        self.length = length
        self._store = dict(columns)
        self._owned = set(columns)      # buffers this table may write into
        self._shared: set[str] = set()  # buffers whose first ``length`` rows another table sees

    def __contains__(self, name: str) -> bool:
        return name in self._store

    def names(self) -> list[str]:
        return list(self._store)

    def view(self, name: str) -> np.ndarray:
        """Read-only view of a column that stays as it is: the next write to the column copies it first."""
        view = self.peek(name)
        self._shared.add(name)
        return view

    def peek(self, name: str) -> np.ndarray:
        """Read-only view of a column for immediate use inside the molecule (later writes may change it)."""
        view = self._store[name][:self.length]
        view.flags.writeable = False
        return view

    def raw(self, name: str) -> np.ndarray:
        """Return the buffer itself, for reading single values quickly (callers must not write)."""
        return self._store[name]

    def copy(self) -> "_Table":
        """Return a table that shares every buffer with this one."""
        other = _Table(self._store, self.length)
        other._owned = set()
        other._shared = set(self._store)
        self._shared.update(self._store)
        return other

    def writable(self, name: str) -> np.ndarray:
        """Buffer of a column that only this table sees, copied first when needed."""
        if name not in self._owned or name in self._shared:
            self._store[name] = self._store[name].copy()
            self._owned.add(name)
            self._shared.discard(name)
        return self._store[name]

    def replace(self, name: str, values: np.ndarray) -> None:
        """Give a column new contents of the table's length (the array is taken, not copied)."""
        self._store[name] = values
        self._owned.add(name)
        self._shared.discard(name)

    def drop(self, name: str) -> None:
        del self._store[name]
        self._owned.discard(name)
        self._shared.discard(name)

    def append(self, rows: dict[str, np.ndarray], count: int, defaults: Mapping[str, ColumnSpec]) -> None:
        """Add ``count`` rows; columns missing from ``rows`` get their default."""
        if count == 0:
            return
        new_length = self.length + count
        for name, buffer in self._store.items():
            if name not in self._owned or len(buffer) < new_length:
                grown = np.empty((max(new_length, 2 * len(buffer), 16), *buffer.shape[1:]), buffer.dtype)
                grown[:self.length] = buffer[:self.length]
                self._store[name] = buffer = grown
                self._owned.add(name)
                self._shared.discard(name)
            values = rows.get(name)
            if values is None:
                buffer[self.length:new_length] = defaults[name].default
            else:
                buffer[self.length:new_length] = values
        self.length = new_length

    def take(self, indices: np.ndarray) -> "_Table":
        """Return a new table made of some rows (copied)."""
        return _Table({name: buffer[:self.length][indices] for name, buffer in self._store.items()},
                      len(indices))


def _empty_atoms() -> _Table:
    return _Table({name: np.empty((0,) + ((spec.width,) if spec.width > 1 else ()), spec.dtype)
                   for name, spec in ATOM_COLUMNS.items()}, 0)


def _empty_bonds() -> _Table:
    return _Table({"pairs": np.empty((0, 2), np.intp), "orders": np.empty(0, np.float64)}, 0)


def _column_array(spec: ColumnSpec, values: Any, count: int) -> np.ndarray:
    """Values for ``count`` rows of a column, as an array of the column's type."""
    shape = (count, spec.width) if spec.width > 1 else (count,)
    if values is None:
        array = np.empty(shape, spec.dtype)
        array[...] = spec.default
        return array
    if spec.dtype == object:
        array = np.empty(count, object)
        array[:] = list(values) if not isinstance(values, np.ndarray) else values
        return array
    return np.asarray(values, dtype=spec.dtype).reshape(shape).copy()


class _Field[T]:
    """A field of :class:`AtomView`: reads and writes one column at the view's row."""

    def __init__(self, column: str, doc: str) -> None:
        self.column = column
        self.convert = _FIELD_VALUES.get(column)
        self.__doc__ = doc

    @overload
    def __get__(self, view: None, owner: type) -> "_Field[T]": ...

    @overload
    def __get__(self, view: "AtomView", owner: type) -> T: ...

    def __get__(self, view: "AtomView | None", owner: type) -> "T | _Field[T]":
        if view is None:
            return self
        table = view._mol._atoms
        if view.index >= table.length:
            raise IndexError(f"atom {view.index} out of range for {table.length} atoms")
        value = table._store[self.column][view.index]
        return value if self.convert is None else self.convert(value)

    def __set__(self, view: "AtomView", value: T) -> None:
        view._mol._set_value(self.column, view.index, value)


class AtomView:
    """One atom of a molecule: its fields, read from and written to the molecule's columns.

    Views are made on demand by :attr:`Molecule.atoms` and refer to a row
    number, so a view taken before atoms are deleted may point at another
    atom afterwards. The fields are those of :class:`Atom`.
    """

    __slots__ = ("_mol", "index")

    symbol = _Field[str]("number", "Element symbol.")
    geom = _Field[str | None]("geom", "Ideal arrangement of the bonds.")
    nsites = _Field[int | None]("nsites", "Number of bonds the atom should carry.")
    valence = _Field[float | None]("valence", "Bond-order sum the atom should reach.")
    inferred = _Field[bool]("inferred", "Whether the bonding pattern was guessed from a file.")
    formal_charge = _Field[int | None]("formal_charge", "Formal charge given by a file, or None.")
    isotope = _Field[int]("isotope", "Mass number of a specific isotope, 0 for the natural mixture.")
    name = _Field[str]("name", "Atom name from the file.")
    altloc = _Field[str]("altloc", "Alternate location indicator.")
    resname = _Field[str]("resname", "Residue name.")
    resseq = _Field[int]("resseq", "Residue number.")
    icode = _Field[str]("icode", "Insertion code of the residue.")
    chain = _Field[str]("chain", "Chain identifier.")
    hetero = _Field[bool]("hetero", "HETATM flag.")
    ss = _Field[str]("ss", "Secondary structure of the residue.")

    def __init__(self, mol: "Molecule", index: int) -> None:
        self._mol = mol
        self.index = index

    def __repr__(self) -> str:
        return f"<atom {self.index}: {self.symbol}>"

    def record(self) -> Atom:
        """Return the fields as a free-standing :class:`Atom`."""
        return Atom(self.symbol, **{name: getattr(self, name) for name in RECORD_FIELDS})

    def get(self, name: str) -> Any:
        """Value of any column for this atom, including properties added by plugins."""
        return self._mol._value(name, self.index)

    def set(self, name: str, value: Any) -> None:
        """Change a column for this atom."""
        self._mol._set_value(name, self.index, value)


AtomLike = Atom | AtomView
"""An atom record or an atom of a molecule: anything with the fields of :class:`Atom`."""


class AtomList(Sequence[AtomView]):
    """The atoms of a molecule as a sequence of :class:`AtomView`."""

    __slots__ = ("_mol",)

    def __init__(self, mol: "Molecule") -> None:
        self._mol = mol

    def __len__(self) -> int:
        return len(self._mol)

    @overload
    def __getitem__(self, i: SupportsIndex) -> AtomView: ...

    @overload
    def __getitem__(self, i: slice) -> list[AtomView]: ...

    def __getitem__(self, i: SupportsIndex | slice) -> AtomView | list[AtomView]:
        n = len(self._mol)
        if isinstance(i, slice):
            return [AtomView(self._mol, k) for k in range(*i.indices(n))]
        k = i.__index__()
        if k < 0:
            k += n
        if not 0 <= k < n:
            raise IndexError(f"atom {i} out of range for {n} atoms")
        return AtomView(self._mol, k)

    def __iter__(self) -> Iterator[AtomView]:
        mol = self._mol
        return (AtomView(mol, k) for k in range(len(mol)))

    def __setitem__(self, i: int, atom: Atom | AtomView) -> None:
        """Give an atom the fields of another record; it then counts as a new atom."""
        self._mol.set_atom(int(i), atom)


class _Graph:
    """Neighbour lists of a bond table as arrays (compressed sparse rows), shared by copies.

    The neighbours of atom ``i`` are ``nbr[start[i]:start[i + 1]]`` in
    ascending order, and ``row`` holds the bond row of each. It covers the
    first ``bonds`` rows of the bond table and the first ``atoms`` atoms;
    bonds added later are kept in a small overlay by the molecule until the
    arrays are rebuilt.
    """

    __slots__ = ("atoms", "bonds", "nbr", "row", "start")

    def __init__(self, pairs: np.ndarray, atoms: int) -> None:
        bonds = len(pairs)
        ends = np.concatenate([pairs[:, 0], pairs[:, 1]])
        others = np.concatenate([pairs[:, 1], pairs[:, 0]])
        order = np.lexsort((others, ends))
        self.nbr = others[order]
        self.row = np.concatenate([np.arange(bonds), np.arange(bonds)])[order]
        self.start = np.zeros(atoms + 1, dtype=np.intp)
        np.cumsum(np.bincount(ends, minlength=atoms), out=self.start[1:])
        self.atoms = atoms
        self.bonds = bonds

    def of(self, i: int) -> tuple[np.ndarray, np.ndarray]:
        """Neighbours of atom ``i`` and the bond row of each (empty for atoms added later)."""
        if i >= self.atoms:
            return self.nbr[:0], self.row[:0]
        lo, hi = self.start[i], self.start[i + 1]
        return self.nbr[lo:hi], self.row[lo:hi]


class BondMap(Mapping[Bond, float]):
    """The bonds of a molecule as a read-only mapping from ``(i, j)``, ``i < j``, to bond order.

    Change bonds with :meth:`Molecule.add_bond`, :meth:`Molecule.remove_bond`
    and :meth:`Molecule.set_bonds`; for array work use
    :attr:`Molecule.bond_pairs` and :attr:`Molecule.bond_orders`.
    """

    __slots__ = ("_mol",)

    def __init__(self, mol: "Molecule") -> None:
        self._mol = mol

    def __getitem__(self, key: Bond) -> float:
        row = self._mol._bond_row(key)
        if row is None:
            raise KeyError(key)
        return float(self._mol._bonds.raw("orders")[row])

    def __contains__(self, key: object) -> bool:
        return isinstance(key, tuple) and self._mol._bond_row(key) is not None

    def __iter__(self) -> Iterator[Bond]:
        return iter(self.keys())

    def __len__(self) -> int:
        return self._mol._bonds.length

    def keys(self) -> list[Bond]:  # type: ignore[override]
        """Bond keys in the order the bonds were made."""
        return [(i, j) for i, j in self._mol.bond_pairs.tolist()]

    def values(self) -> list[float]:  # type: ignore[override]
        """Bond orders, in the order of :meth:`keys`."""
        return self._mol.bond_orders.tolist()

    def items(self) -> list[tuple[Bond, float]]:  # type: ignore[override]
        """``((i, j), order)`` pairs, in the order the bonds were made."""
        return list(zip(self.keys(), self.values(), strict=True))

    def __repr__(self) -> str:
        return repr(dict(self.items()))


class Molecule:
    """Atoms, coordinates, bonds and pinned measurements.

    Attributes
    ----------
    measurements : list of tuple of int
        Pinned distances, angles and dihedrals as 2-4 atom indices.
    property_sources : dict
        For properties worked out from the structure (such as partial
        charges): what they were worked out from, so a cache can tell
        whether they still hold. Copies and undo snapshots keep it.
    """

    def __init__(self) -> None:
        self._atoms = _empty_atoms()
        self._bonds = _empty_bonds()
        self._properties: dict[str, ColumnSpec] = {}   # columns added beyond ATOM_COLUMNS
        self._topology = _new_token()
        self.measurements: list[tuple[int, ...]] = []
        self.property_sources: dict[str, Any] = {}
        self._graph: _Graph | None = None   # neighbour arrays of the bonds, built when first needed
        self._added: dict[int, list[tuple[int, int]]] | None = None   # (neighbour, row) of bonds after it

    def __len__(self) -> int:
        return self._atoms.length

    def __repr__(self) -> str:
        return f"<Molecule {self.formula() or 'empty'}: {len(self)} atoms, {self._bonds.length} bonds>"

    # -- making molecules ---------------------------------------------------------------------

    @classmethod
    def from_atoms(cls, atoms: Iterable[Atom | AtomView], coords: ArrayLike) -> Self:
        """Build a molecule from atom records and coordinates in one step.

        Parameters
        ----------
        atoms : iterable of Atom
            Atom records (or atoms of another molecule).
        coords : array_like
            Coordinates, shape ``(n_atoms, 3)``.

        Returns
        -------
        Molecule
            Molecule without bonds. Much faster than repeated :meth:`add_atom`.
        """
        atoms = list(atoms)
        columns: dict[str, Any] = {"number": [symbol_number(a.symbol) for a in atoms]}
        for name in RECORD_FIELDS:
            columns[name] = [_to_column(name, getattr(a, name)) for a in atoms]
        return cls.from_columns(coords, **columns)

    @classmethod
    def from_columns(cls, coords: ArrayLike, **columns: Any) -> Self:
        """Build a molecule from whole columns, the fastest way to make a large one.

        Parameters
        ----------
        coords : array_like
            Coordinates, shape ``(n_atoms, 3)``.
        **columns
            Values of other columns (see :data:`ATOM_COLUMNS`), one per atom;
            ``symbol`` may be given instead of ``number``. Columns left out
            take their default.

        Returns
        -------
        Molecule
            Molecule without bonds.
        """
        xyz = np.array(coords, dtype=np.float64).reshape(-1, 3)
        n = len(xyz)
        if "symbol" in columns:
            columns["number"] = [symbol_number(s) for s in columns.pop("symbol")]
        unknown = set(columns) - set(ATOM_COLUMNS) - {"coords", "uid"}
        if unknown:
            raise ValueError(f"no such atom column: {', '.join(sorted(unknown))}")
        store = {}
        for name, spec in ATOM_COLUMNS.items():
            if name == "coords":
                store[name] = xyz
            elif name == "uid":
                store[name] = _new_uids(n)
            else:
                values = columns.get(name)
                if values is not None and name in _STRINGS:
                    values = [sys.intern(str(v)) for v in values]
                store[name] = _column_array(spec, values, n)
                if len(store[name]) != n:
                    raise ValueError(f"column {name} has {len(store[name])} values for {n} atoms")
        mol = cls()
        mol._atoms = _Table(store, n)
        return mol

    def copy(self) -> "Molecule":
        """Return an independent copy; it shares storage with this molecule until either one changes."""
        mol = Molecule()
        mol._atoms = self._atoms.copy()
        mol._bonds = self._bonds.copy()
        mol._properties = dict(self._properties)
        mol._topology = self._topology
        mol.measurements = list(self.measurements)
        mol.property_sources = dict(self.property_sources)
        mol._graph = self._graph   # the same bonds: the arrays never change, so they can be shared
        return mol

    # -- columns ------------------------------------------------------------------------------

    @property
    def topology(self) -> int:
        """A number that changes whenever atoms, elements or bonds change, and only then.

        Copies share it while they have the same atoms and bonds, so it
        makes a good cache key for anything worked out from the bond graph.
        """
        return self._topology

    def _topology_changed(self) -> None:
        self._topology = _new_token()

    @property
    def atoms(self) -> AtomList:
        """Per-atom access: ``mol.atoms[i].symbol``, ``for a in mol.atoms``."""
        return AtomList(self)

    @property
    def coords(self) -> np.ndarray:
        """Cartesian coordinates in angstrom, shape ``(n_atoms, 3)`` (read-only).

        Assign a whole new array to move atoms (``mol.coords = mol.coords + shift``),
        or use :meth:`set_positions` for some of them.
        """
        return self._atoms.view("coords")

    @coords.setter
    def coords(self, value: ArrayLike) -> None:
        xyz = np.array(value, dtype=np.float64).reshape(-1, 3)
        if len(xyz) != len(self):
            raise ValueError(f"{len(xyz)} positions for {len(self)} atoms")
        self._atoms.replace("coords", xyz)

    def set_positions(self, indices: ArrayLike, positions: ArrayLike) -> None:
        """Move some atoms.

        Parameters
        ----------
        indices : array_like of int
            Atoms to move.
        positions : array_like
            Their new coordinates, shape ``(len(indices), 3)`` (or ``(3,)`` for all of them).
        """
        idx = np.asarray(indices, dtype=np.intp).reshape(-1)
        self._atoms.writable("coords")[idx] = np.asarray(positions, dtype=np.float64)

    @property
    def numbers(self) -> np.ndarray:
        """Atomic numbers, shape ``(n_atoms,)`` (read-only; 0 for dummy atoms)."""
        return self._atoms.view("number")

    @property
    def symbols(self) -> list[str]:
        """Element symbols of all atoms."""
        return _SYMBOLS[self._atoms.view("number")].tolist()

    @property
    def uids(self) -> np.ndarray:
        """Identifier of each atom, unique in this process and kept through edits and copies."""
        return self._atoms.view("uid")

    def column(self, name: str) -> np.ndarray:
        """All values of a column, read-only: ``"symbol"``, any field of :class:`Atom`, or a property.

        Missing values are -1 for ``nsites`` and NaN for ``valence``.
        """
        if name == "symbol":
            return _SYMBOLS[self._atoms.view("number")]
        if name not in self._atoms:
            raise KeyError(f"no atom column {name!r}")
        return self._atoms.view(name)

    def set_column(self, name: str, values: Any, atoms: ArrayLike | None = None) -> None:
        """Change a column for all atoms, or for some.

        Parameters
        ----------
        name : str
            Column (``"symbol"``, a field of :class:`Atom` or a property).
        values : scalar or array_like
            New values, one per atom given (or one for all of them).
        atoms : array_like of int, optional
            Atoms to change; all when omitted.
        """
        if name == "symbol":
            name = "number"
            values = (symbol_number(values) if isinstance(values, str)
                      else [symbol_number(v) for v in values])
        elif name in _STRINGS:
            values = (sys.intern(str(values)) if isinstance(values, str)
                      else np.array([sys.intern(str(v)) for v in values], dtype=object))
        elif name in ("nsites", "valence", "formal_charge") and not (
                isinstance(values, np.ndarray) and values.dtype != object):
            values = (_to_column(name, values) if values is None or np.isscalar(values)
                      else [_to_column(name, v) for v in values])
        if name not in self._atoms or name in ("uid", "coords"):
            raise KeyError(f"no atom column {name!r} to set")
        if name == "number":
            self._topology_changed()
        buffer = self._atoms.writable(name)
        if atoms is None:
            buffer[:len(self)] = values
        else:
            buffer[np.asarray(atoms, dtype=np.intp)] = values

    def mask(self, name: str, values: Iterable[Any]) -> np.ndarray:
        """Boolean mask of the atoms whose column ``name`` holds one of ``values``.

        Examples
        --------
        >>> water = mol.mask("resname", {"HOH", "WAT"})
        >>> carbons = mol.mask("symbol", {"C"})
        """
        wanted = set(values)
        if name == "symbol":
            numbers = np.zeros(len(_SYMBOLS), dtype=bool)
            numbers[[symbol_number(s) for s in wanted]] = True
            return numbers[self.numbers]
        column = self.column(name)
        if column.dtype != object:
            return np.isin(column, list(wanted))
        # few distinct strings: decide each once
        distinct = {value: value in wanted for value in set(column.tolist())}
        return np.fromiter((distinct[v] for v in column.tolist()), dtype=bool, count=len(column))

    def _value(self, name: str, i: int) -> Any:
        if not 0 <= i < len(self):
            raise IndexError(f"atom {i} out of range for {len(self)} atoms")
        return _from_column(name, self._atoms.raw(name)[i])

    def _set_value(self, name: str, i: int, value: Any) -> None:
        if not 0 <= i < len(self):
            raise IndexError(f"atom {i} out of range for {len(self)} atoms")
        if name == "uid":
            raise KeyError("atom identifiers are given by the molecule and cannot be set")
        if name == "number":
            value = symbol_number(value) if isinstance(value, str) else int(value)
            self._topology_changed()
        elif name not in self._properties:
            value = _to_column(name, value)
        self._atoms.writable(name)[i] = value

    # -- properties added by readers and plugins ------------------------------------------------

    @property
    def properties(self) -> list[str]:
        """Names of the columns added with :meth:`add_property`."""
        return list(self._properties)

    def add_property(self, name: str, dtype: DTypeLike = np.float64, default: Any = math.nan,
                     values: ArrayLike | None = None) -> None:
        """Add a per-atom column of your own, kept through edits, copies and undo.

        Parameters
        ----------
        name : str
            Column name; must not be one of :data:`ATOM_COLUMNS`.
        dtype : numpy dtype, optional
            Type of the values (``object`` for strings).
        default : optional
            Value of atoms added later, and of every atom unless ``values`` is given.
        values : array_like, optional
            Initial value of each atom.
        """
        if name in ATOM_COLUMNS or name == "symbol":
            raise ValueError(f"{name!r} is a built-in column")
        spec = ColumnSpec(np.dtype(dtype), default)
        self._properties[name] = spec
        self._atoms.replace(name, _column_array(spec, values, len(self)))

    def set_property(self, name: str, values: ArrayLike, atoms: ArrayLike | None = None,
                     dtype: DTypeLike = np.float64, default: Any = math.nan) -> None:
        """Set a property for all atoms or some, adding the column first if there is none.

        Parameters
        ----------
        name : str
            Property name.
        values : scalar or array_like
            New values, one per atom given (or one for all of them).
        atoms : array_like of int, optional
            Atoms to change; all when omitted.
        dtype, default : optional
            Type and default of the column if it has to be added.
        """
        if name not in self._properties:
            self.add_property(name, dtype, default)
        buffer = self._atoms.writable(name)
        if atoms is None:
            buffer[:len(self)] = values
        else:
            buffer[np.asarray(atoms, dtype=np.intp)] = values

    def property_spec(self, name: str) -> ColumnSpec:
        """Type and default of a property column."""
        return self._properties[name]

    def has_property(self, name: str) -> bool:
        """Whether a property column of this name exists."""
        return name in self._properties

    def remove_property(self, name: str) -> None:
        """Remove a column added with :meth:`add_property` (nothing happens if there is none)."""
        if self._properties.pop(name, None) is not None:
            self._atoms.drop(name)
        self.property_sources.pop(name, None)

    def _specs(self) -> dict[str, ColumnSpec]:
        return {**ATOM_COLUMNS, **self._properties}

    # -- bonds --------------------------------------------------------------------------------

    @property
    def bonds(self) -> BondMap:
        """Bond orders keyed by ``(i, j)`` with ``i < j``; 1.5 means aromatic (read-only mapping)."""
        return BondMap(self)

    @property
    def bond_pairs(self) -> np.ndarray:
        """Bonded atom pairs, ``i < j`` in each row, shape ``(n_bonds, 2)`` (read-only)."""
        return self._bonds.view("pairs")

    @property
    def bond_orders(self) -> np.ndarray:
        """Order of each bond in :attr:`bond_pairs` (read-only)."""
        return self._bonds.view("orders")

    @staticmethod
    def key(i: int, j: int) -> Bond:
        """Return the dictionary key of the bond between atoms ``i`` and ``j``."""
        return (i, j) if i < j else (j, i)

    def _graph_now(self, complete: bool = False) -> _Graph:
        """Neighbour arrays, rebuilt when missing, when many bonds came after them, or when ``complete``."""
        graph = self._graph
        later = self._bonds.length - (graph.bonds if graph is not None else 0)
        if (graph is None or (later and (complete or later > max(64, graph.bonds // 8)))
                or (complete and graph.atoms != len(self))):
            graph = self._graph = _Graph(self._bonds.peek("pairs"), len(self))
            self._added = None
        return graph

    def _later_bonds(self) -> dict[int, list[tuple[int, int]]]:
        """Neighbour and bond row of the bonds made after the neighbour arrays, by atom."""
        graph = self._graph_now()
        if self._added is None:
            added: dict[int, list[tuple[int, int]]] = {}
            pairs = self._bonds.peek("pairs")[graph.bonds:].tolist()
            for row, (i, j) in enumerate(pairs, graph.bonds):
                added.setdefault(i, []).append((j, row))
                added.setdefault(j, []).append((i, row))
            self._added = added
        return self._added

    def _bond_row(self, key: Bond) -> int | None:
        i, j = key
        if i < 0 or i >= len(self):
            return None
        nbr, rows = self._graph_now().of(i)
        if len(nbr):
            k = int(np.searchsorted(nbr, j))
            if k < len(nbr) and nbr[k] == j:
                return int(rows[k])
        for other, row in self._later_bonds().get(i, ()):
            if other == j:
                return row
        return None

    def _neighbour_rows(self, i: int) -> tuple[list[int], list[int]]:
        """Neighbours of atom ``i`` (ascending) and the bond row of each."""
        nbr, rows = self._graph_now().of(i)
        later = self._later_bonds().get(i)
        if not later:
            return nbr.tolist(), rows.tolist()
        pairs = sorted([*zip(nbr.tolist(), rows.tolist(), strict=True), *later])
        return [j for j, _ in pairs], [row for _, row in pairs]

    def _reset_bonds(self, pairs: np.ndarray, orders: np.ndarray) -> None:
        self._bonds = _Table({"pairs": pairs, "orders": orders}, len(orders))
        self._graph = None
        self._added = None
        self._topology_changed()

    # -- editing ------------------------------------------------------------------------------

    def add_atom(self, atom: Atom | AtomView | str, position: ArrayLike) -> int:
        """Append an atom.

        Parameters
        ----------
        atom : Atom, AtomView or str
            Atom record (or an atom of another molecule, whose fields are
            copied), or an element symbol.
        position : array_like
            Cartesian position, shape ``(3,)``.

        Returns
        -------
        int
            Index of the new atom.
        """
        record = Atom(atom) if isinstance(atom, str) else atom
        row: dict[str, Any] = {"number": symbol_number(record.symbol), "uid": _new_uids(1),
                               "coords": np.asarray(position, dtype=np.float64).reshape(1, 3)}
        for name in RECORD_FIELDS:
            row[name] = _to_column(name, getattr(record, name))
        self._atoms.append(row, 1, self._specs())
        self._topology_changed()
        return len(self) - 1

    def set_atom(self, i: int, atom: Atom | AtomView) -> None:
        """Give atom ``i`` the element and fields of a record; it then counts as a new atom.

        Its position and bonds stay; properties added with :meth:`add_property` keep their values.
        """
        if not 0 <= i < len(self):
            raise IndexError(f"atom {i} out of range for {len(self)} atoms")
        self._atoms.writable("number")[i] = symbol_number(atom.symbol)
        for name in RECORD_FIELDS:
            self._atoms.writable(name)[i] = _to_column(name, getattr(atom, name))
        self._atoms.writable("uid")[i] = _new_uids(1)[0]
        self._topology_changed()

    def add_bond(self, i: int, j: int, order: float = 1.0) -> None:
        """Create a bond, or change the order of an existing one."""
        if i == j:
            return
        key = self.key(int(i), int(j))
        row = self._bond_row(key)
        if row is not None:
            if self._bonds.raw("orders")[row] != float(order):
                self._bonds.writable("orders")[row] = float(order)
                self._topology_changed()
            return
        if key[0] < 0 or key[1] >= len(self):
            raise IndexError(f"bond {key} out of range for {len(self)} atoms")
        self._bonds.append({"pairs": np.array([key]), "orders": np.array([float(order)])}, 1, {})
        if self._added is not None:
            row = self._bonds.length - 1
            self._added.setdefault(key[0], []).append((key[1], row))
            self._added.setdefault(key[1], []).append((key[0], row))
        self._topology_changed()

    def remove_bond(self, i: int, j: int) -> None:
        """Delete the bond between ``i`` and ``j`` if there is one."""
        key = self.key(int(i), int(j))
        row = self._bond_row(key)
        if row is None:
            return
        keep = np.ones(self._bonds.length, dtype=bool)
        keep[row] = False
        self._reset_bonds(self._bonds.peek("pairs")[keep].copy(), self._bonds.peek("orders")[keep].copy())

    def set_bonds(self, bonds: Mapping[Bond, float] | Iterable[tuple[Bond, float]]) -> None:
        """Replace all bonds (a mapping of ``(i, j)`` to order, or such pairs)."""
        items = list(bonds.items()) if isinstance(bonds, Mapping) else list(bonds)
        merged: dict[Bond, float] = {}
        for (i, j), order in items:
            if i != j:
                merged[self.key(int(i), int(j))] = float(order)
        pairs = np.array(list(merged), dtype=np.intp).reshape(-1, 2)
        if len(pairs) and (pairs.min() < 0 or pairs.max() >= len(self)):
            raise IndexError(f"bond out of range for {len(self)} atoms")
        self._reset_bonds(pairs, np.array(list(merged.values()), dtype=np.float64))

    def set_bond_arrays(self, pairs: ArrayLike, orders: ArrayLike) -> None:
        """Replace all bonds from arrays, the fast way for readers and other bulk work.

        Parameters
        ----------
        pairs : array_like of int
            Bonded atom pairs, shape ``(n_bonds, 2)``, in either order; a pair
            given twice keeps its first place and its last order, and pairs
            of an atom with itself are dropped.
        orders : array_like of float
            Order of each bond (1.5 for aromatic).
        """
        ij = np.asarray(pairs, dtype=np.intp).reshape(-1, 2)
        values = np.asarray(orders, dtype=np.float64).reshape(-1)
        if len(values) != len(ij):
            raise ValueError(f"{len(values)} bond orders for {len(ij)} pairs")
        keep = ij[:, 0] != ij[:, 1]
        ij, values = np.sort(ij[keep], axis=1), values[keep]
        if len(ij) and (ij.min() < 0 or ij.max() >= len(self)):
            raise IndexError(f"bond out of range for {len(self)} atoms")
        # one row per pair, where it first appears, with the order it is given last (as set_bonds does)
        _, first = np.unique(ij, axis=0, return_index=True)
        _, last = np.unique(ij[::-1], axis=0, return_index=True)
        order = np.argsort(first)
        self._reset_bonds(ij[first[order]].copy(), values[len(ij) - 1 - last[order]].copy())

    def adjacency(self) -> tuple[np.ndarray, np.ndarray]:
        """Neighbours of every atom as two arrays, for vectorised work over the bond graph.

        Returns
        -------
        start : numpy.ndarray
            Shape ``(n_atoms + 1,)``; the neighbours of atom ``i`` are
            ``neighbours[start[i]:start[i + 1]]``.
        neighbours : numpy.ndarray
            Neighbour indices, ascending for each atom (read-only).
        """
        graph = self._graph_now(complete=True)
        start, nbr = graph.start[:], graph.nbr[:]
        start.flags.writeable = False
        nbr.flags.writeable = False
        return start, nbr

    def bond_order(self, i: int, j: int) -> float:
        """Order of the ``i``-``j`` bond, 0 if they are not bonded."""
        row = self._bond_row(self.key(int(i), int(j)))
        return 0.0 if row is None else float(self._bonds.raw("orders")[row])

    def remove_atoms(self, indices: Iterable[int]) -> list[int]:
        """Delete atoms with their bonds and measurements.

        Parameters
        ----------
        indices : iterable of int
            Atoms to delete.

        Returns
        -------
        list of int
            New index of every old atom, -1 for deleted ones.
        """
        n = len(self)
        drop = np.zeros(n, dtype=bool)
        idx = np.fromiter((int(i) for i in indices), dtype=np.intp)
        drop[idx[(idx >= 0) & (idx < n)]] = True
        keep = np.flatnonzero(~drop)
        mapping = np.full(n, -1, dtype=np.intp)
        mapping[keep] = np.arange(len(keep))
        self._atoms = self._atoms.take(keep)
        pairs, orders = self.bond_pairs, self.bond_orders
        kept = (mapping[pairs] >= 0).all(axis=1) if len(pairs) else np.zeros(0, dtype=bool)
        self._reset_bonds(mapping[pairs[kept]].reshape(-1, 2), orders[kept].copy())
        old = mapping.tolist()
        self.measurements = [tuple(old[i] for i in m) for m in self.measurements
                             if all(i < n and old[i] >= 0 for i in m)]
        return old

    def subset(self, indices: Iterable[int]) -> "Molecule":
        """Return a new molecule made of some atoms and the bonds among them.

        Parameters
        ----------
        indices : iterable of int
            Atoms to keep; they are renumbered in ascending order.

        Returns
        -------
        Molecule
            The extracted molecule.
        """
        idx = np.unique(np.fromiter((int(i) for i in indices), dtype=np.intp))
        position = np.full(len(self), -1, dtype=np.intp)
        position[idx] = np.arange(len(idx))
        mol = Molecule()
        mol._atoms = self._atoms.take(idx)
        mol._properties = dict(self._properties)
        pairs, orders = self.bond_pairs, self.bond_orders
        kept = (position[pairs] >= 0).all(axis=1) if len(pairs) else np.zeros(0, dtype=bool)
        mol._reset_bonds(position[pairs[kept]].reshape(-1, 2), orders[kept].copy())
        pos = position.tolist()
        mol.measurements = [tuple(pos[i] for i in m) for m in self.measurements
                            if all(i < len(pos) and pos[i] >= 0 for i in m)]
        return mol

    def extend(self, other: "Molecule") -> int:
        """Append the atoms and bonds of another molecule.

        Parameters
        ----------
        other : Molecule
            Molecule to append (it is copied; its atoms get new identifiers).

        Returns
        -------
        int
            Index of the first appended atom.
        """
        offset = len(self)
        for name, spec in other._properties.items():
            if name not in self._properties:
                self.add_property(name, spec.dtype, spec.default)
        count = len(other)
        rows = {name: other._atoms.peek(name) for name in other._atoms.names() if name in self._atoms}
        rows["uid"] = _new_uids(count)
        self._atoms.append(rows, count, self._specs())
        self._topology_changed()
        if other._bonds.length:
            pairs = other._bonds.peek("pairs") + offset
            self._bonds.append({"pairs": pairs, "orders": other._bonds.peek("orders")}, len(pairs), {})
            self._added = None   # worked out again from the new rows (or the arrays rebuilt) when needed
        return offset

    # -- graph queries ------------------------------------------------------------------------

    def neighbors(self, i: int) -> list[int]:
        """Atoms bonded to ``i``, sorted."""
        return self._neighbour_rows(i)[0]

    def degree(self, i: int) -> int:
        """Return the number of bonds of atom ``i``."""
        nbr, _ = self._graph_now().of(i)
        return len(nbr) + len(self._later_bonds().get(i, ()))

    def degrees(self) -> np.ndarray:
        """Return the number of bonds of every atom."""
        return np.bincount(self._bonds.peek("pairs").reshape(-1), minlength=len(self))

    def order_sum(self, i: int, exclude: Iterable[int] = ()) -> float:
        """Sum of the bond orders of atom ``i``, skipping bonds to ``exclude``."""
        skip = set(exclude)
        orders = self._bonds.raw("orders")
        nbrs, rows = self._neighbour_rows(i)
        return float(sum(orders[row] for j, row in zip(nbrs, rows, strict=True) if j not in skip))

    def hydrogens_on(self, i: int) -> list[int]:
        """Terminal hydrogens bonded to atom ``i``."""
        numbers = self._atoms.raw("number")
        return [j for j in self.neighbors(i) if numbers[j] == 1 and self.degree(j) == 1]

    def terminal_neighbors(self, i: int) -> list[int]:
        """Neighbours of atom ``i`` that have no other bonds."""
        return [j for j in self.neighbors(i) if self.degree(j) == 1]

    def component(self, start: int, blocked: Iterable[int] = ()) -> set[int]:
        """Atoms reachable from ``start`` without passing through ``blocked`` (breadth-first).

        Parameters
        ----------
        start : int
            First atom.
        blocked : iterable of int, optional
            Atoms the search may not enter.

        Returns
        -------
        set of int
            Reached atoms, including ``start``.
        """
        graph = self._graph_now(complete=True)
        nbr, first = graph.nbr.tolist(), graph.start.tolist()
        stop = set(blocked)
        seen = {start}
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in nbr[first[i]:first[i + 1]]:
                if j not in seen and j not in stop:
                    seen.add(j)
                    queue.append(j)
        return seen

    def components(self) -> list[set[int]]:
        """Return the connected components (separate molecules), ordered by their first atom."""
        labels = self.component_labels()
        groups: dict[int, set[int]] = {}
        for i, label in enumerate(labels.tolist()):
            groups.setdefault(label, set()).add(i)
        return list(groups.values())

    def component_labels(self) -> np.ndarray:
        """Return the connected component of each atom, numbered from 0 in the order of their first atoms."""
        parent = np.arange(len(self))
        a, b = self._bonds.peek("pairs").T
        while len(a):
            # every atom points at the root of its tree; join the trees each bond connects,
            # always under the smaller root, until no bond joins two trees
            low, high = np.minimum(parent[a], parent[b]), np.maximum(parent[a], parent[b])
            joins = low != high
            if not joins.any():
                break
            np.minimum.at(parent, high[joins], low[joins])
            while not np.array_equal(parent, parent[parent]):
                parent = parent[parent]
        # the roots are the first atoms of the components, so their order is kept
        _, labels = np.unique(parent, return_inverse=True)
        return labels.reshape(-1)

    def ring_bonds(self) -> set[Bond]:
        """Bonds that lie in a ring.

        A bond is in a ring exactly when it is not a bridge of the bond graph;
        bridges are found with Tarjan's low-link algorithm in O(atoms + bonds).

        Returns
        -------
        set of tuple of int
            Keys of the ring bonds.
        """
        graph = self._graph_now(complete=True)
        nbr, first = graph.nbr.tolist(), graph.start.tolist()
        n = len(self)
        order = [-1] * n
        low = [0] * n
        bridges: set[Bond] = set()
        counter = 0
        for root in range(n):
            if order[root] != -1:
                continue
            order[root] = low[root] = counter
            counter += 1
            stack = [(root, -1, iter(nbr[first[root]:first[root + 1]]))]
            while stack:
                v, parent, children = stack[-1]
                for w in children:
                    if w == parent:
                        continue
                    if order[w] == -1:
                        order[w] = low[w] = counter
                        counter += 1
                        stack.append((w, v, iter(nbr[first[w]:first[w + 1]])))
                        break
                    low[v] = min(low[v], order[w])
                else:
                    stack.pop()
                    if stack:
                        u = stack[-1][0]
                        low[u] = min(low[u], low[v])
                        if low[v] > order[u]:
                            bridges.add(self.key(u, v))
        return {k for k in self.bonds if k not in bridges}

    def residue_key(self, i: int) -> tuple[str, int, str, str] | None:
        """Return what tells the residue of atom ``i`` apart: chain, number, insertion code and name.

        These are the columns of :data:`RESIDUE_KEY`; every part of MILES
        that groups atoms into residues uses them. None for an atom without
        residue data.
        """
        if not self._atoms.raw("resname")[i]:
            return None
        chain, resseq, icode, resname = (self._atoms.raw(name)[i] for name in RESIDUE_KEY)
        return str(chain), int(resseq), str(icode), str(resname)

    def residue_atoms(self, i: int) -> list[int]:
        """Atoms in the same residue as atom ``i`` (its whole molecule if it has no residue)."""
        key = self.residue_key(i)
        if key is None:
            return sorted(self.component(i))
        mask = np.ones(len(self), dtype=bool)
        for name, value in zip(RESIDUE_KEY, key, strict=True):
            mask &= self._atoms.peek(name) == value
        return np.flatnonzero(mask).tolist()

    def chain_atoms(self, i: int) -> list[int]:
        """Residue atoms in the same chain as atom ``i`` (its molecule if it has no residue)."""
        if not self._atoms.raw("resname")[i]:
            return sorted(self.component(i))
        chain = self._atoms.raw("chain")[i]
        mask = (self.column("chain") == chain) & (self.column("resname") != "")
        return np.flatnonzero(mask).tolist()

    # -- properties of the whole structure ------------------------------------------------------

    def formula(self) -> str:
        """Molecular formula in Hill order, e.g. ``"C2H6O"``."""
        counts = np.bincount(self.numbers, minlength=len(_SYMBOLS))
        present = {str(_SYMBOLS[z]): int(c) for z, c in enumerate(counts) if c}
        first = ["C", "H"] if "C" in present else []
        order = [s for s in first if s in present] + sorted(s for s in present if s not in first)
        return "".join(s if present[s] == 1 else f"{s}{present[s]}" for s in order)

    def electron_count(self, charge: int = 0) -> int:
        """Return the number of electrons for a given total charge."""
        return int(self.numbers.sum(dtype=np.int64)) - charge

    def masses(self) -> np.ndarray:
        """Atomic masses in g/mol, shape ``(n_atoms,)``."""
        return _MASSES[self.numbers]

    def center_of_mass(self) -> np.ndarray:
        """Mass-weighted centre, shape ``(3,)``."""
        if not len(self):
            return np.zeros(3)
        m = self.masses()
        return (self.coords * m[:, None]).sum(axis=0) / m.sum()

    def centroid(self, indices: Iterable[int] | None = None) -> np.ndarray:
        """Geometric centre of all atoms, or of the given ones."""
        if not len(self):
            return np.zeros(3)
        return self.coords.mean(axis=0) if indices is None else self.coords[list(indices)].mean(axis=0)
