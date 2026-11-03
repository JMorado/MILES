"""How each atom of a document is shown: its own representation and whether focus mode hides it.

This is view state, kept apart from the chemistry in :class:`~miles.core.molecule.Molecule`.
An :class:`AtomDisplay` is immutable: every change returns a new one, so
undo steps can keep old ones for free. It follows the atoms by their
identifiers (:attr:`Molecule.uids <miles.core.molecule.Molecule.uids>`),
so it stays with the right atoms when atoms are added, deleted or
reordered by an edit.
"""

from collections.abc import Iterable, Iterator
from typing import Any, Self

import numpy as np
from numpy.typing import ArrayLike

from miles.core.molecule import RESIDUE_KEY, Molecule, Representation


def _frozen(array: np.ndarray) -> np.ndarray:
    array.flags.writeable = False
    return array


class AtomDisplay:
    """Representation override and focus-mode hiding of each atom.

    Parameters
    ----------
    uids : numpy.ndarray
        Identifiers of the atoms, in the molecule's order.
    rep : array_like of str, optional
        Representation of each atom (a :class:`Representation` value);
        ``""`` means the view's style. Default: all ``""``.
    hidden : array_like of bool, optional
        Atoms hidden by focus mode. Default: none.
    """

    __slots__ = ("hidden", "rep", "uids")

    def __init__(self, uids: ArrayLike, rep: Iterable[str] | None = None,
                 hidden: ArrayLike | None = None) -> None:
        self.uids = _frozen(np.array(uids, dtype=np.int64))
        n = len(self.uids)
        reps = np.full(n, "", dtype=object)
        if rep is not None:
            reps[:] = [str(r) for r in rep]
        self.rep = _frozen(reps)
        self.hidden = _frozen(np.zeros(n, dtype=bool) if hidden is None else np.array(hidden, dtype=bool))
        if len(self.rep) != n or len(self.hidden) != n:
            raise ValueError("one representation and one hidden flag are needed per atom")

    @classmethod
    def plain(cls, mol: Molecule) -> Self:
        """Every atom in the view's style, none hidden."""
        return cls(mol.uids)

    def __len__(self) -> int:
        return len(self.uids)

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, AtomDisplay) and np.array_equal(self.uids, other.uids)
                and np.array_equal(self.rep, other.rep) and np.array_equal(self.hidden, other.hidden))

    __hash__ = None  # type: ignore[assignment]

    # -- changes ------------------------------------------------------------------------------

    def with_rep(self, rep: str | Iterable[str], atoms: Iterable[int] | None = None) -> "AtomDisplay":
        """Return a copy with other representations.

        Parameters
        ----------
        rep : str or array_like of str
            One representation for all given atoms, or one per atom.
        atoms : iterable of int, optional
            Atoms to change; all when omitted.
        """
        reps = self.rep.copy()
        values = str(rep) if isinstance(rep, str) else [str(r) for r in rep]
        if atoms is None:
            reps[:] = values
        else:
            reps[np.fromiter((int(i) for i in atoms), dtype=np.intp)] = values
        return AtomDisplay(self.uids, reps, self.hidden)

    def with_hidden(self, hidden: ArrayLike) -> "AtomDisplay":
        """Return a copy where the atoms in a boolean mask (one value per atom) are hidden."""
        return AtomDisplay(self.uids, self.rep, hidden)

    def aligned(self, mol: Molecule) -> "AtomDisplay":
        """Return the display of the same atoms after an edit of the molecule.

        Atoms keep their settings wherever they moved to. New atoms are
        shown in the view's style, except that a new hydrogen takes the
        representation of the atom it is bonded to, and a new atom in a
        residue that was already there takes that residue's representation
        (so a mutated side chain looks like the rest of its residue).

        Parameters
        ----------
        mol : Molecule
            The molecule after the edit.

        Returns
        -------
        AtomDisplay
            This object if the atoms are the same, otherwise a new one.
        """
        uids = mol.uids
        if len(uids) == len(self.uids) and np.array_equal(uids, self.uids):
            return self
        n = len(uids)
        rep = np.full(n, "", dtype=object)
        hidden = np.zeros(n, dtype=bool)
        if len(self.uids) and n:
            order = np.argsort(self.uids, kind="stable")
            ranked = self.uids[order]
            position = np.clip(np.searchsorted(ranked, uids), 0, len(ranked) - 1)
            found = ranked[position] == uids
            old = order[position[found]]
            rep[found] = self.rep[old]
            hidden[found] = self.hidden[old]
            new = ~found
            if new.any() and (rep[found] != "").any():
                self._inherit(mol, rep, found, new)
        return AtomDisplay(uids, rep, hidden)

    @staticmethod
    def _inherit(mol: Molecule, rep: np.ndarray, found: np.ndarray, new: np.ndarray) -> None:
        """Give new hydrogens and new residue atoms the representation of their neighbours (in place)."""
        pairs = mol.bond_pairs
        hydrogen = new & (mol.numbers == 1)
        if hydrogen.any() and len(pairs):
            for h_end, other in ((0, 1), (1, 0)):
                rows = hydrogen[pairs[:, h_end]] & found[pairs[:, other]]
                rep[pairs[rows, h_end]] = rep[pairs[rows, other]]
        residue_atoms = np.flatnonzero(new & ~hydrogen & (mol.column("resname") != ""))
        if not len(residue_atoms):
            return
        columns = [mol.column(name) for name in RESIDUE_KEY]

        def keys(atoms: np.ndarray) -> list[tuple[Any, ...]]:
            return list(zip(*(column[atoms].tolist() for column in columns), strict=True))

        wanted = set(keys(residue_atoms))
        styled = np.flatnonzero(found & (rep != ""))
        known: dict[tuple[Any, ...], str] = {}
        for i, key in zip(styled.tolist(), keys(styled), strict=True):
            if key in wanted:
                known.setdefault(key, rep[i])
        for i, key in zip(residue_atoms.tolist(), keys(residue_atoms), strict=True):
            style = known.get(key)
            if style is not None:
                rep[i] = style

    # -- saving ---------------------------------------------------------------------------------

    def to_lists(self) -> dict[str, list[Any]]:
        """Representation and hidden flag of each atom, as JSON-ready lists."""
        return {"rep": self.rep.tolist(), "hidden": self.hidden.tolist()}

    @classmethod
    def from_lists(cls, mol: Molecule, data: dict[str, list[Any]] | None) -> Self:
        """Return the display of a freshly read molecule from :meth:`to_lists` data.

        The display is plain (view style, nothing hidden) when the data is
        missing or has the wrong number of atoms.
        """
        if not data or len(data.get("rep", ())) != len(mol) or len(data.get("hidden", ())) != len(mol):
            return cls(mol.uids)
        return cls(mol.uids, data["rep"], data["hidden"])

    def counts(self) -> dict[str, int]:
        """How many atoms have each representation override (``""`` for none)."""
        values, numbers = np.unique(self.rep.astype(str), return_counts=True)
        return dict(zip(values.tolist(), numbers.tolist(), strict=True))


def representations(mol: Molecule, display: AtomDisplay | None, default: str) -> np.ndarray:
    """Representation of every atom: its own, or ``default``, with focus-mode hiding as ``HIDDEN``."""
    if display is None or len(display) != len(mol):
        return np.full(len(mol), str(default), dtype=object)
    reps = np.where(display.rep == Representation.DEFAULT, str(default), display.rep).astype(object)
    reps[display.hidden] = str(Representation.HIDDEN)
    return reps


class DisplayEditor:
    """A molecule whose atoms also have writable ``rep`` and ``hidden``, for editing a display.

    :meth:`miles.api.DocumentView.edit_style` hands one to its function, so
    code can write ``mol.atoms[i].rep = "tube"`` and ``mol.atoms[i].hidden =
    True`` as it did when these lived in the molecule. Everything else
    (fields, coordinates, measurements) is the molecule's own.

    Parameters
    ----------
    mol : Molecule
        The structure.
    display : AtomDisplay
        Its display before the edit.
    """

    def __init__(self, mol: Molecule, display: AtomDisplay) -> None:
        object.__setattr__(self, "_mol", mol)
        object.__setattr__(self, "_uids", display.uids)
        object.__setattr__(self, "rep", display.rep.copy())
        object.__setattr__(self, "hidden", display.hidden.copy())

    @property
    def atoms(self) -> "_EditorAtoms":
        """The atoms, with ``rep`` and ``hidden`` besides the molecule's fields."""
        return _EditorAtoms(self)

    def display(self) -> AtomDisplay:
        """Return the display as edited."""
        return AtomDisplay(self._uids, self.rep, self.hidden)

    def __len__(self) -> int:
        return len(self._mol)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._mol, name)

    def __setattr__(self, name: str, value: Any) -> None:
        setattr(self._mol, name, value)


class _EditorAtoms:
    __slots__ = ("_editor",)

    def __init__(self, editor: DisplayEditor) -> None:
        self._editor = editor

    def __len__(self) -> int:
        return len(self._editor)

    def __getitem__(self, i: int) -> "_EditorAtom":
        n = len(self._editor)
        k = int(i) + n if int(i) < 0 else int(i)
        if not 0 <= k < n:
            raise IndexError(f"atom {i} out of range for {n} atoms")
        return _EditorAtom(self._editor, k)

    def __iter__(self) -> Iterator["_EditorAtom"]:
        return (_EditorAtom(self._editor, k) for k in range(len(self._editor)))


class _EditorAtom:
    __slots__ = ("_editor", "index")

    def __init__(self, editor: DisplayEditor, index: int) -> None:
        object.__setattr__(self, "_editor", editor)
        object.__setattr__(self, "index", index)

    @property
    def rep(self) -> str:
        return str(self._editor.rep[self.index])

    @rep.setter
    def rep(self, value: str) -> None:
        self._editor.rep[self.index] = str(value)

    @property
    def hidden(self) -> bool:
        return bool(self._editor.hidden[self.index])

    @hidden.setter
    def hidden(self, value: bool) -> None:
        self._editor.hidden[self.index] = bool(value)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._editor._mol.atoms[self.index], name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name in ("rep", "hidden"):
            object.__setattr__(self, name, value)
        else:
            setattr(self._editor._mol.atoms[self.index], name, value)
