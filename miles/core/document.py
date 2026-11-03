"""Documents: a molecule with its display, selection, trajectory, results, file and undo history.

This module has no interface code, so documents work in plain scripts as
well as in the window; :class:`miles.ui.document.Document` adds Qt signals.

A document announces each change with a :class:`ChangeEvent` that says
what kind of thing changed (atoms and bonds, positions, display, per-atom
properties, frames, named selections, results) and, when known, which
atoms. Listeners use the kind to redo only the work that depends on it:
moving atoms keeps the bond graph, restyling them keeps the chemistry.

Undo keeps whole snapshots of the document: molecule, display, selection,
trajectory, frame shown, results, charge model and named selections. They
cost little: a molecule copy shares its columns until one side changes
(see :mod:`miles.core.molecule`), and the display, trajectory and results
are immutable, so an undo step holds memory only for what its edit
changed. Each step remembers the kind of its edit, so undoing a move does
not make the window rework the bond graph.

What is shown is not always what the document holds. Playback blends two
frames, and an animation (a vibration, say) shows frames of its own; both
are view state. Before any edit, undo step or save the document puts its
own coordinates back, so blends and animations never reach the frames,
the undo history or a file.

The selection is kept by atom identifier (:attr:`Molecule.uids
<miles.core.molecule.Molecule.uids>`), so it stays on the same atoms when
others are added or deleted. Undo and redo keep the selection on its
atoms; when none of them is left, they bring back the selection of the
step they return to, so undoing the deletion of selected atoms selects
them again.

Partial charges are a property column of the molecule, ``"charge"``,
worked out on demand for the document's :attr:`Document.charge_model` and
kept until the atoms or bonds change. Results read from a file (grids,
normal modes, orbitals, charges, dipole) belong to the structure they were
computed for (:mod:`miles.core.results`); after an edit that changes the
atoms or bonds they are stale until undo returns to that structure.
"""

import contextlib
import os
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, replace
from enum import StrEnum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from miles.core import selection
from miles.core.charges import ChargeModel, partial_charges
from miles.core.chem import ChemError
from miles.core.display import AtomDisplay
from miles.core.frames import Trajectory
from miles.core.jobinfo import JobInfo
from miles.core.molecule import Molecule, Representation
from miles.core.results import Results
from miles.core.selection import Context
from miles.core.volume import Volume

UNDO_LIMIT = 200
"""Undo steps kept; with shared columns this bounds the count, not the memory, of old states."""

CHARGE = "charge"
"""Name of the property column that holds the partial charges of the chosen model."""


class Change(StrEnum):
    """What an edit changed."""

    STRUCTURE = "structure"     # atoms or bonds
    BONDS = "bonds"             # bond orders only
    COORDS = "coords"           # positions only
    STYLE = "style"             # per-atom display (representation, hidden) or measurements
    PROPERTIES = "properties"   # per-atom properties such as partial charges, or the charge model
    FRAME = "frame"             # another trajectory frame (or a point between two) is shown
    FRAMES = "frames"           # the frames were replaced or dropped, or an animation started or stopped
    NAMES = "names"             # named selections were defined or removed
    RESULTS = "results"         # grids were added or removed, or the results became stale or valid again


TOPOLOGY = frozenset({Change.STRUCTURE, Change.BONDS})
"""Kinds after which anything worked out from the atoms and bonds must be worked out again."""

_KEEPS_LAYOUT = frozenset({Change.COORDS, Change.FRAME, Change.NAMES, Change.RESULTS})
"""Kinds after which :attr:`Document.layout_revision` stays: nothing drawn depends on them but positions."""

_MOVES_ATOMS = frozenset({Change.STRUCTURE, Change.BONDS, Change.COORDS})
"""Kinds after which the frame shown takes the molecule's coordinates."""

_RANK = {Change.STRUCTURE: 6, Change.BONDS: 5, Change.PROPERTIES: 4, Change.STYLE: 3, Change.COORDS: 2,
         Change.FRAMES: 1, Change.NAMES: 1, Change.RESULTS: 1, Change.FRAME: 0}
"""When one undo step changed several kinds of thing, undo announces the one with the highest rank."""


@dataclass(frozen=True)
class ChangeEvent:
    """One change of a document.

    Attributes
    ----------
    kind : Change
        What changed.
    atoms : numpy.ndarray or None
        Indices of the atoms concerned, when the editor said; None means
        any of them may have changed.
    properties : tuple of str
        Names of the per-atom properties that changed (for ``PROPERTIES``).
    """

    kind: Change
    atoms: np.ndarray | None = None
    properties: tuple[str, ...] = ()


@dataclass
class _Step:
    """State an undo step goes back to, and what the edit after it changed."""

    mol: Molecule
    display: AtomDisplay
    selected: tuple[int, ...]
    trajectory: Trajectory | None
    frame: int
    results: Results
    charge_model: ChargeModel
    names: dict[str, str]
    kind: Change | None = None


class Document:
    """A molecule with its display, selection, file name, title, trajectory, results and undo stack.

    Edits go through :meth:`edit` or :meth:`changing` (or :meth:`checkpoint`
    and :meth:`notify` around changes made in place, as a drag does), which
    record an undo step, bump :attr:`revision` and announce a
    :class:`ChangeEvent` to :attr:`listeners`. An edit that fails is undone
    before the error reaches the caller.

    Attributes
    ----------
    mol : Molecule
        The structure. Undo, redo and a failed edit replace this object, so
        read it from the document each time rather than keeping it.
    info : JobInfo
        Title, charge, multiplicity, saved view and reader notes. Its
        ``volumes``, ``vibrations``, ``orbitals``, ``charges`` and ``dipole``
        show :attr:`results` while they belong to the structure, and are
        empty while they are stale.
    results : Results
        Grids, normal modes, orbitals, charges and dipole read with the
        structure, stamped with the topology they belong to.
    trajectory : Trajectory or None
        The document's frames (an animation is not one of them).
    revision : int
        Changes with every change of the document.
    layout_revision : int
        Changes with every change except atoms moving (``COORDS`` edits and
        trajectory frames), named selections and results: a cache key for
        what depends on the atoms, bonds, residues and display but not on
        positions.
    last_change : ChangeEvent or None
        The latest change, with the atoms and properties it concerns.

    Parameters
    ----------
    mol : Molecule, optional
        Initial structure.
    info : JobInfo, optional
        What a reader returned with the structure: title, charge, frames,
        results and, for sessions, the saved display, view and selection.
        The document takes copies; ``info`` itself is left as it is.
    path : str, optional
        File the structure came from.
    """

    def __init__(self, mol: Molecule | None = None, info: JobInfo | None = None,
                 path: str | None = None) -> None:
        self.mol = mol if mol is not None else Molecule()
        given = info if info is not None else JobInfo()
        self.info = replace(given, frames=None, frame_labels=None, energies=None, display=None,
                            view=dict(given.view), notes=list(given.notes))
        self.path = path
        self.modified = False
        self.revision = 0
        self.layout_revision = 0
        self.last_change: ChangeEvent | None = None
        self.listeners: list[Callable[[ChangeEvent], Any]] = []
        self.display = AtomDisplay.from_lists(self.mol, given.display)
        self._undo: list[_Step] = []
        self._redo: list[_Step] = []
        self.results = Results.from_info(given, self.mol)
        self._results_valid = True
        self._sync_info_results()
        self.trajectory: Trajectory | None = None
        if given.frames is not None:
            found = Trajectory.from_array(given.frames, given.frame_labels, given.energies)
            if found is not None and found.atom_count == len(self.mol):
                self.trajectory = found
        self._frame = self._shown_frame()
        self._animation: Trajectory | None = None
        self._animation_frame = 0
        self._rest: np.ndarray | None = None   # the document's coordinates while others are shown
        self.charge_model = ChargeModel.FILE if self.results.charges is not None else ChargeModel.GASTEIGER
        self._charge_failure: tuple[Any, str] = (None, "")
        self.selections: dict[str, str] = dict(given.view.get("selections") or {})
        self._selected: tuple[int, ...] = ()
        self._selected_cache: tuple[Any, list[int]] = (None, [])
        self._uid_order: tuple[Any, np.ndarray, np.ndarray] = (None, np.empty(0, np.intp),
                                                                np.empty(0, np.int64))
        saved = given.view.get("selection")
        if isinstance(saved, list) and all(isinstance(i, int) and 0 <= i < len(self.mol) for i in saved):
            self.select(saved)
        self._seen = (len(self.mol), self.mol.topology)

    @property
    def name(self) -> str:
        """File name, title, or "Untitled"."""
        if self.path:
            return os.path.basename(self.path)
        return self.info.title or "Untitled"

    def set_path(self, path: str) -> None:
        """Record that the document was saved to ``path``."""
        self.path = path
        self.modified = False
        self._renamed()

    # -- announcing changes (the window's document turns these into signals) -------------------

    def _announce(self, event: ChangeEvent) -> None:
        self.last_change = event
        for listener in list(self.listeners):
            listener(event)

    def _renamed(self) -> None:
        """Tell listeners that the name or the modified flag changed (the window's document does)."""

    # -- editing ------------------------------------------------------------------------------

    @contextlib.contextmanager
    def changing(self, kind: Change = Change.STRUCTURE) -> Iterator[Molecule]:
        """Change the document in a ``with`` block, as one undoable step.

        If the block raises, everything is put back as it was before it and
        the error goes on to the caller; nothing is left half changed.

        Examples
        --------
        >>> with doc.changing(Change.COORDS) as mol:
        ...     mol.set_positions([0], [[0.0, 0.0, 1.0]])
        """
        redo = list(self._redo)
        self.checkpoint()
        try:
            yield self.mol
        except BaseException:
            self._roll_back(redo)
            raise
        self.notify(kind)

    def edit[T](self, func: Callable[..., T], *args: Any, kind: Change = Change.STRUCTURE,
                **kwargs: Any) -> T:
        """Run ``func(mol, *args, **kwargs)`` as one undoable step and return its result.

        ``kind`` says what the edit changes. An edit that adds or removes
        atoms, or changes elements or bonds, counts as a change of the
        structure whatever ``kind`` says.
        """
        with self.changing(kind) as mol:
            return func(mol, *args, **kwargs)

    def edit_display(self, func: Callable[[AtomDisplay], AtomDisplay]) -> None:
        """Replace the display with ``func(display)`` as one undoable step."""
        with self.changing(Change.STYLE):
            self.display = func(self.display)

    def set_representation(self, rep: str, atoms: Iterable[int] | None = None) -> None:
        """Draw some atoms (or all) in a representation, as one undoable step (``""``: the view's style)."""
        chosen = None if atoms is None else list(atoms)
        self.edit_display(lambda display: display.with_rep(rep, chosen))

    def set_hidden(self, hidden: ArrayLike) -> None:
        """Hide the atoms in a boolean mask (one value per atom) and show the rest, as one undoable step."""
        self.edit_display(lambda display: display.with_hidden(hidden))

    def checkpoint(self, frames: bool = False) -> None:
        """Record an undo step before changing the document in place.

        The document first puts its own coordinates back if a blend of two
        frames or an animation is shown, so the step holds the document's
        state, not what happened to be on screen.

        Parameters
        ----------
        frames : bool, optional
            Accepted for older callers; every step keeps the trajectory.
        """
        self._settle(announce=False)
        self._undo.append(self._state(copy=True))
        if len(self._undo) > UNDO_LIMIT:
            self._undo.pop(0)
        self._redo.clear()

    def _state(self, copy: bool, kind: Change | None = None) -> _Step:
        return _Step(self.mol.copy() if copy else self.mol, self.display, self._selected, self.trajectory,
                     self._frame, self.results, self.charge_model, self.selections, kind)

    def _restore(self, step: _Step) -> None:
        self.mol, self.display, self._selected = step.mol, step.display, step.selected
        self.trajectory, self._frame, self.results = step.trajectory, step.frame, step.results
        self.charge_model, self.selections = step.charge_model, step.names

    def _roll_back(self, redo: list[_Step]) -> None:
        """Put back the state of the last checkpoint after an edit failed, and forget the step."""
        self._restore(self._undo.pop())
        self._redo = redo
        self._sync_info_results()
        self.revision += 1
        self.layout_revision += 1
        self._seen = (len(self.mol), self.mol.topology)
        self._announce(ChangeEvent(Change.STRUCTURE))   # the molecule object was replaced

    def _step(self, source: list[_Step], target: list[_Step]) -> None:
        """Go to the state at the end of ``source``, keeping the current one on ``target``."""
        animation_ended = self._end_animation()
        self._settle(announce=False)
        step = source.pop()
        target.append(self._state(copy=False, kind=step.kind))
        frames_changed = animation_ended or step.trajectory is not self.trajectory
        current = self._selected
        self._restore(step)
        kept = self._resolve(current)
        if kept:   # the selection made since stays on its atoms; only an empty one takes the step's
            self._selected = tuple(self.mol.uids[np.asarray(kept, dtype=np.intp)].tolist())
        if frames_changed:
            self._announce(ChangeEvent(Change.FRAMES))
        self._changed(step.kind or Change.STRUCTURE, None, (), history=True)

    def notify(self, kind: Change = Change.STRUCTURE, atoms: Iterable[int] | None = None,
               properties: Iterable[str] = ()) -> None:
        """Announce a change made in place.

        Parameters
        ----------
        kind : Change, optional
            What changed; a change of atoms, elements or bonds counts as
            ``STRUCTURE`` whatever is given.
        atoms : iterable of int, optional
            The atoms concerned, if known.
        properties : iterable of str, optional
            The per-atom properties that changed.
        """
        self._changed(Change(kind), atoms, properties, history=False)

    def _changed(self, kind: Change, atoms: Iterable[int] | None, properties: Iterable[str],
                 history: bool) -> None:
        shape = (len(self.mol), self.mol.topology)
        if shape != self._seen and kind not in TOPOLOGY:
            kind = Change.STRUCTURE   # the editor said less than it did
        self._seen = shape
        if kind in _MOVES_ATOMS:
            self._rest = None   # the coordinates now shown are the edited ones
            if self._animation is not None:   # it was made for the geometry before the edit
                self._animation = None
                self._announce(ChangeEvent(Change.FRAMES))
        self.revision += 1
        if kind not in _KEEPS_LAYOUT:
            self.layout_revision += 1
        self.display = self.display.aligned(self.mol)
        if not history and self._undo:
            step = self._undo[-1]
            if step.kind is None or _RANK[kind] > _RANK[step.kind]:
                step.kind = kind
        self._follow_trajectory(kind)
        was_valid = self._results_valid
        self._sync_info_results()
        if not self.modified:
            self.modified = True
            self._renamed()
        which = None if atoms is None else np.fromiter((int(i) for i in atoms), dtype=np.intp)
        self._announce(ChangeEvent(kind, which, tuple(properties)))
        if self._results_valid != was_valid and kind != Change.RESULTS:
            self._announce(ChangeEvent(Change.RESULTS))

    def can_undo(self) -> bool:
        """Whether there is a step to undo."""
        return bool(self._undo)

    def can_redo(self) -> bool:
        """Whether there is a step to redo."""
        return bool(self._redo)

    def undo(self) -> None:
        """Go back one step."""
        if self._undo:
            self._step(self._undo, self._redo)

    def redo(self) -> None:
        """Go forward one step."""
        if self._redo:
            self._step(self._redo, self._undo)

    def session_info(self, view: dict[str, Any] | None = None) -> JobInfo:
        """Everything a session keeps besides the molecule.

        That is the frames, display, view, selection and named selections,
        and the results while they belong to the structure (stale ones
        cannot belong to it again once the history is gone). The
        document's own coordinates are put back first if others are shown.
        """
        self._settle(announce=True)
        view = dict(view or {})
        if self.selections:
            view["selections"] = dict(self.selections)
        if self.selected:
            view["selection"] = list(self.selected)
        traj = self.trajectory
        results = self.results if self.results.belong_to(self.mol) else Results()
        return replace(self.info, frames=None if traj is None else traj.array,
                       frame_labels=None if traj is None or traj.labels is None else list(traj.labels),
                       energies=None if traj is None or traj.energies is None else list(traj.energies),
                       view=view, display=self.display.to_lists(), notes=list(self.info.notes),
                       **results.as_info_fields())

    # -- the selection ------------------------------------------------------------------------

    @property
    def selected(self) -> list[int]:
        """Selected atoms, as indices in the order they were chosen.

        The document keeps the atoms' identifiers, so the selection stays
        on the same atoms when others are added or deleted; atoms deleted
        since drop out. Undo and redo keep it, or bring back the step's
        own selection when none of its atoms is left. The list is shared:
        do not change it in place.
        """
        key = (self._selected, id(self.mol), self.mol.topology, len(self.mol))
        if self._selected_cache[0] != key:
            self._selected_cache = (key, self._resolve(self._selected))
        return self._selected_cache[1]

    def select(self, atoms: Iterable[int]) -> None:
        """Select atoms by index, in the order given (repeats are dropped); not an undo step.

        Raises
        ------
        IndexError
            If an index is out of range.
        """
        chosen = list(dict.fromkeys(int(i) for i in atoms))
        n = len(self.mol)
        if any(not 0 <= i < n for i in chosen):
            raise IndexError(f"Atom indices must be between 0 and {n - 1}")
        uids = self.mol.uids
        self._selected = tuple(uids[np.asarray(chosen, dtype=np.intp)].tolist()) if chosen else ()

    def _resolve(self, wanted: tuple[int, ...]) -> list[int]:
        """Return the indices of the atoms with these identifiers, in order; missing ones drop out."""
        if not wanted or not len(self.mol):
            return []
        key = (id(self.mol), self.mol.topology, len(self.mol))
        if self._uid_order[0] != key:
            uids = self.mol.uids
            order = np.argsort(uids, kind="stable")
            self._uid_order = (key, order, uids[order])
        _, order, ordered = self._uid_order
        ids = np.fromiter(wanted, dtype=np.int64, count=len(wanted))
        positions = np.minimum(np.searchsorted(ordered, ids), len(ordered) - 1)
        found = ordered[positions] == ids
        return order[positions[found]].tolist()

    # -- selecting with expressions -----------------------------------------------------------

    def selection_context(self, selected: Iterable[int] | None = None) -> Context:
        """Return what expressions may refer to: selection, atoms not drawn, named selections, charges.

        Parameters
        ----------
        selected : iterable of int, optional
            Atoms the word ``selected`` stands for; the document's selection by default.
        """
        hidden = self.display.hidden | (self.display.rep == str(Representation.HIDDEN))
        chosen = self.selected if selected is None else selected
        return Context(selected=chosen, hidden=hidden, names=self.selections, columns={CHARGE: self.charges})

    def find(self, expression: str, selected: Iterable[int] | None = None) -> np.ndarray:
        """Return the atoms an expression of the selection language picks (:mod:`miles.core.selection`).

        Parameters
        ----------
        expression : str
            For example ``"byres within 5 of ligand"``.
        selected : iterable of int, optional
            Atoms the word ``selected`` stands for; the document's selection by default.

        Returns
        -------
        numpy.ndarray
            Indices (0-based, ascending).

        Raises
        ------
        SelectionError
            If the expression is not valid; the error says where.
        """
        return selection.find(self.mol, expression, self.selection_context(selected))

    def define(self, name: str, expression: str | None) -> None:
        """Name an expression so later ones can use the name as a word; None removes the name.

        The expression is kept as text and worked out each time the name is
        used, so it follows edits. Defining or removing a name is an undoable
        step; session files keep the names.

        Raises
        ------
        ValueError
            If the name is not a plain word or is a keyword, or the expression
            is not valid (:class:`~miles.core.selection.SelectionError`).
        """
        if expression is None:
            if name not in self.selections:
                return
            names = {key: value for key, value in self.selections.items() if key != name}
        else:
            selection.check_name(name)
            names = {**self.selections, name: expression}
            selection.evaluate(self.mol, name, replace(self.selection_context(), names=names))
        self.checkpoint()
        self.selections = names
        self.notify(Change.NAMES)

    # -- results read with the structure ------------------------------------------------------

    def _sync_info_results(self) -> None:
        """Show the results in :attr:`info` while they belong to the structure, and nothing while stale."""
        if self.results.is_empty() and not self.results.belong_to(self.mol):
            self.results = Results(self.mol.topology)   # no results: they simply follow the structure
        valid = self.results.belong_to(self.mol)
        shown = self.results if valid else Results()   # fresh lists: nothing kept depends on them
        for name, value in shown.as_info_fields().items():
            setattr(self.info, name, value)
        self._results_valid = valid

    def stale_results(self) -> list[str]:
        """Return what the results hold while they belong to another structure; empty when they apply."""
        return [] if self._results_valid else self.results.names()

    def add_volume(self, volume: Volume) -> int:
        """Add a grid to the results, as one undoable step; returns its index.

        If the results are stale, the grid starts new results for the
        structure as it is now (undo brings the old ones back).
        """
        base = self.results if self.results.belong_to(self.mol) else Results(self.mol.topology)
        self.checkpoint()
        self.results = base.with_volumes([*base.volumes, volume], self.mol.topology)
        self.notify(Change.RESULTS)
        return len(self.results.volumes) - 1

    def keep_volume(self, volume: Volume) -> int:
        """Keep a grid worked out from the document's own results, such as an orbital, with the results.

        Unlike :meth:`add_volume` it is not an edit: no undo step, and the
        document is not marked modified. Sessions save it with the other
        grids, so it need not be worked out again; an undo of an earlier
        edit drops it, like everything that came after that edit.

        Returns
        -------
        int
            Index of the grid among the results.
        """
        base = self.results if self.results.belong_to(self.mol) else Results(self.mol.topology)
        self.results = base.with_volumes([*base.volumes, volume], self.mol.topology)
        self._sync_info_results()
        self._announce(ChangeEvent(Change.RESULTS))
        return len(self.results.volumes) - 1

    def remove_volume(self, index: int) -> Volume:
        """Remove a grid from the results, as one undoable step, and return it.

        Raises
        ------
        IndexError
            If there is no such grid among the results that apply.
        """
        volumes = self.info.volumes
        if not 0 <= index < len(volumes):
            raise IndexError(f"There is no grid {index}")
        removed = volumes[index]
        self.checkpoint()
        self.results = self.results.with_volumes([v for k, v in enumerate(volumes) if k != index],
                                                 self.mol.topology)
        self.notify(Change.RESULTS)
        return removed

    # -- partial charges ----------------------------------------------------------------------

    def file_charges_apply(self) -> bool:
        """Whether the file gave charges and the atoms and bonds are still those of the file."""
        return self.results.charges is not None and self.results.belong_to(self.mol)

    def set_charge_model(self, model: str) -> None:
        """Take partial charges from the file or work them out with another model, as one undoable step."""
        chosen = ChargeModel(model)
        if chosen == self.charge_model:
            return
        self.checkpoint()
        self.charge_model = chosen
        self.notify(Change.PROPERTIES, properties=(CHARGE,))

    def charges(self) -> np.ndarray | None:
        """Partial charges of the chosen model, or None if there are none (see :meth:`charge_error`).

        They are worked out once and kept in the molecule's ``"charge"``
        column until the atoms or bonds change.
        """
        return self._charge_result()[0]

    def charge_error(self) -> str:
        """Why the chosen model has no charges for this structure; empty when it has."""
        return self._charge_result()[1]

    def _charge_result(self) -> tuple[np.ndarray | None, str]:
        mol, model = self.mol, self.charge_model
        key = (mol.topology, str(model))
        if mol.has_property(CHARGE) and mol.property_sources.get(CHARGE) == key:
            return mol.column(CHARGE), ""
        if self._charge_failure[0] == key:
            return None, self._charge_failure[1]
        if model == ChargeModel.FILE:
            if not self.file_charges_apply() or self.results.charges is None:
                self._charge_failure = (key, "The charges from the file no longer match the atoms")
                return None, self._charge_failure[1]
            values = np.asarray(self.results.charges, dtype=float)
        else:
            try:
                values = partial_charges(mol, model)
            except ChemError as exc:
                self._charge_failure = (key, str(exc))
                return None, self._charge_failure[1]
        mol.set_property(CHARGE, values)
        mol.property_sources[CHARGE] = key
        return mol.column(CHARGE), ""

    # -- trajectory frames ------------------------------------------------------------------

    def _shown(self) -> Trajectory | None:
        """Return the frames playback runs through: an animation if one is shown, else the document's own."""
        return self._animation if self._animation is not None else self.trajectory

    def _shown_frame(self) -> int:
        """Index of the frame the molecule shows (quantum chemistry outputs show their last geometry)."""
        traj = self.trajectory
        if traj is None:
            return 0
        matches = [k for k in range(len(traj)) if np.allclose(traj.coords(k), self.mol.coords)]
        return matches[-1] if matches else 0

    @property
    def frames(self) -> np.ndarray | None:
        """Coordinates of every frame shown, shape ``(n_frames, n_atoms, 3)``, read-only, or None."""
        shown = self._shown()
        return None if shown is None else shown.array

    @property
    def frame_labels(self) -> list[str] | None:
        """Text of each frame shown, or None."""
        shown = self._shown()
        return None if shown is None or shown.labels is None else list(shown.labels)

    @property
    def frame_energies(self) -> list[float] | None:
        """Energy of each frame shown in kcal/mol, or None unless every frame has one."""
        shown = self._shown()
        return None if shown is None or shown.energies is None else list(shown.energies)

    @property
    def frames_loop(self) -> bool:
        """Whether the frames shown are a periodic motion, so playback blends the last into the first."""
        shown = self._shown()
        return shown is not None and shown.loop

    @property
    def frame_count(self) -> int:
        """Number of frames shown (0 without a trajectory)."""
        shown = self._shown()
        return 0 if shown is None else len(shown)

    @property
    def frame(self) -> int:
        """Index of the frame shown."""
        return self._animation_frame if self._animation is not None else self._frame

    @frame.setter
    def frame(self, k: int) -> None:
        self.set_frame(k)

    @property
    def animating(self) -> bool:
        """Whether an animation is shown instead of the document's own frames (:meth:`show_animation`)."""
        return self._animation is not None

    def set_frames(self, frames: ArrayLike | None, labels: list[str] | None = None,
                   loop: bool = False, energies: list[float] | None = None) -> None:
        """Replace the trajectory (e.g. with conformers) and show its first frame.

        This is not an undo step by itself: callers that want one call
        :meth:`checkpoint` first (``DocumentView.set_frames`` does).

        Parameters
        ----------
        frames : array_like or None
            Coordinates, shape ``(n_frames, n_atoms, 3)``; None or an empty
            array removes the trajectory.
        labels : list of str, optional
            Text for each frame.
        loop : bool, optional
            The motion is periodic, so playback may blend the last frame into the first.
        energies : list of float, optional
            Energy of each frame in kcal/mol; the timeline plots them.

        Raises
        ------
        ValueError
            If the frames do not have one position for each atom.
        """
        found = None if frames is None else Trajectory.from_array(frames, labels, energies, loop)
        if found is not None and found.atom_count != len(self.mol):
            raise ValueError(f"The frames have {found.atom_count} atoms; the structure has {len(self.mol)}")
        self._animation = None
        self._settle(announce=False)
        self.trajectory = found
        self._frame = 0
        if found is not None and len(self.mol):
            self.mol.coords = found.coords(0)
        self.revision += 1
        self.layout_revision += 1
        self._announce(ChangeEvent(Change.FRAMES))

    def show_animation(self, frames: ArrayLike, labels: list[str] | None = None, loop: bool = True) -> None:
        """Play other frames, such as a vibration, without changing the document.

        The animation takes the place of the document's frames in playback
        until :meth:`stop_animation`, an edit that moves atoms, or undo.
        It is view state: edits, undo steps and saved files never see its
        coordinates, and the document's own come back when it stops.

        Parameters
        ----------
        frames : array_like
            Coordinates, shape ``(n_frames, n_atoms, 3)``; an empty array stops the animation.
        labels : list of str, optional
            Text for each frame.
        loop : bool, optional
            The motion is periodic (the default), so playback blends the last frame into the first.

        Raises
        ------
        ValueError
            If the frames do not have one position for each atom.
        """
        found = Trajectory.from_array(frames, labels, None, loop)
        if found is None:
            self.stop_animation()
            return
        if found.atom_count != len(self.mol):
            raise ValueError(f"The frames have {found.atom_count} atoms; the structure has {len(self.mol)}")
        if self._rest is None:
            self._rest = np.array(self.mol.coords)
        self._animation, self._animation_frame = found, 0
        self.mol.coords = found.coords(0)
        self.revision += 1
        self.layout_revision += 1
        self._announce(ChangeEvent(Change.FRAMES))

    def stop_animation(self) -> None:
        """Stop the animation and show the document's own coordinates and frames again."""
        if self._end_animation():
            self.revision += 1
            self.layout_revision += 1
            self._announce(ChangeEvent(Change.FRAMES))

    def _end_animation(self) -> bool:
        """Drop the animation and put the document's coordinates back, without announcing it."""
        if self._animation is None:
            return False
        self._animation, self._animation_frame = None, 0
        self._settle(announce=False)
        return True

    def show_geometry(self, coords: ArrayLike) -> None:
        """Show other coordinates without changing the document, e.g. a preview.

        The document's own coordinates come back with the next edit, undo,
        save, frame or :meth:`settle`; use :meth:`edit` to change them.
        """
        shown = np.array(coords, dtype=float)
        if shown.shape != (len(self.mol), 3):
            raise ValueError(f"Coordinates must have the shape {(len(self.mol), 3)}, not {shown.shape}")
        if self._rest is None:
            self._rest = np.array(self.mol.coords)
        self.mol.coords = shown
        self.revision += 1
        self._announce(ChangeEvent(Change.FRAME))

    def settle(self) -> None:
        """Show the document's own coordinates again if a blend or other coordinates are shown."""
        self._settle(announce=True)

    def _settle(self, announce: bool) -> None:
        if self._rest is None:
            return
        self.mol.coords = self._rest
        self._rest = None
        if announce:
            self.revision += 1
            self._announce(ChangeEvent(Change.FRAME))

    def set_frame(self, k: int) -> None:
        """Show frame ``k`` (wrapping around)."""
        shown = self._shown()
        if shown is None or not len(self.mol):
            return
        k = int(k) % len(shown)
        if self._animation is not None:
            self._animation_frame = k
        else:
            self._frame = k
            self._rest = None
        self.mol.coords = shown.coords(k)
        self.revision += 1
        self._announce(ChangeEvent(Change.FRAME))

    def set_position(self, position: float) -> None:
        """Show a point between two frames, blending their coordinates (smooth playback).

        ``position`` counts frames from 0; its whole part is the current
        frame. Between the last frame and the first, coordinates are blended
        only for periodic motion; a trajectory jumps back to the start. A
        blend is view state: the next edit, undo or save puts the frame's
        own coordinates back first.
        """
        shown = self._shown()
        if shown is None or not len(self.mol):
            return
        k, t, coords = shown.blend(position)
        if self._animation is not None:
            self._animation_frame = k
        else:
            self._frame = k
            self._rest = shown.coords(k) if t else None
        self.mol.coords = coords
        self.revision += 1
        self._announce(ChangeEvent(Change.FRAME))

    def _follow_trajectory(self, kind: Change) -> None:
        """Keep an edit of the positions in the frame shown; drop the frames when atoms come or go."""
        traj = self.trajectory
        if traj is None:
            return
        if traj.atom_count != len(self.mol):
            self.trajectory, self._frame = None, 0   # undo brings them back with the atoms
            self._announce(ChangeEvent(Change.FRAMES))
        elif kind in _MOVES_ATOMS and not np.array_equal(traj.coords(self._frame), self.mol.coords):
            self.trajectory = traj.with_frame(self._frame, self.mol.coords)
