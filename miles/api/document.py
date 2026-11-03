"""One tab as scripts and plugins see it: :class:`DocumentView`."""

import contextlib
import math
from collections.abc import Callable, Iterable, Iterator
from typing import TYPE_CHECKING, Any, TypeVar

import numpy as np

from miles.core.molecule import Molecule
from miles.core.patterns import match_pattern
from miles.files.common import JobInfo

if TYPE_CHECKING:
    from PySide6.QtGui import QImage

    from miles.core.display import AtomDisplay
    from miles.core.document import ChangeEvent
    from miles.core.results import Results
    from miles.core.volume import Volume
    from miles.render.scene import Style
    from miles.ui.canvas import Canvas
    from miles.ui.mainwindow import MainWindow

from miles.api.picture import Picture

T = TypeVar("T")

_VIEW = "_miles_document_view"   # attribute of a canvas holding its DocumentView


def document_view(window: "MainWindow", canvas: "Canvas") -> "DocumentView":
    """Return the handle of a tab: the same object every time for the same tab.

    Handles live as long as their tab, so plugins can use them as keys of
    dictionaries and weak dictionaries.
    """
    view = getattr(canvas, _VIEW, None)
    if not isinstance(view, DocumentView):
        view = DocumentView(window, canvas)
        setattr(canvas, _VIEW, view)
    return view


class DocumentView:
    """One open structure: molecule, selection, frames and view.

    :attr:`App.current`, :attr:`App.documents` and the callbacks of plugins
    hand out the same object for the same tab, as long as the tab is open,
    so a handle can be the key of a dictionary.
    """

    def __init__(self, window: "MainWindow", canvas: "Canvas") -> None:
        self._window = window
        self._canvas = canvas
        self._caches: dict[str, dict[str, Any]] = {}

    def __eq__(self, other: object) -> bool:
        return isinstance(other, DocumentView) and other._canvas is self._canvas

    def __hash__(self) -> int:
        return id(self._canvas)

    def __repr__(self) -> str:
        return f"<{self.name}: {len(self.molecule)} atoms>"

    @property
    def canvas(self) -> "Canvas":
        """The Qt canvas of the tab, for experiments in the console.

        It is internal: it may change in any release, so plugins use the
        methods of this class instead.
        """
        return self._canvas

    @property
    def name(self) -> str:
        """Tab name."""
        return self._canvas.doc.name

    @property
    def path(self) -> str | None:
        """File the structure came from or was saved to."""
        return self._canvas.doc.path

    @property
    def molecule(self) -> Molecule:
        """The structure; change it inside :meth:`edit` so the change can be undone."""
        return self._canvas.mol

    @property
    def info(self) -> JobInfo:
        """What came with the file: title, charge, grids, orbitals, vibrations..."""
        return self._canvas.doc.info

    @property
    def revision(self) -> int:
        """A number that changes whenever the document does (atoms, coordinates, frames, styles)."""
        return self._canvas.doc.revision

    @property
    def layout_revision(self) -> int:
        """A number raised by every change except moved atoms and frames.

        Unlike :attr:`revision`, it stays the same while atoms are dragged
        or a trajectory plays, so work that does not depend on positions
        (bond graphs, residue lists) can be kept until it changes.
        """
        return self._canvas.doc.layout_revision

    @property
    def display(self) -> "Style":
        """The tab's display options (read them here, change them with :meth:`style`)."""
        return self._canvas.display

    def activate(self) -> None:
        """Make this the current tab."""
        self._window.tabbar.setCurrentIndex(self._window.canvases.index(self._canvas))

    def cache(self, name: str) -> dict[str, Any]:
        """Return a dictionary kept for this tab under a name, and dropped when the tab closes.

        It is for working data that is not saved (a plugin gets its own with
        :meth:`~miles.plugins.PluginContext.document_cache`); data to save in
        session files belongs in
        :meth:`~miles.plugins.PluginContext.document_state`.
        """
        return self._caches.setdefault(name, {})

    def forget_cache(self, name: str) -> None:
        """Drop the dictionary :meth:`cache` keeps under a name."""
        self._caches.pop(name, None)

    # -- editing ----------------------------------------------------------------------------

    def edit(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        """Change the molecule as one undoable step.

        Parameters
        ----------
        func : callable
            Called as ``func(molecule, *args, **kwargs)``; it changes the
            molecule in place, e.g. :func:`miles.add_hydrogens <miles.core.building.add_hydrogens>`.
        *args, **kwargs
            Passed on to ``func``, except ``kind``.
        kind : str, optional
            What the edit changes, a :class:`~miles.core.document.Change`
            value: ``"structure"`` (the default: atoms or bonds),
            ``"bonds"`` (bond orders only) or ``"coords"`` (positions only),
            so that what does not depend on it is not worked out again.

        Returns
        -------
        object
            What ``func`` returns.
        """
        return self._canvas.doc.edit(func, *args, **kwargs)

    @contextlib.contextmanager
    def editing(self) -> Iterator[Molecule]:
        """Context manager for an undoable change made in several statements.

        Examples
        --------
        >>> with doc.editing() as mol:
        ...     mol.coords = mol.coords + [0.0, 0.0, 1.0]
        ...     mol.atoms[0].symbol = "N"

        If the block raises, the molecule is put back as it was and the
        error goes on; no half-made change is left behind.
        """
        with self._canvas.doc.changing() as mol:
            yield mol

    def set_molecule(self, molecule: Molecule) -> None:
        """Replace the structure with another, as one undoable step (e.g. the result of background work)."""
        doc = self._canvas.doc
        doc.checkpoint()
        doc.mol = molecule
        doc.notify()

    def edit_style(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        """Change per-atom display settings (``rep``, ``hidden``) or measurements as one undoable step.

        ``func`` gets the molecule with ``rep`` and ``hidden`` writable on
        its atoms (a :class:`~miles.core.display.DisplayEditor`):
        ``doc.edit_style(lambda mol: setattr(mol.atoms[0], "rep", "tube"))``.
        The structure does not count as changed, so nothing that depends on
        the atoms or bonds is worked out again. For whole groups of atoms,
        :meth:`set_representation`, :meth:`hide` and :meth:`show` are simpler.
        """
        from miles.core.display import DisplayEditor
        from miles.ui.document import Change

        doc = self._canvas.doc
        with doc.changing(Change.STYLE):
            editor = DisplayEditor(doc.mol, doc.display)
            result = func(editor, *args, **kwargs)
            doc.display = editor.display()
        return result

    # -- partial charges ----------------------------------------------------------------------

    @property
    def charge_model(self) -> str:
        """Where partial charges come from: ``"file"``, ``"gasteiger"`` or ``"mmff94"``."""
        return str(self._canvas.doc.charge_model)

    def set_charge_model(self, model: str) -> None:
        """Take partial charges from the file or work them out with another model."""
        self._canvas.doc.set_charge_model(model)

    def charges(self) -> np.ndarray | None:
        """Partial charge of each atom from the chosen model, or None (see :meth:`charge_error`).

        They are kept in the molecule's ``"charge"`` property column until the atoms or bonds change.
        """
        return self._canvas.doc.charges()

    def charge_error(self) -> str:
        """Why the chosen model has no charges for this structure; empty when it has."""
        return self._canvas.doc.charge_error()

    def file_charges_apply(self) -> bool:
        """Whether the file gave charges and the atoms and bonds are still those of the file."""
        return self._canvas.doc.file_charges_apply()

    # -- display ------------------------------------------------------------------------------

    def edit_display(self, func: Callable[["AtomDisplay"], "AtomDisplay"]) -> None:
        """Replace the display with ``func(display)`` as one undoable step.

        Examples
        --------
        >>> doc.edit_display(lambda d: d.with_rep("tube", pocket).with_rep("ball", ligand))
        """
        self._canvas.doc.edit_display(func)

    @property
    def last_change(self) -> "ChangeEvent | None":
        """The latest change of the document: its kind and, when known, the atoms and properties concerned.

        A plugin's :meth:`~miles.plugins.PluginContext.on_document_changed`
        callback can read it to skip work, e.g. when only some atoms moved.
        """
        return self._canvas.doc.last_change

    @property
    def atom_display(self) -> "AtomDisplay":
        """Each atom's own representation and whether focus mode hides it (read-only arrays)."""
        return self._canvas.doc.display

    def set_representation(self, rep: str, atoms: Iterable[int] | None = None) -> None:
        """Draw atoms (all when omitted) in a representation, as one undoable step.

        Parameters
        ----------
        rep : str
            A :class:`~miles.core.molecule.Representation` value, such as
            ``"tube"``; ``""`` gives the atoms back the view's style.
        atoms : iterable of int, optional
            Atoms to change.
        """
        self._canvas.doc.set_representation(rep, None if atoms is None else [int(i) for i in atoms])

    def hide(self, atoms: Iterable[int]) -> None:
        """Hide atoms (focus mode) as one undoable step; the others keep their state."""
        hidden = self._canvas.doc.display.hidden.copy()
        hidden[[int(i) for i in atoms]] = True
        self._canvas.doc.set_hidden(hidden)

    def show(self, atoms: Iterable[int] | None = None) -> None:
        """Show hidden atoms again (all when omitted) as one undoable step."""
        hidden = self._canvas.doc.display.hidden.copy()
        hidden[slice(None) if atoms is None else [int(i) for i in atoms]] = False
        self._canvas.doc.set_hidden(hidden)

    def undo(self) -> None:
        """Undo the last change."""
        self._canvas.doc.undo()

    # -- selection --------------------------------------------------------------------------

    @property
    def selection(self) -> list[int]:
        """Selected atom indices (0-based)."""
        return list(self._canvas.selection)

    @selection.setter
    def selection(self, atoms: Iterable[int]) -> None:
        self.select(atoms)

    def select(self, what: str | Iterable[int] | None = None, smarts: str | None = None,
               add: bool = False, *, atoms: Iterable[int] | None = None) -> list[int]:
        """Select atoms with an expression, by index or by a SMARTS pattern.

        Parameters
        ----------
        what : str or iterable of int, optional
            An expression of the selection language, such as
            ``"byres within 5 of ligand"`` (the selection language of the user guide),
            or atom indices (0-based). Nothing clears the selection.
        smarts : str, optional
            Pattern whose matches are selected, e.g. ``"c1ccccc1"``
            (the same as the expression ``smarts "c1ccccc1"``).
        add : bool, optional
            Add to the current selection instead of replacing it.
        atoms : iterable of int, optional
            Atom indices, as a keyword.

        Returns
        -------
        list of int
            The selection afterwards.

        Raises
        ------
        ValueError
            If an index is out of range, the expression is not valid
            (:class:`~miles.core.selection.SelectionError`, which says where)
            or the pattern is invalid.
        """
        chosen = set(self.selection) if add else set()
        n = len(self.molecule)
        if isinstance(what, str):
            chosen.update(self.find(what))
        elif what is not None:
            atoms = [*(atoms or ()), *what]
        if atoms is not None:
            wanted = [int(i) for i in atoms]
            if any(not 0 <= i < n for i in wanted):
                raise ValueError(f"Atom indices must be between 0 and {n - 1}")
            chosen.update(wanted)
        if smarts is not None:
            chosen.update(match_pattern(self.molecule, smarts).atoms)
        self._canvas.set_selection(sorted(chosen))
        return self.selection

    def find(self, expression: str, selected: Iterable[int] | None = None) -> list[int]:
        """Return the atoms an expression picks, without selecting them.

        Parameters
        ----------
        expression : str
            For example ``"protein and name CA"`` or ``"water within 3.5 of ligand"``.
        selected : iterable of int, optional
            Atoms the word ``selected`` stands for; the current selection by default.

        Returns
        -------
        list of int
            Atom indices (0-based, ascending).

        Raises
        ------
        ValueError
            If the expression is not valid (:class:`~miles.core.selection.SelectionError`).
        """
        chosen = self._canvas.selection if selected is None else list(selected)
        return self._canvas.doc.find(expression, chosen).tolist()

    def count(self, expression: str) -> int:
        """Return the number of atoms an expression picks."""
        return len(self.find(expression))

    def define(self, name: str, expression: str | None) -> None:
        """Name an expression, so later ones can use the name as a word (None removes the name).

        ``doc.define("pocket", "byres within 5 of ligand")`` lets you write
        ``doc.select("pocket and not hydrogen")``. The expression is worked
        out again each time the name is used, so it follows edits; session
        files keep the names.

        Raises
        ------
        ValueError
            If the name is not a plain word, is a keyword, or the expression is not valid.
        """
        self._canvas.doc.define(name, expression)

    @property
    def selections(self) -> dict[str, str]:
        """Named selections of this document: each name and its expression (a copy)."""
        return dict(self._canvas.doc.selections)

    def selection_with_hydrogens(self) -> list[int]:
        """Return the selected atoms with their hydrogens, or every atom when nothing is selected."""
        return self._canvas.selected_with_hydrogens()

    # -- results read with the structure -------------------------------------------------------

    @property
    def results(self) -> "Results":
        """Grids, normal modes, orbitals, charges and dipole read with the structure.

        They belong to the structure they were computed for: after an edit
        that changes the atoms or bonds they are stale, and this property
        gives empty results until undo returns to that structure
        (:meth:`stale_results` says what is waiting).
        """
        from miles.core.results import Results

        doc = self._canvas.doc
        return doc.results if not doc.stale_results() else Results(doc.mol.topology)

    def stale_results(self) -> list[str]:
        """Return what the results hold while they belong to another structure; empty when they apply."""
        return self._canvas.doc.stale_results()

    def add_volume(self, volume: "Volume") -> int:
        """Add a grid to the results as one undoable step; returns its index."""
        return self._canvas.doc.add_volume(volume)

    def keep_volume(self, volume: "Volume") -> int:
        """Keep a grid worked out from the results, such as an orbital, without an undo step.

        Use it for grids made for display from what the document already
        holds; grids the user adds or makes go through :meth:`add_volume`.

        Returns
        -------
        int
            The grid's index in ``results.volumes``.
        """
        return self._canvas.doc.keep_volume(volume)

    def remove_volume(self, index: int) -> "Volume":
        """Remove a grid from the results as one undoable step, and return it."""
        return self._canvas.doc.remove_volume(index)

    # -- frames -----------------------------------------------------------------------------

    @property
    def frames(self) -> np.ndarray | None:
        """Coordinates of every frame, shape ``(n_frames, n_atoms, 3)``, or None."""
        return self._canvas.doc.frames

    @property
    def frame_count(self) -> int:
        """Number of frames (1 without a trajectory)."""
        return self._canvas.doc.frame_count

    @property
    def frame_labels(self) -> list[str] | None:
        """Text of each frame, e.g. conformer energies or scan values."""
        return self._canvas.doc.frame_labels

    @property
    def frame_energies(self) -> list[float] | None:
        """Energy of each frame in kcal/mol, if known."""
        return self._canvas.doc.frame_energies

    def set_frames(self, frames: np.ndarray | None, labels: list[str] | None = None,
                   energies: list[float] | None = None) -> None:
        """Replace the trajectory as one undoable step and show its first frame.

        Parameters
        ----------
        frames : numpy.ndarray or None
            Coordinates, shape ``(n_frames, n_atoms, 3)``; None removes the frames.
        labels : list of str, optional
            Text of each frame.
        energies : list of float, optional
            Energy of each frame in kcal/mol (shown in the timeline plot).
        """
        from miles.ui.document import Change

        doc = self._canvas.doc
        with doc.changing(Change.COORDS):
            doc.set_frames(frames, labels, energies=energies)

    def replace_frame_coordinates(self, frames: np.ndarray) -> None:
        """Change the coordinates of every frame as one undoable step.

        Labels, energies and the frame shown stay as they are, which suits
        operations such as aligning or centring the frames.

        Parameters
        ----------
        frames : numpy.ndarray
            New coordinates, the same shape as :attr:`frames`.

        Raises
        ------
        ValueError
            If the shape differs from the current frames.
        """
        from miles.ui.document import Change

        doc = self._canvas.doc
        doc.stop_animation()
        own = doc.trajectory
        new = np.asarray(frames, dtype=float)
        if own is None or new.shape != own.array.shape:
            raise ValueError("The new coordinates must have the shape of the current frames")
        frame = doc.frame
        labels = None if own.labels is None else list(own.labels)
        energies = None if own.energies is None else list(own.energies)
        with doc.changing(Change.COORDS):
            doc.set_frames(new, labels, own.loop, energies)
            doc.set_frame(frame)

    def show_frames(self, frames: np.ndarray | None, labels: list[str] | None = None,
                    loop: bool = True) -> None:
        """Play other frames, such as a vibration, without changing the document.

        The animation takes the place of the document's frames in playback
        and the timeline until ``show_frames(None)``, an edit that moves
        atoms, or undo. It is view state: edits, undo steps and saved files
        never see its coordinates, and the document's own frames and
        geometry come back when it stops. :meth:`set_frames` changes the
        document's frames instead.

        Parameters
        ----------
        frames : numpy.ndarray or None
            Coordinates, shape ``(n_frames, n_atoms, 3)``; None stops the animation.
        labels : list of str, optional
            Text of each frame.
        loop : bool, optional
            The motion is periodic (the default), so playback blends the last frame into the first.
        """
        doc = self._canvas.doc
        if frames is None:
            doc.stop_animation()
        else:
            doc.show_animation(frames, labels, loop=loop)

    @property
    def animating(self) -> bool:
        """Whether frames from :meth:`show_frames` are shown instead of the document's own."""
        return self._canvas.doc.animating

    def show_geometry(self, coords: np.ndarray) -> None:
        """Show other coordinates without changing the document, e.g. a preview.

        The document's own coordinates come back with the next edit, undo,
        save or frame change; :meth:`edit` changes them.
        """
        self._canvas.doc.show_geometry(np.asarray(coords, dtype=float))

    @property
    def frame(self) -> int:
        """Index of the frame shown."""
        return self._canvas.doc.frame

    @frame.setter
    def frame(self, index: int) -> None:
        self._canvas.doc.set_frame(int(index))

    # -- view -------------------------------------------------------------------------------

    def style(self, **options: Any) -> None:
        """Change display options, as the View menu does.

        Parameters
        ----------
        **options
            Fields of :class:`~miles.render.scene.Style`, e.g.
            ``mode="tube"``, ``labels="symbol"``, ``show_hydrogens=False``,
            or the name of a plugin's toggle, e.g. ``hydrogen_bonds=True``,
            ``dipole=True``.

        Raises
        ------
        AttributeError
            If an option does not exist.
        """
        display = self._canvas.display
        manager = getattr(self._window, "plugins", None)
        toggles = manager.toggles if manager is not None else {}
        for name in options:
            if not hasattr(display, name) and name not in toggles:
                raise AttributeError(f"No display option {name!r}")
        converted = {}
        for name, value in options.items():
            if not hasattr(display, name):
                continue
            current = getattr(display, name)
            converted[name] = type(current)(value) if isinstance(current, str) else value
        if converted:
            self._canvas.set_style(**converted)
        for name, value in options.items():
            if not hasattr(display, name):
                toggles[name].set(self, bool(value))

    def fit(self, atoms: Iterable[int] | None = None) -> None:
        """Zoom and centre so the whole structure is in view, or zoom smoothly to some atoms."""
        if atoms is None:
            self._canvas.fit(animate=False)
        else:
            self._canvas.zoom_to(int(i) for i in atoms)

    def shown_atoms(self) -> np.ndarray:
        """Boolean mask of the atoms shown: drawn, or standing for their residue in a cartoon (CA or P)."""
        return self._canvas.shown_atoms()

    def drawn_atoms(self) -> np.ndarray:
        """Boolean mask of the atoms drawn as spheres (not hidden, not part of a cartoon)."""
        return self._canvas.drawn_atoms()

    def play(self) -> None:
        """Start playing the frames in a loop in the timeline (when this is the tab in front)."""
        if self._window.current_canvas() is self._canvas:
            self._window.stage.timeline.play()

    def turn(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None:
        """Turn the view about the screen axes, in degrees."""
        self._canvas.camera.turn(math.radians(x), math.radians(y), math.radians(z))
        self._canvas.redraw()

    @property
    def rotation(self) -> np.ndarray:
        """Orientation of the view, a 3 x 3 rotation matrix (a copy; assign a new one to turn the view)."""
        return np.array(self._canvas.camera.rotation, dtype=float)

    @rotation.setter
    def rotation(self, matrix: np.ndarray) -> None:
        value = np.asarray(matrix, dtype=float)
        if value.shape != (3, 3):
            raise ValueError("The rotation must be a 3 x 3 matrix")
        self._canvas.camera.rotation = value
        self._canvas.redraw()

    def atom_reps(self) -> np.ndarray:
        """Representation of each atom: its own, or the view's style.

        Returns
        -------
        numpy.ndarray
            One :class:`~miles.core.molecule.Representation` value per
            atom, such as ``"ball"`` or ``"cartoon"``; ``"hidden"`` for atoms
            hidden in focus mode or given that representation. Hydrogens
            hidden by **View > Hydrogens** keep theirs (see
            :meth:`drawn_atoms`).
        """
        from miles.render.scene import atom_reps

        return atom_reps(self._canvas.mol, self._canvas.display, self._canvas.doc.display)

    def render(self, scale: float = 1.0, background: str = "white",
               ambient_occlusion: bool = False) -> "QImage":
        """Draw the view as a picture, with what plugins add to it, as exported pictures are drawn.

        Parameters
        ----------
        scale : float, optional
            Pixels per point of the view on screen.
        background : {"white", "transparent", "theme"}, optional
            Background of the picture.
        ambient_occlusion : bool, optional
            Shade the atoms by how buried they are (OpenGL only).

        Returns
        -------
        QImage
            The picture, the size of the view times ``scale``.
        """
        from miles.render.theme import Background

        return self._canvas.render_export(scale, Background(background), ambient_occlusion)

    def export_image(self, path: str, scale: float = 2.0, background: str = "white") -> None:
        """Save a picture of the view.

        Parameters
        ----------
        path : str
            PNG file to write.
        scale : float, optional
            Pixels per screen pixel.
        background : {"white", "transparent", "theme"}, optional
            Background of the picture.

        Raises
        ------
        OSError
            If the file cannot be written.
        """
        from miles.render.theme import Background

        image = self._canvas.render_export(scale, Background(background), False)
        if not image.save(path):
            raise OSError(f"Could not write {path}")

    def image(self, scale: float = 1.0, background: str = "white", crop: bool = True) -> "Picture":
        """Return a picture of the view, which IPython and Jupyter show below the cell.

        Parameters
        ----------
        scale : float, optional
            Shown size relative to the view on screen. The picture has
            twice as many pixels, so it stays sharp on high-resolution
            screens.
        background : {"white", "transparent", "theme"}, optional
            Background of the picture.
        crop : bool, optional
            Trim the empty space around the structure, as **File > Export
            Image...** does.

        Returns
        -------
        Picture
            The picture; ``save(path)`` writes it as PNG.

        Raises
        ------
        RuntimeError
            If the view has no size yet (the window was never shown).
        """
        from miles.render.theme import Background
        from miles.ui.export_dialog import crop_image

        kind = Background(background)
        if self._canvas.width() < 1 or self._canvas.height() < 1:
            raise RuntimeError("The view has no size yet; show the window first")
        image = self._canvas.render_export(scale * Picture.DENSITY, kind, False)
        if crop:
            image = crop_image(image, kind == Background.TRANSPARENT)
        return Picture(image)

    def save(self, path: str | None = None) -> None:
        """Save the structure (to ``path``, or to the file it came from).

        The format follows the extension; ``.miles`` keeps everything, the
        view included.

        Raises
        ------
        OSError
            If saving fails.
        """
        doc = self._canvas.doc
        if path is not None:
            doc.path = path
        if not doc.path or not self._window.save_document(doc):
            raise OSError("Could not save; give a path")
