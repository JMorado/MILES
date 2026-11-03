The Molecule Model
##################

Everything in MILES works on one data structure,
:class:`~miles.core.molecule.Molecule`, defined in
``core/molecule.py``, held by a :class:`~miles.core.document.Document`
together with how its atoms are shown. This chapter describes how a
molecule stores its atoms and bonds, the rules for changing it, the
properties readers and plugins can add, and how documents keep display
state, undo history and change events.

Atoms are rows of columns
*************************

A molecule is a table. Each per-atom quantity is one NumPy array with one
entry per atom, so work on the whole structure (drawing, selecting,
measuring, writing files) is array arithmetic rather than a loop over
Python objects. The built-in columns are listed in
:data:`~miles.core.molecule.ATOM_COLUMNS`:

.. list-table::
   :header-rows: 1
   :widths: 18 18 64

   * - Column
     - Type
     - Holds
   * - ``number``
     - uint8
     - Atomic number; 0 is a dummy atom. ``mol.symbols`` and
       ``mol.column("symbol")`` give the element symbols.
   * - ``coords``
     - float, (n, 3)
     - Position in angstrom (:attr:`~miles.core.molecule.Molecule.coords`).
   * - ``geom``, ``nsites``, ``valence``
     - object, int16, float
     - How the atom is meant to bond: its coordination geometry (a
       :class:`~miles.core.geometry.Geometry` such as ``tetrahedral``), the
       number of bonds it should have and the bond-order sum it should
       reach. A carbonyl oxygen is ``linear`` with one site and valence 2.
       Unset values are ``None``, -1 and NaN.
   * - ``inferred``
     - bool
     - The three fields above were guessed from a file (see :doc:`chemistry`),
       so Add Hydrogens falls back to standard valences.
   * - ``formal_charge``
     - int8
     - Formal charge as a file or SMILES stated it; -128
       (:data:`~miles.core.molecule.UNSET_CHARGE`, ``None`` through an
       atom view) when nothing stated it. Unset charges are worked out from
       the bonds when needed (:func:`~miles.core.charges.formal_charges`).
   * - ``isotope``
     - uint16
     - Mass number of a specific isotope (13 for carbon-13); 0 for the
       natural mixture.
   * - ``name``, ``altloc``
     - str, str
     - Atom name (``CA``) and alternate location indicator (``A``, ``B``)
       from PDB-like files.
   * - ``resname``, ``resseq``, ``icode``, ``chain``, ``hetero``
     - str, int, str, str, bool
     - Residue data: residue name (``ALA``), number, insertion code (the
       ``A`` of residue ``52A``), chain and whether the record was
       ``HETATM``. Atoms built in MILES have empty residue fields, except
       hydrogens added to a residue, which take them from their parent.
   * - ``ss``
     - str
     - Secondary structure of the residue: ``H`` helix, ``E`` strand or empty.
   * - ``uid``
     - int64
     - An identifier unique in the running program, kept through edits and
       copies (:attr:`~miles.core.molecule.Molecule.uids`). New atoms get new
       ones, and so do atoms pasted or appended from another molecule.

Text columns hold Python strings; the values are interned, so a hundred
thousand ``"CA"`` names cost one string and a hundred thousand pointers.

There are three ways to read atoms:

.. code-block:: python

   mol.coords, mol.numbers, mol.column("resname")   # whole columns, read-only arrays
   mol.mask("resname", {"HOH", "WAT"})                # boolean mask of matching atoms
   mol.atoms[i].resname                               # one atom, through a view

``mol.atoms`` is a sequence of :class:`~miles.core.molecule.AtomView`
objects made on demand. A view reads and writes the columns at its row, so
``mol.atoms[i].symbol = "N"`` changes the molecule. Views are meant for
readable code that touches a few atoms; code that visits every atom should
read columns, which is tens of times faster. The free-standing
:class:`~miles.core.molecule.Atom` dataclass carries the same fields
between molecules: :meth:`~miles.core.molecule.Molecule.add_atom` takes
one (or a view, or a symbol), and
:meth:`~miles.core.molecule.Molecule.from_atoms` builds a molecule from a
list of them. The fastest way to make a large molecule is
:meth:`~miles.core.molecule.Molecule.from_columns`, which takes whole
columns.

Storage shared between copies
*****************************

:meth:`~miles.core.molecule.Molecule.copy` does not copy any array. The
copy and the original share every column until one of them changes it; the
one that changes a column first takes a private copy of that column only.
Appending atoms writes into spare room at the end of a column, which other
copies cannot see, so building atom by atom does not copy whole columns
either.

This makes copies cost the same for ten atoms as for a million, and it is
what undo relies on (see below). It comes with one rule: **never write into
an array a molecule handed out.** Arrays returned by ``coords``,
``numbers``, :meth:`~miles.core.molecule.Molecule.column`,
``bond_pairs`` and ``bond_orders`` are read-only views, so a write such as
``mol.coords[0] += 1`` raises an error instead of silently changing a
shared column. Change a molecule through its methods:

.. code-block:: python

   mol.coords = mol.coords + shift              # replace a whole column
   mol.set_positions(atoms, new_positions)      # move some atoms
   mol.set_column("ss", "H", helix_atoms)       # set a field for some atoms
   mol.atoms[i].name = "CB"                     # one field of one atom

An array handed out stays as it was: the molecule marks the column as
shared, so its next write copies the column first (the old array is the
snapshot, the molecule moves on). The atom identifiers in ``uid`` cannot
be set at all, so two atoms never share one.

Bonds
*****

Bonds are a second table: :attr:`~miles.core.molecule.Molecule.bond_pairs`
holds the atom pairs ``(i, j)`` with ``i < j`` and
:attr:`~miles.core.molecule.Molecule.bond_orders` their orders, in the
order the bonds were made. :attr:`~miles.core.molecule.Molecule.bonds`
presents the same data as a read-only mapping from pair to order, for code
that looks bonds up (``mol.bonds[(i, j)]``, ``(i, j) in mol.bonds``).

Neighbour lists are built from the table when first needed, as arrays in
compressed sparse rows: the neighbours of atom ``i``, in ascending order,
with the bond row of each. ``neighbors(i)``, ``degree(i)``, ``order_sum(i)``
and :meth:`~miles.core.molecule.Molecule.bond_order` read them, in time
proportional to the number of neighbours. The arrays never change, so
copies (and undo snapshots) share them for as long as their bonds are the
same. Bonds added one at a time go to a small overlay instead of a rebuild,
which keeps building atom by atom cheap; the arrays are rebuilt when the
overlay grows. :meth:`~miles.core.molecule.Molecule.adjacency` hands the
arrays out for vectorised work over the bond graph.

Change bonds through :meth:`~miles.core.molecule.Molecule.add_bond` (which
also changes the order of an existing bond),
:meth:`~miles.core.molecule.Molecule.remove_bond`,
:meth:`~miles.core.molecule.Molecule.set_bonds` and, for readers and other
bulk work, :meth:`~miles.core.molecule.Molecule.set_bond_arrays`.

Bond orders are floats. 1, 2 and 3 are single, double and triple bonds; 1.5
marks an aromatic bond, drawn as one solid and one dashed line; an order
between 0 and 1 is drawn as a dashed line and can be used for coordination or
partial bonds.

The topology token
==================

:attr:`~miles.core.molecule.Molecule.topology` is a number that changes
whenever atoms are added or removed, an element changes or a bond changes,
and only then. Copies share it while their atoms and bonds are the same.
That makes it the natural cache key for anything worked out from the bond
graph: rings, perceived chemistry, force-field types, partial charges.
Moving atoms or renaming them leaves it alone.

Properties
**********

Readers and plugins can give atoms data of their own as extra columns:

.. code-block:: python

   mol.set_property("b_factor", values)              # adds the column if needed
   mol.add_property("site", dtype=object, default="")
   mol.column("b_factor"), mol.properties, mol.has_property("site")
   mol.atoms[i].get("b_factor")

A property column behaves like a built-in one: it follows its atoms when
atoms are deleted or reordered, new atoms get the column's default (NaN for
floats), copies and undo steps share it, and session files save it. PDB
files give ``occupancy`` and ``b_factor`` (read from the records' fixed
columns, since Open Babel does not hand them back).

Some properties are worked out from the structure rather than read or set
by the user, partial charges being the main one. For these,
``mol.property_sources[name]`` records what they were worked out from, so a
cache can tell whether they still hold; copies and undo snapshots keep it
with the column, and sessions leave such columns out because they can be
worked out again.

Partial charges
===============

The charges of the document's chosen model (``Document.charge_model``:
from the file, Gasteiger or MMFF94) live in the ``"charge"`` column. They
are worked out the first time something asks for them
(:meth:`Document.charges <miles.core.document.Document.charges>`) and kept
with the source ``(topology, model)``; moving atoms keeps them, any change
of atoms or bonds makes the next request work them out again, and undo
brings back the charges that went with the old structure. Charges from the
file apply only while the topology is the one the file described.

Indices and their stability
***************************

Atoms are identified by their index, their row in the table. Selections,
measurements and highlights store indices, so the rules for when indices
change matter.

* Adding atoms appends them; existing indices do not change.
* Attaching a fragment re-uses the index of the hydrogen it replaces for the
  fragment's anchor atom (which then counts as a new atom with a new
  ``uid``) and appends the rest. A selection that contained the hydrogen now
  contains the new atom, which is usually what the user expects.
* Removing atoms renumbers everything after them.
  :meth:`~miles.core.molecule.Molecule.remove_atoms` returns the map from
  old to new indices (``-1`` for removed atoms) and updates the bonds,
  properties and pinned measurements itself. Code that keeps indices
  elsewhere, such as the canvas selection, must either use that map or drop
  what it kept.

Code that must follow atoms through arbitrary edits keeps their ``uid``
instead, as the display and the selection do: the document keeps the
selection as identifiers (:attr:`Document.selected
<miles.core.document.Document.selected>` turns them into indices), so it
stays on the same atoms when others are added or deleted, and never lands
on another atom after an undo.

Graph queries
*************

A few graph algorithms are used throughout the program.

**Connected parts.** :meth:`~miles.core.molecule.Molecule.component`
returns the atoms reachable from a start atom by breadth-first search. It
takes a set of blocked atoms that the search may not enter. Blocking one end
of a bond and starting from the other gives the part of the molecule on that
side of the bond, which is how bond lengths, angles and dihedrals decide what
moves (see :doc:`building`).
:meth:`~miles.core.molecule.Molecule.component_labels` labels every atom
with its connected molecule at once, on the bond arrays.

**Ring bonds.** :meth:`~miles.core.molecule.Molecule.ring_bonds` finds
the bonds that lie in a ring. A bond is in a ring exactly when it is not a
bridge of the bond graph, so the method finds all bridges with Tarjan's
depth-first search in time linear in the number of atoms and bonds, and
returns the remaining bonds. It is used, for example, to decide that a bond
with an "aromatic" length outside a ring is really a single bond.

**Residues and chains.** A residue is told apart by its chain, number,
insertion code and name, the columns of
:data:`~miles.core.molecule.RESIDUE_KEY`;
:meth:`~miles.core.molecule.Molecule.residue_key` gives them for one atom.
Every part of MILES that groups atoms into residues uses this key, so
residues ``52`` and ``52A`` stay apart everywhere.
:meth:`~miles.core.molecule.Molecule.residue_atoms` and
:meth:`~miles.core.molecule.Molecule.chain_atoms` return the atoms that
share an atom's residue or chain. For atoms without residue data they
return the atom's whole connected molecule, so double-click selection works
for built molecules too.

Pinned measurements
*******************

``Molecule.measurements`` is a list of 2-4 atom index tuples for distances,
angles and dihedrals that should stay on screen. The list is kept up to date
when atoms are removed and is saved in session files.

Display state
*************

How each atom is drawn is view state, not chemistry, so it is not in the
molecule. A document keeps it in an
:class:`~miles.core.display.AtomDisplay`: a representation override per
atom (ball and stick, tube, wire, space fill, cartoon or hidden; empty
means the view's style) and whether focus mode hides the atom.

An ``AtomDisplay`` never changes; ``with_rep`` and ``with_hidden`` return
new ones, which is what lets undo keep old ones for nothing. It is keyed by
atom ``uid``: after any edit the document calls
:meth:`~miles.core.display.AtomDisplay.aligned`, which carries each atom's
settings to its new index and gives new atoms the view's style, with two
exceptions that make edits look natural. A new hydrogen takes the
representation of the atom it is bonded to, and a new atom in a residue
that was already there takes that residue's representation, so a mutated
side chain looks like the rest of its residue.

The renderer receives the display next to the molecule
(:func:`~miles.render.scene.atom_reps`, :func:`~miles.render.scene.build_scene`,
:class:`~miles.render.picking.Picker`), and session files store it under
their own key. Plugins read it as ``doc.atom_display`` and change it with
``doc.set_representation``, ``doc.hide``, ``doc.show`` or
``doc.edit_display`` (see :doc:`plugins`).

Documents
*********

:class:`miles.core.document.Document` holds a molecule, its display, its
selection, its trajectory, the results read with it, its file and title,
the chosen charge model, the named selections and the undo history. It
needs no Qt, so scripts, tests and other programs can edit and undo
documents without a window; :class:`miles.ui.document.Document` is the
same class with Qt signals for the widgets.

A reader returns a molecule and a :class:`~miles.core.jobinfo.JobInfo`.
The document takes copies of what it needs and leaves the ``JobInfo`` as
it was, so the same reader result can make several documents.

Trajectories and playback
=========================

The frames live in a :class:`~miles.core.frames.Trajectory`, apart from
the molecule: one read-only array per frame, with optional labels and
energies. A trajectory is never changed in place; changing a frame makes a
new one that shares the other frames, so an undo step that keeps the old
trajectory costs one frame. The molecule shows one frame at a time. An
edit that moves atoms is kept in the frame shown; an edit that adds or
removes atoms drops the frames, and undo brings them back with the atoms.

What is on screen is not always the document's own state. Playback blends
two frames (:meth:`~miles.core.document.Document.set_position`), an
animation such as a vibration shows frames of its own
(:meth:`~miles.core.document.Document.show_animation`), and
:meth:`~miles.core.document.Document.show_geometry` previews other
coordinates. These are view state: before any edit, undo step or save the
document puts its own coordinates back, so a blend or an animation never
reaches the frames, the history or a file. An edit that moves atoms ends
an animation, since it was made for the old geometry.

Results
=======

Grids, normal modes, orbitals, charges and the dipole moment read with a
structure live in :class:`~miles.core.results.Results`, stamped with the
topology token of the molecule they belong to. After an edit that changes
the atoms or bonds they are stale: the document keeps them (so undo makes
them apply again) but nothing uses them, ``doc.info`` shows them empty, a
``results`` event is announced and the canvas says so. Saving a session
leaves stale results out, since nothing could make them apply again.
Grids are added and removed with
:meth:`~miles.core.document.Document.add_volume` and
:meth:`~miles.core.document.Document.remove_volume`, as undoable steps.

Undo
====

Every edit goes through :meth:`~miles.core.document.Document.edit` (or
:meth:`~miles.core.document.Document.changing`, its ``with`` form), which
records an undo step, calls the editing function on the molecule and
announces the change. If the function raises, the document is put back as
it was and the step forgotten before the error goes on, so no edit is left
half made. An undo step is a snapshot of the whole document: a copy of the
molecule and the display, selection, trajectory, frame shown, results,
charge model and named selections. Because copies share columns and the
other parts are immutable, a snapshot costs almost nothing when it is
taken, and the edit that follows copies only the columns it changes. Moving one atom of a 100 000-atom structure therefore
adds one coordinate array to the history and nothing else. The document
keeps up to :data:`~miles.core.document.UNDO_LIMIT` steps; with shared
columns this limits their number rather than their memory, which grows
with what the edits changed. Keeping snapshots instead of inverse
operations keeps every editing function free of undo logic.

Each step remembers what its edit changed, so undoing a move announces
``coords`` and keeps the caches that do not depend on positions. Undo and
redo keep the current selection on its atoms; only when none of them is
left do they bring back the selection the step had, so undoing the
deletion of selected atoms selects them again. Defining a named selection
and choosing a charge model are undoable steps too.

Change events
=============

Each change is announced as a :class:`~miles.core.document.ChangeEvent`
to the document's ``listeners`` (and as the Qt signals ``changed`` and
``changed_detail`` in the window). Its ``kind`` says what changed:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Kind
     - Meaning
   * - ``structure``
     - Atoms or bonds added or removed. An edit that changes atoms,
       elements or bonds counts as this whatever its editor said.
   * - ``bonds``
     - Bond orders only.
   * - ``coords``
     - Positions only.
   * - ``style``
     - The display (representations, hidden atoms) or measurements.
   * - ``properties``
     - Per-atom properties, such as the charges after another charge model is chosen.
   * - ``frame``
     - Another trajectory frame, a blend of two, or other coordinates are shown.
   * - ``frames``
     - The frames were replaced or dropped, or an animation started or stopped.
   * - ``names``
     - Named selections were defined or removed.
   * - ``results``
     - Grids were added or removed, or the results became stale or apply again.

``atoms`` and ``properties`` say which atoms and properties changed when
the code that made the change said so; ``Document.last_change`` keeps the
latest event for code that is told only the kind, such as plugin callbacks.
Undo and redo announce the kind of the step they undo or redo.

Two revision counters turn events into cache keys. ``revision`` changes
with every change. ``layout_revision`` changes with every change except
atoms moving (``coords`` edits and frames), named selections and results. The canvas keeps everything
that depends on the atoms, bonds and display but not on positions (the
representation of each atom, radii, colours, which residues form the
cartoon, the hit tester) under the layout revision, so playing a trajectory
or dragging atoms only works out new positions. The sequence strip does the
same.
