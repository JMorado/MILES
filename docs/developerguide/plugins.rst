.. _dg-plugins:

Writing Plugins
###############

A plugin adds commands, file formats, panels and reactions to events
without changing the code of MILES. It is plain Python: a module with a
``register(context)`` function, either a single file dropped into the
plugin folder or a package installed with ``pip`` that declares an entry
point. Plugins work through the same scripting interface as the Python
console (:class:`~miles.plugin.App` and
:class:`~miles.plugin.DocumentView`), so anything that can be typed in
the console can be made into a menu command. A plugin imports from two
modules: :mod:`miles` for the chemistry and :mod:`miles.plugin` for the
window (:ref:`dg-plugin-surface`).

Many of the built-in features are written this way, in
``miles/builtin_plugins`` (:ref:`dg-builtin-plugins`), using only those two
modules, like any other plugin.
Each is a good place to start, and ``examples/plugins/template`` is a
complete plugin package to copy (:ref:`dg-plugin-template`).

A first plugin
**************

Save this as ``hello.py`` in the plugin folder and start MILES:

.. code-block:: python

   """Count the atoms of the structure."""

   NAME = "Hello"
   VERSION = "0.1"


   def count_atoms(app):
       doc = app.current
       app.message(f"{doc.name} has {len(doc.molecule)} atoms")


   def register(context):
       context.add_command("Tools > Count Atoms", count_atoms, shortcut="Ctrl+Alt+Shift+C")

**Tools > Count Atoms** appears at the end of the Tools menu, in the
command palette (:kbd:`Shift Cmd P`) and in ``app.commands()`` in the
console, and ``app.run("Count Atoms")`` runs it from a script. Help >
Plugins... lists the plugin with its version and the commands it added.

The plugin folder is shown at the bottom of Help > Plugins..., and **Open
Plugin Folder** there creates it and opens it in the file browser. It is

- macOS: ``~/Library/Application Support/miles/MILES/plugins``
- Linux: ``~/.local/share/miles/MILES/plugins``
- Windows: ``%APPDATA%\miles\MILES\plugins``

The environment variable ``MILES_PLUGIN_DIR`` names another folder,
which is handy while developing a plugin. Every ``.py`` file in the folder
is a plugin, and so is every folder with an ``__init__.py`` (a package, for
a plugin of several modules); names starting with ``_`` or ``.`` are left
alone, so helper modules can live next to the plugins.

The module
**********

MILES imports the module and calls ``register(context)`` once, at
start-up, after the window is built. These module attributes are optional:

``NAME``
   Name shown in Help > Plugins (default: the file or entry point name).
``VERSION``
   Version shown next to the name.
``DESCRIPTION``
   One line describing the plugin (default: the first line of the module
   docstring).

``REQUIRES_API``
   Version of the plugin interface the plugin needs (see below).

``register`` should only set things up; work that takes time belongs in
the commands. Imports of heavy libraries can also go inside the commands,
so that start-up stays quick.

.. _dg-plugin-surface:
.. _dg-stability:

What to import, and what stays stable
=====================================

A plugin has two modules of MILES to import from:

``import miles``
   The chemistry library: molecules and their files, the selection
   language, geometry and neighbour searches, elements, the Open Babel
   bridge and force fields, building and perception, protein residues, and
   the data of quantum chemistry outputs (vibrations, orbitals, grids,
   charges). :doc:`../api/library` lists it by purpose.

``from miles.plugin import ...``
   Working with the window: the context and what it creates (``PluginContext``,
   ``Setting``, ``Toggle``, ``Choice``, ``PluginLayer``,
   ``BackgroundResults``), the application and its documents (``App``,
   ``DocumentView``, ``Severity`` of messages, the ``Change`` kinds, the
   display of atoms), the drawing types and ``color``, menu items and
   structure sources, panels (``FloatingPanel``, ``DocumentPanel``) and
   dialogs (``PluginDialog``), the theme, and text helpers such as
   ``plural`` and ``atom_name``. :doc:`../api/plugin` lists it by purpose.

``from miles.plugin.widgets import ...`` and ``from miles.plugin.render import ...``
   Two provisional submodules: the controls, colours and plot menus the
   built-in panels are made of (:doc:`../api/plugin_widgets`), and the
   meshes, arrows and surfaces they draw (:doc:`../api/plugin_render`).

Each name has one home: chemistry in :mod:`miles`, the window in
:mod:`miles.plugin`. A plugin's algorithms belong to the plugin: a
relaxed scan, a Ramachandran classification or a protonation table is part
of its feature, not of the library, and the plugin's own modules import
each other freely.

The promise is this:

* :mod:`miles` and :mod:`miles.plugin` are the stable surface, versioned
  together. A release that changes them in a way that breaks plugins raises
  the version of the plugin interface, ``PLUGIN_API``; a plugin that sets
  ``REQUIRES_API`` to a higher number fails to load with a message saying
  which version it needs.
* :mod:`miles.plugin.widgets`, :mod:`miles.plugin.render` and
  :mod:`miles.testing` are public but provisional: they follow the panels,
  the renderer and the window of MILES, and a release may change them
  without raising ``PLUGIN_API``. A plugin that uses them may need small
  changes; the release notes say which.
* Everything else, ``miles.core`` included, is internal and may change in
  any release. ``app.window`` and ``doc.canvas`` exist for experiments in the
  console; a plugin that uses them can break with any release.

The built-in plugins keep to this: ``tests/test_plugin_boundary.py``
parses them, the code they share and the examples, and fails on imports of
any other module of MILES (or of another plugin, relative imports
included), on names imported from :mod:`miles` that are not public, on uses
of ``.window`` and ``.canvas``, and on underscore names of objects that are
not their own. It also checks the other way round: every name of
:mod:`miles` is used by MILES itself, apart from a short list that scripts
need, so the library does not collect code that only one plugin uses.

What a plugin can add
*********************

The context passed to ``register`` is a
:class:`~miles.plugins.PluginContext`. Its ``app`` attribute is the
running :class:`~miles.api.App`.

Commands
========

:meth:`~miles.plugins.PluginContext.add_command` takes a menu path,
the function to call and, optionally, a shortcut, an ``enabled`` test and
the item to place the command after or before:

.. code-block:: python

   def has_selection(app):
       return app.current is not None and bool(app.current.selection)

   context.add_command("Tools > Charges > Sum Over Selection", sum_charges,
                       enabled=has_selection, after="Compute")

Menus and submenus in the path are created when missing; a new top-level
menu goes before Help, and a path with no menu (just ``"Sum"``) goes to a
**Plugins** menu. ``after`` names the item to follow and ``before`` the
item to precede; a command whose item is missing goes to the end of the
menu. Anchoring to one of the built-in items with ``before`` keeps a
group of commands in order even when some of the plugins that add them are
switched off. The callback receives the ``App``; the ``enabled`` test is
called when the menu opens, greys the item out when it returns False, and
is checked again when the command runs from a shortcut, the palette or
``app.run``. A shortcut that another command already has is not set, and
Help > Plugins says so.

Commands run on the interface thread. A command that takes more than a
second or so should say what it is doing with ``app.message``.

File formats
============

:meth:`~miles.plugins.PluginContext.add_file_reader` and
:meth:`~miles.plugins.PluginContext.add_file_writer` register functions
for extensions in the format registry
(the registry of file formats), in front of Open Babel:

.. code-block:: python

   def read_toy(path):          # -> (Molecule, JobInfo)
       ...

   def write_toy(path, molecule, info):
       ...

   context.add_file_reader("toy", read_toy)
   context.add_file_writer("toy", write_toy)

The extension is then listed in **File > Open...**, such files open from
there, by dropping them on the window and from the command line, and
**File > Save As...** offers the format. A reader raises :class:`~miles.files.common.FileFormatError` with a
message for the user when the file is not what it expects; any other
exception is turned into such an error naming the plugin, so a bug shows
the usual "Could not open file" message instead of an error dialog.
Readers that give a structure without bonds can call
:func:`miles.rebond <miles.core.perception.rebond>` to work them out from the
distances.

Panels
======

:meth:`~miles.plugins.PluginContext.add_panel` puts a widget in a
floating panel on the right of the window, below the built-in panels,
and adds a command that shows and hides it (**View > title** by default):

.. code-block:: python

   from PySide6.QtWidgets import QLabel

   label = QLabel()
   context.add_panel("Charges", lambda app: label, shortcut="Ctrl+Alt+Shift+Q")

The panel has a title row with a chevron that folds it and a **Close**
link, can be dragged elsewhere, and its place and fold are remembered like
those of the other panels. It starts hidden. The widget is built once, by
the function given, and takes the app's colours from the panel's
stylesheet.

Every panel opens and closes the same way: ``open_panel()``,
``close_panel()`` and ``toggle()``, and its ``opened`` and ``closed``
signals say when it appears and goes, for a panel that should only work
while it is seen.

A panel that needs more, such as buttons in its title row, can be a
:class:`~miles.plugin.FloatingPanel` subclass of the plugin's own.
:meth:`~miles.plugins.PluginContext.add_floating_panel` puts it in the same
column, or with ``place="bottom"`` above the timeline, under a name of the
plugin's choice, which is also the settings key of its place and fold; the
plugin adds the command that opens it:

.. code-block:: python

   panel = NotesPanel()
   context.add_floating_panel("notes", panel)
   context.add_command("Tools > Notes", lambda app: panel.toggle(), shortcut="Ctrl+Alt+N")

Most panels are about the tab in front: they show its data, follow it when
another tab comes to the front, are worked out again when it changes and
hide when it has nothing to show. A
:class:`~miles.plugin.DocumentPanel` does this bookkeeping, and
:meth:`~miles.plugins.PluginContext.add_document_panel` puts it in the
window with its command. The subclass says what the panel needs and how to
show it:

.. code-block:: python

   from miles.plugin import Change, DocumentPanel

   class ChargePanel(DocumentPanel):
       """Partial charges of the tab in front, in a QLabel (``self.text``) built in __init__."""

       refresh_on = frozenset({Change.STRUCTURE, Change.PROPERTIES})   # changes that matter
       unavailable = "No partial charges for this structure"

       def has_content(self, doc):
           return doc.charges() is not None

       def show_contents(self, doc):
           self.text.setText("  ".join(f"{q:+.2f}" for q in doc.charges()))

   panel = context.add_document_panel("charges", ChargePanel(), "Tools > Charge Table",
                                      shortcut="Ctrl+Alt+Q")

The command opens the panel, or closes it when it is open, and says
``unavailable`` when the tab in front has nothing to show. A panel that
should open by itself when a tab has something to show (like the
vibrations panel) sets ``wanted = True``; ``follows_selection`` makes it
refresh when the selection changes; ``closing()`` is called when the user
closes it and ``set_theme()`` when the appearance changes. Data the panel
keeps for each tab goes in ``store(doc)``, which is dropped with the tab,
and what it should remember in session files in ``saved(doc)``. The
vibrations, Ramachandran, properties and grid panels are made this way.

Dialogs
=======

A :class:`~miles.plugin.PluginDialog` is a dialog in the colours of the
interface, owned by the window, with the margins and buttons of the
built-in dialogs. Put the fields in ``body`` and add the row of buttons
with ``buttons()``; ``run()`` shows the dialog, returns whether it was
accepted and deletes it afterwards, so read the fields right after it:

.. code-block:: python

   from PySide6.QtWidgets import QSpinBox
   from miles.plugin import PluginDialog

   def ask_count(app):
       dialog = PluginDialog(app, "Repeat")
       steps = QSpinBox()
       dialog.body.addWidget(steps)
       dialog.buttons("Repeat")
       if dialog.run():
           repeat(app.current, steps.value())

A dialog that checks what was typed overrides ``accept`` and calls
``super().accept()`` only when it is right, as Build from Sequence does.

Events
======

These methods register functions that receive a
:class:`~miles.api.DocumentView`:

:meth:`~miles.plugins.PluginContext.on_document_opened`
   a new tab: a file opened, a structure downloaded, a new document;
:meth:`~miles.plugins.PluginContext.on_current_changed`
   another tab came to the front;
:meth:`~miles.plugins.PluginContext.on_selection_changed`
   the selection of a tab changed;
:meth:`~miles.plugins.PluginContext.on_structure_changed`
   atoms, bonds or coordinates of a tab changed (also while atoms are
   dragged, so keep the function quick);
:meth:`~miles.plugins.PluginContext.on_frame_changed`
   a tab shows another frame of its trajectory, or its frames were
   replaced or removed (many times a second during playback);
:meth:`~miles.plugins.PluginContext.on_document_changed`
   anything changed, with a second argument, a :class:`~miles.plugin.Change`,
   saying what:

   ``Change.STRUCTURE``
      atoms added or removed (or a new molecule);
   ``Change.BONDS``
      bonds or bond orders;
   ``Change.COORDS``
      positions only;
   ``Change.STYLE``
      representations, hidden atoms or display options;
   ``Change.PROPERTIES``
      per-atom properties such as charges, names or residues, and the
      charge model;
   ``Change.FRAME``
      another frame of the trajectory is shown;
   ``Change.FRAMES``
      the frames were replaced or removed;
   ``Change.NAMES``
      named selections;
   ``Change.RESULTS``
      grids added or removed, or results that became stale or apply again.

   The two events above are filtered forms of this one; a panel uses it to
   skip work that a kind of change does not affect. A ``Change`` compares
   equal to its text (``"coords"``), so older plugins keep working.

:meth:`~miles.plugins.PluginContext.on_theme_changed` registers a
function that receives the :class:`~miles.render.theme.Theme`. It is
called once straight away, so a panel that draws its own colours is right
from the start, and again when the appearance changes:

.. code-block:: python

   context.on_theme_changed(plot.set_theme)
   context.on_frame_changed(lambda doc: plot.refresh() if doc == context.app.current else None)

:meth:`~miles.plugins.PluginContext.on_document_closed` registers a
function called just before a tab closes, with the tab still open, for a
plugin that keeps things of its own for each tab outside
:meth:`~miles.plugins.PluginContext.document_cache`.
:meth:`~miles.plugins.PluginContext.on_unload` registers a function called
when the plugin is switched off or reloaded, for what MILES does not track
itself, such as a timer or a thread of the plugin's own; everything added
through the context is taken out without it.

The document passed may not be the current one; compare it with
``context.app.current`` (or the document a panel follows) when only the
visible tab matters. A tab has one :class:`~miles.api.DocumentView` for as
long as it is open, so the handles compare equal and can be dictionary
keys. The first, empty tab exists before plugins load, so a plugin that
keeps something per document should also look at ``context.app.documents``
in ``register``.

Menus and toggles
=================

:meth:`~miles.plugins.PluginContext.add_menu` creates a submenu at a
chosen place and :meth:`~miles.plugins.PluginContext.add_separator`
draws a line in a menu. ``after`` and ``before`` take an item's title, or
the name of one of the separators that mark where the built-in plugins
put their commands: ``"scripts"`` in File (before Run Script...), and
``"protein tools"`` and ``"tool choice"`` in Tools. A tuple of anchors is
tried in turn, so a command can follow another plugin's command and still
find its place when that plugin is switched off.

:meth:`~miles.plugins.PluginContext.add_toggle` adds a checkable
command whose state belongs to each tab, like **View > Hydrogen Bonds**.
The check mark follows the tab in front, the state is saved in session
files, and scripts set it by name with ``doc.style``:

.. code-block:: python

   def switched(doc, on):
       layer.update(doc)        # draw again with the new setting

   axes = context.add_toggle("View > Principal Axes", "principal_axes", switched,
                             after="Hydrogen Bonds")
   ...
   if axes.value(doc):           # in a drawing function
       ...

``default`` gives the state of a tab that was never switched, or a function
working it out. Sessions written before a feature became a plugin are read
too: a toggle whose state was saved with the view options under the same
name takes it from there.

:meth:`~miles.plugins.PluginContext.add_choice` adds a group of checkable
commands of which one is chosen per tab, like the colourings in
**View > Surface**. Its value is saved like a toggle's, ``changed`` is
called with the tab and the new value, ``enabled`` greys out options that
do not apply, and ``read`` and ``write`` keep the value elsewhere, as the
charge model of **View > Charges** lives in the document:

.. code-block:: python

   coloring = context.add_choice("View > Shading", "shading",
                                 {"flat": "Flat", "smooth": "Smooth"},
                                 lambda doc, value: layer.update(doc), default="smooth")
   ...
   if coloring.value(doc) == "flat":   # in a drawing function
       ...

Drawing over the molecule
=========================

:meth:`~miles.plugins.PluginContext.add_annotations` registers a
function that returns :class:`~miles.render.layers.Annotations`: lines
and text labels anchored to atoms (an index) or to points in angstrom.
MILES projects them and draws them on screen, onto exported pictures
and onto every movie frame:

.. code-block:: python

   from miles.plugin import Annotations, Label, Line

   def contacts(doc, drawing):
       if not doc.selection:
           return None
       first = doc.selection[0]
       others = [i for i in doc.selection[1:]]
       return Annotations(
           lines=[Line(first, i, "#3A9BDC", dashes=(4, 2), drawn_only=True) for i in others],
           labels=[Label(f"{i + 1}", i, offset=(0, -12)) for i in others])

   layer = context.add_annotations(contacts)

The function is called again when the structure, the selection, the frame
shown, the display options or the theme change, and when the plugin calls
the returned layer's :meth:`~miles.plugins.PluginLayer.update`; otherwise
its result is reused, so drawing stays fast. Colours can be given as a name
(``"#3A9BDC"``, ``"red"``), a ``QColor`` or a tuple, of integers from 0 to
255 or of floats from 0 to 1; :func:`~miles.plugin.color` reads any of
them, and :func:`~miles.plugin.atom_name` gives the label of an atom with
its residue (``A/ASP25:OD1``). ``drawing.export`` is True for
pictures and movies, where results must be exact; ``drawing.theme`` has the
colours. A line with ``drawn_only`` (or a label with it) is left out when an
atom it is anchored to is not shown; a line with ``ribbon`` ends on the
cartoon ribbon when its atom is a backbone atom drawn only as part of a
cartoon, as hydrogen bonds do. ``Label.away_from`` puts a label just past
the end of a line on screen, whatever the view, as the size of the dipole
moment is placed beyond the tip of its arrow.
:class:`~miles.render.layers.Arc` marks an angle at a vertex, the
shorter way round. Lines, labels and arcs take a ``tip``, shown when the
cursor rests on them.

Drawing in 3D
=============

:meth:`~miles.plugins.PluginContext.add_scene_items` registers a
function that returns :class:`~miles.render.layers.SceneItems`, drawn
with the molecule by both renderers, so atoms in front hide them:

* ``arrows``: rows from :func:`~miles.plugin.render.arrow_rows`, such as
  the dipole and the displacements of a vibrational mode;
* ``meshes``: :class:`~miles.render.layers.MeshItem` objects, triangle
  meshes with an opacity, drawn like the molecular surface (see-through
  below 1), such as the isosurfaces of grids and orbitals;
* ``atom_colors``: an RGB colour per atom in place of the element colours,
  as Color by Charge does.

.. code-block:: python

   from miles.plugin import SceneItems
   from miles.plugin.render import arrow_rows

   def pointer(doc, drawing):
       mol = doc.molecule
       if not len(mol):
           return None
       start = mol.center_of_mass()
       return SceneItems(arrows=arrow_rows([start], [start + (0, 0, 3)], 0.08,
                                           drawing.theme.accent))

   context.add_scene_items(pointer)

The same caching applies. Return the same mesh object while it has not
changed (keep it in a dictionary keyed by what it depends on), so it is not
sent to the graphics card again. A mesh that takes long to work out
belongs on a worker thread (see `Slow work`_). While atoms are dragged,
``drawing.dragging`` is True and the function is called again when the
drag ends, so a slow mesh can be kept until then, as the molecular
surface is.

Status text
===========

:meth:`~miles.plugins.PluginContext.add_status` adds text to the
status line after the formula, mass and point group, as the dipole moment
and the hydrogen bond count do. The function is called a moment after the
structure changes, not during atom drags or playback:

.. code-block:: python

   context.add_status(lambda doc: f"{len(doc.selection)} selected" if doc.selection else None)

Right-click menus
=================

:meth:`~miles.plugins.PluginContext.add_context_items` adds items to
the canvas's right-click menu. The function receives the document and a
:class:`~miles.ui.canvas.extensions.MenuTarget` (the atom or bond
clicked) and returns :class:`~miles.ui.canvas.extensions.MenuItem`
objects; an item with ``items`` is a submenu. ``section`` puts them after
the items for the atom (``"atom"``), the bond (``"bond"``) or the selection
(``"selection"``):

.. code-block:: python

   from miles.plugin import MenuItem

   def items(doc, target):
       if target.atom is None:
           return []
       return [MenuItem("Select Neighbours",
                        lambda: doc.select(doc.molecule.neighbors(target.atom), add=True))]

   context.add_context_items(items, "atom")

Structure sources
=================

:meth:`~miles.plugins.PluginContext.add_structure_source` adds a
:class:`~miles.ui.lookup.StructureSource`: somewhere structures can be
fetched from by a code or a name. Sources serve the search box of the
fragment window, names and codes given on the command line,
``app.look_up(text)`` and ``app.fetch(text)``. The ``fetch`` function runs
on a worker thread and returns a molecule and its
:class:`~miles.files.common.JobInfo`:

.. code-block:: python

   from miles.plugin import StructureSource

   def fetch_entry(code):            # e.g. a download; raise FileFormatError if not found
       ...
       return molecule, JobInfo(title=code.upper())

   context.add_structure_source(StructureSource(
       "Toy Database", fetch_entry, matches=lambda text: text.lower().startswith("toy"),
       what="toy code", example="toy42"))

A source with ``matches`` handles the texts it accepts before SMILES is
tried; a source without it is the fallback for text that is neither, like
compound names. ``what`` and ``example`` appear in the search box's hint.
A source that gives whole structures (``whole``, the default) opens a new
tab; otherwise the result becomes the fragment to place. A command such
as **File > Open from PDB...** is an ordinary command that asks for the
text and calls ``app.look_up``.

Fragment libraries
==================

:meth:`~miles.plugins.PluginContext.add_fragments` adds the fragments
of a JSON file in the format of the built-in library
(``miles/core/data/fragments.json``): categories that map to lists of
fragments, each with a ``key``, a ``name``, ``symbols``, ``coords`` in
angstrom, ``bonds`` as ``[i, j, order]``, the ``attach`` hydrogen (or
null) and optionally PDB atom ``names``. The categories ``"bio"`` and
``"ligand"`` add to the Biomolecules and Ligands tabs of the fragment
window; any other category becomes a tab of that title:

.. code-block:: python

   from pathlib import Path

   context.add_fragments(Path(__file__).with_name("sugars.json"))

Users can do the same without a plugin: ``.json`` files in the
``fragments`` folder next to the plugin folder are read at start-up.

Settings
========

:meth:`~miles.plugins.PluginContext.add_setting` adds a setting that
users change in the plugin's section of the Settings window
(**Settings...** in Help > Plugins opens it there): a check box for a bool, a number field for an int or a
float (with a range, a step and a unit), a text field for a string. Values are kept with the other
settings, and ``changed`` is called after a change:

.. code-block:: python

   cutoff = context.add_setting("cutoff", "Largest distance", 3.5, minimum=2.5, maximum=4.5,
                                step=0.1, unit="angstrom", changed=lambda value: layer.update())
   ...
   found = search(doc.molecule, cutoff.value)

Hydrogen Bonds (distance and angle), Binding Pocket (radius) and Charges
and Dipole (the length of the dipole arrow) have settings.
:meth:`~miles.plugins.PluginContext.setting_value` reads a setting of
another plugin, with a default for when that plugin is switched off: the
Binding Pocket finds hydrogen bonds with the criteria set for Hydrogen
Bonds, so both show the same bonds.

Slow work
=========

A command that takes more than a moment should not freeze the window.
:meth:`App.run_in_background <miles.api.App.run_in_background>` runs a
function on a worker thread; a bar at the top right says what runs and has
a **Cancel** link. Give the function copies of what it needs, since it
must not touch the window or the documents, and apply the result in
``done``, which runs on the interface thread. Check there that the document
did not change meanwhile:

.. code-block:: python

   def optimise(app):
       doc = app.current
       revision = doc.revision

       def done(mol):
           if doc in app.documents and doc.revision == revision:
               doc.set_molecule(mol)          # one undoable step
           else:
               app.message("The structure changed meanwhile; run it again")

       app.run_in_background(relax, doc.molecule.copy(), message="Optimising", done=done)

Generate Conformers, Scan Coordinate, Build from Sequence and Prepare
Protein work this way. ``app.wait()`` waits for the background work to
end, for scripts and tests.

A drawing or a panel whose result takes long, such as a molecular surface
or the descriptors of a large structure, asks for it with
:meth:`~miles.plugins.PluginContext.background_results` instead of
computing it. ``get`` returns the result for a key that says what it
depends on; a new key starts the work on a worker thread, and meanwhile the
previous result (or None) is drawn. When it arrives, the layer is drawn
again. Pictures and movie frames pass ``now=True`` and get the exact result
at once:

.. code-block:: python

   def scene(doc, drawing):
       key = (doc.revision, drawing.theme.name)
       mesh = meshes.get(doc, key, slow_mesh, doc.molecule.copy(), now=drawing.export)
       return SceneItems(meshes=[MeshItem(mesh, 1.0)]) if mesh is not None else None

   layer = context.add_scene_items(scene)
   meshes = context.background_results(layer, "Computing the mesh")

Only the latest request of a tab is worked out, and a cancelled one is not
started again until something changes. A panel passes a function instead
of a layer, called with the tab when its result arrives.

Keeping data for each tab
=========================

:meth:`~miles.plugins.PluginContext.document_state` returns a
dictionary the plugin keeps for a tab. It is saved in session files and
filled again when a session is opened, so a plugin can remember what a tab
showed; values must be JSON-compatible. Toggles and choices keep their
states there. The Orbitals and Grids plugin keeps the isosurface of each
tab this way. A :class:`~miles.plugin.DocumentPanel` reaches both
dictionaries as ``saved(doc)`` and ``store(doc)``.

:meth:`~miles.plugins.PluginContext.document_cache` returns another
dictionary for working data that is not saved: meshes, search results,
arrays. It may hold anything, and it is dropped when the tab closes and
when the plugin is switched off or reloaded, so the plugin needs no
bookkeeping of its own:

.. code-block:: python

   def scene(doc, drawing):
       cache = context.document_cache(doc)
       if cache.get("revision") != doc.revision:
           cache["mesh"], cache["revision"] = expensive_mesh(doc.molecule), doc.revision
       return SceneItems(meshes=[MeshItem(cache["mesh"], 1.0)])

Working with documents
**********************

The :class:`~miles.api.App` and :class:`~miles.api.DocumentView`
classes are the whole interface a plugin needs in most cases:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - To
     - Use
   * - read the structure
     - ``doc.molecule`` (a :class:`Molecule <miles.core.molecule.Molecule>`; read
       whole columns such as ``mol.coords`` and ``mol.column("resname")``
       for speed), ``doc.info``, ``doc.selection``
   * - keep data of your own on each atom
     - ``mol.set_property("score", values)``, ``mol.column("score")``: a
       property column follows its atoms through edits and is saved in
       sessions (see :doc:`molecule_model`)
   * - partial charges
     - ``doc.charges()``, ``doc.charge_model``, ``doc.set_charge_model(...)``,
       ``doc.charge_error()``, ``doc.file_charges_apply()``
   * - change it, undoably
     - ``doc.edit(func, ...)`` or ``with doc.editing() as mol:``;
       ``doc.edit(func, ..., kind=Change.COORDS)`` when only positions
       change, so that nothing depending on the bonds is worked out again
   * - select atoms
     - ``doc.select(atoms)``, ``doc.select(smarts="[OX2H]")``
   * - work with a trajectory
     - ``doc.frames``, ``doc.frame``, ``doc.set_frames(...)``,
       ``doc.replace_frame_coordinates(...)``; ``doc.show_frames(...)``
       plays an animation and ``doc.show_frames(None)`` stops it, and
       ``doc.show_geometry(coords)`` previews coordinates: both are view
       state that edits, undo and saving never see
   * - use the results read with the structure
     - ``doc.results`` (grids, modes, orbitals, charges, dipole; empty while
       stale), ``doc.stale_results()``; ``doc.add_volume(grid)`` and
       ``doc.remove_volume(index)`` for grids the user adds or removes (one
       undoable step each), ``doc.keep_volume(grid)`` for a grid worked out
       for display, such as an orbital (no undo step)
   * - change the look
     - ``doc.style(mode="tube", labels="symbol")``, ``doc.display``; per
       atom ``doc.set_representation("tube", atoms)``, ``doc.hide(atoms)``,
       ``doc.show()``, ``doc.edit_display(func)``, and ``doc.atom_display``
       to read it (``doc.edit_style(func)`` still works for code that sets
       ``mol.atoms[i].rep`` and ``.hidden``)
   * - know what changed
     - ``doc.last_change`` in an ``on_document_changed`` callback: the
       kind, and the atoms and properties when known
   * - move the view
     - ``doc.fit()``, ``doc.fit(atoms)``, ``doc.turn(y=30)``,
       ``doc.rotation`` (read it, or assign a rotation matrix)
   * - know what is drawn
     - ``doc.shown_atoms()``, ``doc.drawn_atoms()``, ``doc.atom_reps()``,
       ``doc.revision`` (raised by every change) and ``doc.layout_revision``
       (raised only when atoms or bonds change, not when they move): keys
       for cached results
   * - make a picture of the view
     - ``doc.render(scale, background)`` (a ``QImage``, with what plugins
       draw), ``doc.export_image(path)``, ``doc.image()``
   * - play the frames
     - ``doc.play()``
   * - make a new tab
     - ``app.new(molecule, "title")``, ``app.open(path)``,
       ``app.fetch("1crn")``
   * - ask for a file or a folder
     - ``app.ask_save_path(title, name, filters)``,
       ``app.ask_open_path(title, filters)``, ``app.ask_folder(title)``
   * - tell the user
     - ``app.message(text)``, ``app.message(text, severity=Severity.WARNING)``
       for a note that stays longer, ``severity=Severity.ERROR`` (or
       ``app.warn(title, text)``) for a dialog
   * - run another command
     - ``app.run("Edit > Clean Geometry")``
   * - build a dialog in the app's colours
     - :class:`~miles.plugin.PluginDialog` (see `Dialogs`_); ``app.theme``
       with ``theme_palette`` from :mod:`miles.plugin.widgets` for widgets
       drawn natively
   * - copy or export a plot
     - ``app.export_plot(panel, kind)``, as the right-click menus of the
       built-in plots do

Everything in :mod:`miles` works on the molecules these objects hand
out.

Errors
******

A plugin cannot stop MILES from starting. If importing the module or
calling ``register`` raises, the plugin is marked *failed*, the traceback
is printed to the terminal and kept, and the next plugin is loaded. If a
command, panel factory or event function raises later, the error dialog
shows the traceback once (the same error repeated is not shown again) and
Help > Plugins keeps the last one. Commands added before ``register``
failed stay in the menus.

Help > Plugins... lists every plugin found with its version, where it came
from and its state: *loaded*, *failed*, *switched off* or *not loaded
(--no-plugins)*. Selecting a row shows the plugin's description, its
commands and its last error. The check box in the first column switches a
plugin off or on at once; the choice is stored in the settings under
``plugins/disabled`` for later starts. MILES records everything a
plugin adds while it registers (commands, menus, toggles, settings, panels,
drawing layers, menu items, sources, fragments, file formats, event
functions), so switching it off takes all of it out again. **Reload**
unloads the selected plugin and loads it again, with its shortcuts; a
plugin and its own modules are imported anew, so edits to any of its files
count. Notes that are not errors, such as a shortcut already used by
another command, are kept apart in ``PluginInfo.warnings``. A panel added again goes to the end of its column
until the next start.

A plugin may set ``REQUIRES_API`` to the version of the plugin interface it
needs (``PLUGIN_API`` of :mod:`miles.plugin`, now 1). With an older
MILES it is then marked failed with a message saying so, instead of
failing later in a way that is hard to understand.

``miles --no-plugins`` (or the environment variable
``MILES_NO_PLUGINS``) starts without plugins from packages and from
the plugin folder, which helps when a plugin keeps MILES from working.
The built-in plugins still load; switch those off in Help > Plugins.

Packaging a plugin
******************

A plugin to share is best made a package that ``pip`` installs. Declare the
module in the ``miles.plugins`` entry point group of its
``pyproject.toml``:

.. code-block:: toml

   [project]
   name = "miles-charges"
   version = "1.0"

   [project.entry-points."miles.plugins"]
   charges = "miles_charges"

The value names the module that has ``register``, or the function itself
(``"miles_charges.plugin:register"``). Installed into the same
environment as MILES, the plugin loads at the next start; Help >
Plugins shows the package and its version as the source. MILES
declares its own plugins the same way in its ``pyproject.toml``.

Plugins load in a fixed order: the built-in ones, then entry points sorted
by name, then the plugin folder's files sorted by name. A later plugin that
registers a reader for an extension replaces the earlier one.

.. _dg-plugin-testing:

Testing a plugin
****************

:mod:`miles.testing` gives a plugin's tests a real MILES window that
runs without a screen. Load its fixtures with one line in the
``conftest.py`` at the top of the plugin's folder:

.. code-block:: python

   pytest_plugins = ["miles.testing"]

The ``miles_app`` fixture is the :class:`~miles.plugin.App` of a new
window with an empty tab; it has settings of its own and is closed after
the test without asking to save. :func:`~miles.testing.load_plugin`
registers the plugin's module in it as MILES does at start-up, and fails
the test with the traceback if ``register`` raises:

.. code-block:: python

   import my_plugin
   from miles.testing import document_from_smiles, load_plugin

   def test_count_atoms(miles_app):
       info = load_plugin(miles_app, my_plugin)
       document_from_smiles(miles_app, "CCO", "ethanol")
       miles_app.run("Tools > Count Atoms")
       assert not info.error        # no callback of the plugin failed

:func:`~miles.testing.wait_until` keeps the window working until a
condition holds, and ``miles_app.wait()`` until the work in the background
has ended. The other helpers do what a user would:

* :func:`~miles.testing.document_from_smiles` and
  :func:`~miles.testing.document_from_text` open a structure in a new tab;
* :func:`~miles.testing.switch_off`, :func:`~miles.testing.switch_on` and
  :func:`~miles.testing.reload_plugin` do what Help > Plugins does, to test
  that the plugin cleans up after itself;
* :func:`~miles.testing.plugin_setting` returns one of the plugin's
  settings, to read it or ``set`` it;
* :func:`~miles.testing.captured_messages` records what the window tells
  the user, errors included, without dialogs waiting for a click;
* :func:`~miles.testing.drawings` returns what the plugin's layers draw for
  a tab, with ``export=True`` as for a picture.

:mod:`miles.testing` is provisional (see :ref:`dg-stability`). Importing the module
prepares the environment: Qt draws offscreen and the canvas with QPainter,
and the plugin, fragment and autosave folders of daily use are replaced by
empty temporary ones, so tests neither see the plugins installed for daily
use nor leave documents to restore. Only the built-in plugins load besides
the ones a test adds.

.. _dg-plugin-template:

From a template to an installed plugin
**************************************

``examples/plugins/template`` is a complete plugin package, the Atom
Counter: **Tools > Count Atoms** and a toggle that keeps the count in the
status line. It has the package, its ``pyproject.toml`` with the entry
point and a test folder:

.. code-block:: text

   template/
       pyproject.toml              entry point in the miles.plugins group
       conftest.py                 pytest_plugins = ["miles.testing"]
       miles_atom_counter/
           __init__.py             NAME, VERSION, REQUIRES_API, register(context)
       tests/
           test_atom_counter.py

To make a plugin of your own from it:

1. Copy the folder and rename the package (``miles_atom_counter``) and the
   project in ``pyproject.toml``; the entry point names the package.
2. Write the plugin in the package's ``__init__.py``, importing from
   :mod:`miles.plugin`, and its tests in ``tests/``.
3. Run the tests in the environment MILES is installed in, from the
   plugin's folder: ``python -m pytest``. They need no screen.
4. Install the plugin into that environment with
   ``pip install -e . --no-deps`` (``-e`` keeps it pointing at your
   folder, so edits count after **Reload** in Help > Plugins).
5. Start MILES: Help > Plugins lists the plugin with the package as its
   source, and its commands are in the menus.

While writing it, the plugin folder is quicker still: point
``MILES_PLUGIN_DIR`` at a folder holding the package, or link the package
into the plugin folder, and use **Reload** after each change.
``tests/test_plugin_template.py`` runs the template's tests as part of the
test suite of MILES, so the template always works with the current
version.

An example
**********

``examples/plugins/selection_info.py`` in the repository uses the basic
parts of the interface: a panel with the formula, mass, centre of mass and radius
of gyration of the selection that follows the selection, the structure and
the current tab; two undoable commands in a Tools submenu; and a reader and
writer for ``.atoms`` files of ``symbol x y z`` lines.

.. _dg-builtin-plugins:

Built-in plugins
****************

These features of MILES are plugins in ``miles/builtin_plugins``,
written against the same interface as any other: they import only
:mod:`miles` and :mod:`miles.plugin` (besides their own modules, the
shared package below, the standard library, NumPy, PySide6 and Open
Babel), which ``tests/test_plugin_boundary.py`` checks. Each owns its
whole feature, interface and algorithms; code that several of them use
and the rest of MILES does not lives in ``miles/builtin_plugins/shared``,
which is not a plugin and not public. Plugins never import each other. They load first, in the
order of :data:`~miles.plugins.discovery.BUILTIN`, also with ``--no-plugins``,
and each can be switched off in Help > Plugins like a third-party plugin;
the core of MILES does not import them.

Each is self-contained: a module, or a package whose ``__init__.py``
registers it and whose other modules hold the panels, dialogs and
computations that only it uses (for example ``ramachandran/panel.py`` and
``ramachandran/angles.py``). Code that the core needs stays in the core:
partial charges (the charge colours use them), and vibrational modes,
orbitals and grids, which the file readers produce. Code that several
plugins need and the core does not is in the shared package: what a
hydrogen bond is (``shared/hbonds.py``, for Hydrogen Bonds and Binding
Pocket), side-chain replacement (``shared/mutation.py``) and what counts as
a small molecule (``shared/molecules.py``).

They solve the problems every plugin meets with the interface above:

* panels about the tab in front are :class:`~miles.plugin.DocumentPanel`
  objects (vibrations, Ramachandran, properties, grids), each with one
  command that opens and closes it;
* slow results come from :meth:`~miles.plugins.PluginContext.background_results`
  (the molecular surface, and the properties of large structures) or
  :meth:`App.run_in_background <miles.api.App.run_in_background>`
  (conformers, scans, sequences, protein preparation, orbital grids);
* dialogs are :class:`~miles.plugin.PluginDialog` objects that are deleted
  after use (Select, Open from PDB, Build from Sequence, Scan, Superpose,
  Prepare Protein, Export Movie, Combine Grids);
* exclusive options are choices (the surface colouring and the charge
  model), on/off options toggles;
* results that depend only on the bonds are keyed by
  ``doc.layout_revision`` (hydrogen bond donors and acceptors), so moving
  atoms does not work them out again.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Plugin
     - Adds
     - Shows how to
   * - ``drawing/``
     - File > Export 2D Drawing..., Edit > Copy 2D Drawing
     - ask for a file name, use the clipboard
   * - ``movie/``
     - File > Export Movie...
     - draw every frame with ``doc.render``, turning the view with
       ``doc.rotation``
   * - ``structure_lookup/``
     - File > Open from PDB..., PDB codes and PubChem names in the fragment
       search and on the command line
     - add structure sources
   * - ``select_atoms/``
     - Edit > Select...
     - open a dialog in the app's colours, select atoms
   * - ``vibrations/``
     - Tools > Vibrations, the vibrations panel, mode arrows
     - a document panel that opens by itself, plots that export, arrows in
       3D
   * - ``sequence/``
     - Tools > Build from Sequence...
     - open a structure in a new tab
   * - ``conformers``
     - Tools > Generate Conformers
     - replace the frames as one undoable step
   * - ``scan/``
     - Tools > Scan Coordinate...
     - frames with labels and energies for the timeline
   * - ``align_frames``
     - Tools > Align Frames
     - change the coordinates of every frame
   * - ``superpose/``
     - Tools > Superpose...
     - work with two tabs at once
   * - ``prepare_protein/``
     - Tools > Prepare Protein...
     - slow work in the background, then one undoable step
   * - ``binding_pocket``
     - Tools > Show Binding Pocket, the same and Clear Interactions in the
       right-click menu, interaction lines
     - annotations, right-click items, style edits and zooming to atoms
   * - ``ramachandran/``
     - Tools > Ramachandran Plot (:kbd:`Option Cmd R`) and its panel
     - a document panel that follows edits, frames, the selection and the
       theme
   * - ``properties/``
     - Tools > Properties (:kbd:`Cmd I`) and its panel
     - a document panel that skips work moving atoms does not affect, and
       works large structures out in the background
   * - ``molecular_surface``
     - View > Surface: show, hide, colourings and See-Through
     - a slow mesh worked out in the background and kept during drags, a
       choice, the atoms it wraps kept by id
   * - ``charges``
     - View > Charges: Color by Charge, Charge Labels, Dipole Moment and the
       charge model
     - toggles, a choice kept in the document, atom colours, labels, a 3D
       arrow and status text together
   * - ``hydrogen_bonds``
     - View > Hydrogen Bonds (:kbd:`Shift Cmd B`)
     - lines that end on the ribbon, status text, exact pictures while
       reusing slow results on screen
   * - ``mutate_residue``
     - Mutate in the right-click menu of amino acids
     - a right-click submenu
   * - ``grids/``
     - the grid panel: isosurfaces of grids, orbitals and densities
     - meshes, a document panel at the bottom, grids worked out in the
       background, state saved in sessions

Their commands keep the places they had before they became plugins: each
anchors to an item or a named separator of the core menus with ``before``
or ``after``, and the loading order decides the order of commands anchored
to the same item. Their panels are added under the names they had
(``vibrations``, ``ramachandran``, ``properties``, ``grid``), so remembered
places and folds still apply, and the settings they saved in sessions
before (``hydrogen_bonds``, ``charge_colors``, ``volume`` and so on) are
still read.

``examples/plugins/principal_axes.py`` is a smaller example of drawing: a
toggle, arrows along the axes of inertia, labels with the moments and the
kind of rotor in the status line. ``examples/plugins/template`` is the
plugin package to start from (:ref:`dg-plugin-template`).

How it works
************

The machinery is the package :mod:`miles.plugins`:
:mod:`~miles.plugins.discovery` finds plugins and records what is known
of each, :mod:`~miles.plugins.context` is what ``register`` receives,
:mod:`~miles.plugins.settings` holds toggles, choices and settings,
:mod:`~miles.plugins.layers` the drawing layers,
:mod:`~miles.plugins.background` slow results, :mod:`~miles.plugins.dialogs`
the dialog base, :mod:`~miles.plugins.colors` and :mod:`~miles.plugins.text`
the colour forms and atom labels, and :mod:`~miles.plugins.manager` loads
plugins and connects what they register to the window. :mod:`miles.plugin`
and its two submodules only gather the public names from here and
elsewhere; each public name is the name of the object where it is defined.

:class:`~miles.plugins.PluginManager` is created by the main window
after the scripting interface (``window.api``) and the menus.
:func:`~miles.plugins.discover` lists the candidates in loading order,
each with a :class:`~miles.plugins.PluginInfo` and a function that
imports it; files from the plugin folder are imported with
:mod:`importlib.util` under the name ``miles_user_plugins.<stem>``.
Entry points of the built-in plugins are skipped there, since
:data:`~miles.plugins.discovery.BUILTIN` lists them; that list also makes them
load when the package was not installed with its metadata.
:meth:`~miles.plugins.PluginManager.add_module` registers a module
directly, which is how :func:`miles.testing.load_plugin` works.

Commands are ordinary ``QAction`` objects added to the menu bar, which is
why the command palette and :meth:`miles.api.App.commands` see them
without further work. Each callback is wrapped: an exception is recorded in
the plugin's ``PluginInfo`` and handed to ``sys.excepthook``, where the
:class:`~miles.ui.errors.ErrorReporter` shows it. Panels from
``add_panel`` are :class:`~miles.ui.plugin_panel.PluginPanel` objects
named ``plugin-<key>-<title>``; they and the panels from
``add_floating_panel`` are kept by the stage in ``extra_panels``, in
loading order, and laid out after its own panels. The window restores and
saves their arrangement under ``panels/<name>``.

Drawing layers, status text and right-click items go to the window's
:class:`~miles.ui.canvas.extensions.Extensions`, which every canvas
receives; each plugin function is wrapped so that an error is reported once
and the layer is skipped for that tab until the document changes or its
``update()`` is called. Structure
sources go to the list shared by the window's and the fragment window's
:class:`~miles.ui.lookup.StructureLookup`, fragment libraries to the
fragment window, and per-tab state to the canvas's ``plugin_state``, which
sessions save under ``plugins/<key>``. The working data of
``document_cache`` lives on the tab's :class:`~miles.api.DocumentView`
(:meth:`~miles.api.DocumentView.cache`), which
:func:`~miles.api.document_view` makes once per tab and keeps with it, so
it goes when the tab goes. Toggles are checkable commands whose
check mark the manager sets when another tab comes to the front.

Events come from signals of :class:`~miles.ui.mainwindow.window.MainWindow`:
``documentOpened``, ``documentClosed``, ``currentDocumentChanged``,
``documentSelectionChanged``, ``documentStructureChanged`` (emitted for
structure, bond and coordinate changes), ``documentChanged`` (every change,
with its kind; the manager also derives the frame event from it) and
``themeChanged``. The manager hands each plugin the tab's ``DocumentView``
and calls the plugins' functions in turn.
