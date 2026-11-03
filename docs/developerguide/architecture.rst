Architecture
############

This chapter describes how MILES is put together: the layers of the
package, the objects that exist while the program runs and how they relate,
and what happens when the user edits a molecule or the window is redrawn.
Later chapters describe each part in detail.

Layers
******

The package is divided into four layers, shown in Figure 2.1. Each layer may
use the layers below it and never the ones above.

.. raw:: html
   :file: ../_static/diagrams/layers.svg

.. rst-class:: diagram-caption

   Figure 2.1. The layers of the package. Arrows point from a layer to the
   layers it uses; the labels name the external libraries each layer
   imports. NumPy is used everywhere.

**core** is the molecule model and everything that computes on it: geometry,
building, perception, the Clean force field, symmetry, protein residues, the
fragment library and the bridge to Open Babel. It also holds documents: a
molecule with its display state, trajectory and undo history. It imports
neither Qt nor OpenGL (a test makes sure). You can use it from a script or
a notebook, and it is tested without a display.

**files** reads and writes structure files, parses pasted text and
downloads structures. It builds on ``core`` and Open Babel and is also free
of Qt.

**render** decides what to draw and draws it: display styles, the scene
arrays, picking, the camera, cartoon meshes, colour themes, the OpenGL
renderer and the QPainter fallback. It uses Qt only for colours and painting,
never for widgets.

**ui** is the Qt application: the main window, documents with Qt signals,
the canvas and its tools, the panels, popovers, menus and dialogs.

Three modules at the top of the package open the application to code
written outside it. :mod:`miles.api` wraps the main window in a small
scripting interface (the window, its tabs, their molecules, selections and
views) used by the Python console, scripts and plugins.
:mod:`miles.plugins` finds plugins, calls their ``register`` function
and connects what they add (menu commands, file formats, panels and event
callbacks) to the window. Two modules are the public surface, kept stable
for scripts and plugins (:ref:`dg-plugin-surface`): :mod:`miles`, the
chemistry library, which loads its names from ``core`` and ``files`` when
they are first used, and :mod:`miles.plugin`, what a plugin needs to work
with the window. They import ``ui`` only when a window exists, so the lower
layers stay free of them.

Many features are themselves plugins, in :mod:`miles.builtin_plugins`
(:ref:`dg-builtin-plugins`). They sit on top of everything else and use
:mod:`miles` and :mod:`miles.plugin` like a third-party plugin would.
Nothing below them imports them, so each can be switched off.

.. _dg-code-placement:

Where code goes
===============

Code lives in one of three places. ``core`` (with the rest of the
application) holds the app and its foundations: a module belongs there only
if code that is not a plugin uses it, so MILES started with ``--no-plugins``
still needs every core module. ``builtin_plugins/<name>`` holds one feature
each, its interface and its algorithms together; the relaxed scan, the
Ramachandran regions and the protonation table live with their plugins.
``builtin_plugins/shared`` holds what several built-in plugins use and
nothing else does (the hydrogen-bond finder, side-chain mutation, molecular
descriptors); it is not a plugin and not public. Plugins import :mod:`miles`,
:mod:`miles.plugin`, their own package and the shared one, never each other.
``tests/test_plugin_boundary.py`` checks all of this, including that every
core module has a user outside the plugins.

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * - Place
     - What it holds
     - Examples
   * - ``core`` (and ``files``, ``render``, ``ui``, ``api``, ``plugins``)
     - the application and the foundations it and scripts build on
     - the molecule model, documents and undo, selections, geometry and
       neighbours, the Open Babel bridge, building and perception, file
       formats, rendering, the window
   * - ``builtin_plugins/<name>``
     - one feature, its interface and algorithms
     - ``scan`` (dialog and relaxed scan), ``ramachandran`` (panel, phi and
       psi, regions), ``prepare_protein`` (dialog and protonation),
       ``sequence`` (dialog and peptide builder)
   * - ``builtin_plugins/shared``
     - helpers several built-in plugins use
     - ``interactions`` (Hydrogen Bonds, Binding Pocket), ``mutation``
       (Mutate Residue, Prepare Protein, Build from Sequence),
       ``properties`` (Properties, 2D drawings)

The features fall on either side like this:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Core (``core``, ``files``, ``render``, ``ui``, ``api``, ``plugins``)
     - Built-in plugins (``builtin_plugins``)
   * - the molecule model, building, editing, selection, undo, the fragment
       window and the fragment library
     - fragment searches by name and PDB code (structure sources), Open
       from PDB, Build from Sequence, Mutate Residue
   * - every file format, sessions, the format registry
     - 2D drawings, movies
   * - the rendering pipeline, display styles, partial charges and the
       surface and isosurface computations
     - the molecular surface, charge colours and labels, the dipole,
       hydrogen bonds, binding pockets, grid and orbital isosurfaces with
       the grid panel
   * - the timeline, measurements, symmetry, stereo labels, Clean, Add
       Hydrogens
     - conformers, scans, aligning frames, superposing, protein
       preparation, pattern selection
   * - the command palette, the Python console, the plugin manager,
       background tasks, Help > Shortcuts built from the menus
     - the vibrations, Ramachandran and properties panels

What a plugin draws reaches the canvas through extension points rather than
code of its own in the canvas: annotations over the molecule, items drawn
in 3D (arrows, meshes, atom colours), status text, right-click menu items,
per-tab state saved in sessions, structure sources and fragment libraries
(:ref:`dg-plugins`).

The rule that a layer never imports a layer above it keeps each part usable
on its own. The chemistry does not depend on how it is shown, the renderer
does not depend on how the user interacts with it, and the interface is a
thin layer that connects them.

Objects and their relationships
*******************************

Figure 2.2 shows the objects of a running application.

.. raw:: html
   :file: ../_static/diagrams/objects.svg

.. rst-class:: diagram-caption

   Figure 2.2. Objects of a running application. Solid arrows mean
   ownership; dashed arrows mean that one object reads or observes another.

There is one :class:`~miles.ui.mainwindow.window.MainWindow`. It owns the
floating panels, one :class:`~miles.ui.builder.Builder`, one canvas for
each open tab (in :class:`~miles.ui.mainwindow.tabs.DocumentTabs`), the
file controller and the dialogs, and, with OpenGL, one
:class:`~miles.ui.canvas.surfaces.GLView` that draws the molecule of the
tab in front.

The **Builder** holds the editing state shared by all tabs: the active tool,
the current fragment and the mouse mode. Choosing carbon in one tab chooses
it in every tab. The rail, the menus and the canvases observe it through its
signals.

Each **canvas** shows one **Document**. A document owns the current
:class:`~miles.core.molecule.Molecule`, the display state of its atoms
(:class:`~miles.core.display.AtomDisplay`: representations and focus-mode
hiding), the selection, the trajectory (:class:`~miles.core.frames.Trajectory`),
the results read with the structure (:class:`~miles.core.results.Results`)
and an undo stack of earlier snapshots. Note that, whereas a
Molecule is plain data that any function may change, a Document decides
when changes happen and tells everyone else about them. The document's
logic is in :class:`miles.core.document.Document`, which needs no Qt; the
window uses :class:`miles.ui.document.Document`, the same class with Qt
signals. Code in ``core`` that edits works on molecules and knows nothing
about documents; code in ``ui`` changes molecules only through their
documents. :doc:`molecule_model` describes both.

Besides its document, a canvas owns everything about how the molecule looks
in that tab, in parts with one job each: the camera
(:class:`~miles.ui.canvas.view.ViewController`), the display style, a
:class:`~miles.ui.canvas.scene_cache.SceneCache` with the arrays ready to
draw and the picker, what plugins draw (:class:`~miles.ui.canvas.drawing.Drawing`),
input (:class:`~miles.ui.canvas.input.InputController`) and the menu
commands (:class:`~miles.ui.canvas.commands.Commands`). What a click does
is delegated to the active **tool**, one of three stateless objects shared
by all canvases.

The **panels** (the tool rail, the selection panel and the popovers) belong
to the window, not to a canvas. They observe the builder and follow the
canvas of the current tab: switching tabs points the selection panel at the
new canvas.

The life of an edit
*******************

Every change to a molecule follows the same path. Figure 2.3 shows it for
the most common edit, clicking a hydrogen with the build tool to attach a
fragment.

.. raw:: html
   :file: ../_static/diagrams/edit.svg

.. rst-class:: diagram-caption

   Figure 2.3. Calls made when a click attaches a fragment.

1. The canvas receives the mouse release and asks the picker what is under
   the cursor. The answer is an :class:`~miles.render.picking.AtomHit`
   for the hydrogen.
2. The canvas passes the hit to the active tool's ``click`` method. The
   canvas itself has no tool-specific code.
3. The build tool calls :meth:`Document.edit
   <miles.ui.document.Document.edit>` with the ``core`` function to run,
   :func:`~miles.core.building.attach_fragment`, and its arguments.
4. The document puts a snapshot of itself on its undo stack: the
   molecule, display, selection, trajectory and results. The snapshot
   shares the molecule's columns, so this costs nothing until the edit
   changes one of them.
5. The document calls the function, which changes the molecule in place.
   If the function raises, the document goes back to the snapshot and
   forgets the step before the error goes on.
6. The document carries the display over to the changed atoms and
   increments its revision numbers. Caches compare revisions to decide
   whether to rebuild.
7. The document announces a :class:`~miles.core.document.ChangeEvent`;
   in the window it is the ``changed`` signal with the kind of change. The
   canvas schedules a repaint; the selection panel refreshes its values;
   the window marks the document as modified and enables Undo.

Every editing command, from the context menu to the atom list dialog, goes
through step 3, so every edit can be undone and every view stays up to date
without knowing who made the change.

Edits made in many small steps, such as dragging atoms or turning the dial
of the selection panel, record one undo step at the start
(:meth:`~miles.ui.document.Document.checkpoint`), change the molecule
in place as the mouse moves and announce each step with
:meth:`~miles.ui.document.Document.notify`. The whole drag is undone at
once.

Drawing a frame
***************

Drawing has to be fast at two very different rates. The structure changes a
few times a second at most, but the camera and the highlighted atom change
on every mouse move. The pipeline in Figure 2.4 separates the two.

.. raw:: html
   :file: ../_static/diagrams/frame.svg

.. rst-class:: diagram-caption

   Figure 2.4. From molecule to pixels. The scene cache does the expensive
   work only when its inputs change.

When Qt asks a canvas to paint:

1. The canvas asks its scene cache for a frame, passing the molecule, the
   document revision, the display style, the theme, the highlights and the
   selection.
2. If the revision, the style or the theme changed since the last frame, the
   cache calls :func:`~miles.render.scene.build_scene`, which walks the
   atoms and bonds once and produces the sphere and cylinder arrays and, if
   needed, the cartoon mesh. It also recomputes the bounding sphere used for
   depth fading.
3. If only the highlights changed, the cache copies the sphere array and
   sets the highlight column of the affected atoms. If the selection
   changed, it recolours the cartoon vertices of the selected residues.
4. The cache returns a :class:`~miles.ui.canvas.scene_cache.Frame` with
   the arrays and a flag for each array saying whether it changed.
5. With OpenGL, the window's one :class:`~miles.ui.canvas.surfaces.GLView`
   under the canvases uploads the changed arrays to the GPU and draws them
   with the current camera matrices. Rotating the view uploads nothing.
   Without OpenGL, the canvas draws the frame back to front with QPainter.
6. The overlay (labels, selection numbers, status line and messages) is
   painted on top by the canvas itself, with QPainter.

Opening a file
**************

Reading a structure ties the layers together:

1. The window calls :func:`~miles.files.formats.read_file` with the
   path.
2. The format registry looks up the extension. Session files have their own
   reader, and plugins may register readers for other extensions; every
   other extension goes to Open Babel.
3. Open Babel reads the file. :func:`~miles.core.chem.from_openbabel`
   converts the result into a molecule, keeping residue data for formats
   that have it, and infers each atom's bonding pattern from its bonds.
   Further structures in the file become trajectory frames.
4. If the structure contains a protein or nucleic acid,
   :func:`~miles.core.protein.default_representations` works out the
   representations it opens with (the polymer as a cartoon, water hidden)
   and estimates the secondary structure when the file has none.
5. The window wraps the molecule in a new document, gives the document's
   display those representations (unless a session file brought its own),
   creates a canvas for it in a new tab, and fits the view.

Where to extend
***************

The places where the program is meant to grow are explicit:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - To add
     - Change
   * - a file format
     - register a reader and writer in the format registry
   * - a ring or functional group
     - add an entry to the ring table or the group recipes
   * - a tool
     - write a ``CanvasTool`` subclass and add it to ``TOOLS``
   * - a menu command
     - add a ``MainWindow`` method and an entry in ``MENUS``
   * - a representation
     - add a ``Representation`` member and its radius and name
   * - a feature kept apart from the main code
     - write a plugin: a module with ``register(context)``

:doc:`extending` shows each of these with code, and :doc:`plugins`
describes plugins.
