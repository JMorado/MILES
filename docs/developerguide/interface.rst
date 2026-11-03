The Interface
#############

The ``ui`` package is the Qt application. This chapter describes how its
parts fit together: documents and undo, the shared editing state, the canvas
and its tools, the menus, and the floating panels.

Documents and undo
******************

A :class:`~miles.ui.document.Document` owns one molecule, the display
state of its atoms, its file path and title, and its undo history. Its
logic is in :class:`miles.core.document.Document`, which needs no Qt (see
:doc:`molecule_model`); the window's class adds the signals ``changed``
(the kind of change), ``changed_detail`` (the whole
:class:`~miles.core.document.ChangeEvent`) and ``renamed``. Every edit
goes through :meth:`~miles.ui.document.Document.edit`:

.. code-block:: python

   doc.edit(building.add_hydrogens, selection)          # undoable, notifies views
   doc.edit(clean, kind=Change.COORDS)                   # only coordinates changed

``edit(func, *args, kind=...)`` puts a snapshot of the molecule and its
display on the undo stack, calls ``func(mol, *args)``, increments the
document's revisions and emits ``changed`` with a
:class:`~miles.ui.document.Change` saying what kind of change it was. Undo
swaps the current molecule and display with the snapshot on the stack.
Representations and focus-mode hiding change through
``set_representation``, ``set_hidden`` and ``edit_display``, which take
the same path with the kind ``style``.

Because of this, the functions in ``core`` know nothing about undo or about
the interface. They take a molecule and change it, and any of them can be
made an undoable command by passing it to ``edit``.

Edits made in many small steps, such as dragging atoms or turning the dial,
call :meth:`~miles.ui.document.Document.checkpoint` once at the start,
change the molecule in place and call
:meth:`~miles.ui.document.Document.notify` after each step, so the whole
drag is one undo step.

The ``revision`` number is what caches compare: whenever it changes, the
canvas knows the structure may have changed. ``layout_revision`` changes
with everything except atoms moving, so what does not depend on positions
(representations, colours, cartoon residues, the hit tester) is kept while
a trajectory plays or atoms are dragged.

Shared editing state
********************

:class:`~miles.ui.builder.Builder` holds the state that every window
and tab shares: the active :class:`~miles.ui.builder.Tool` (select,
build or delete), the current fragment with its preview rotation, and the
:class:`~miles.ui.builder.MouseMode` (rotate, translate or scale). It
emits ``toolChanged``, ``fragmentChanged`` and ``mouseModeChanged``, which
the rail, the menus, the panels and the canvases listen to.

The builder also interprets typed letters. A letter that is an element symbol
chooses that element; a second letter typed within 0.9 seconds completes a
two-letter symbol. The letters ``r``, ``t`` and ``s`` switch the mouse mode,
as in VMD, except that ``s`` stays sulfur with the build tool. If such a
letter turns out to start a two-letter symbol, the previous mouse mode is
restored. Clicking the canvas ends symbol typing.

The canvas
**********

The canvas is the view of one document: what is shown of it, the selection
(kept by the document), input and the status line. It is one widget,
:class:`~miles.ui.canvas.canvas.Canvas`, that owns parts with one job
each; the parts reach the document and the canvas's state through the
canvas, and the canvas calls them, so no part shares private state with
another. It is a package, ``ui/canvas``, with one module per part.

``canvas.py``
   :class:`~miles.ui.canvas.canvas.Canvas`: the document, the display
   options (``canvas.display``), the theme,
   the selection, picking, toasts and the status line, and the Qt events,
   which it paints itself or hands to ``canvas.input``.
   It declares the canvas signals once: ``selectionChanged``, ``focused``,
   ``nudge``, ``playRequested``, and ``documentChanged`` and
   ``documentRenamed``, which pass on the document's own signals so that
   the window's connections go with the canvas when its tab closes.

``view.py``
   :class:`~miles.ui.canvas.view.ViewController` (``canvas.view``): the
   camera and the smooth moves between views (fit, centre, reset, look
   along an axis, zoom to atoms), with the animation's timer.

``scene_cache.py``
   :class:`~miles.ui.canvas.scene_cache.SceneCache` (``canvas.scenes``):
   the cached scene arrays and picker.

``drawing.py``
   :class:`~miles.ui.canvas.drawing.Drawing` (``canvas.drawing``): the
   frame to draw, highlights, what plugins draw (with its cache) and the
   overlay of exported pictures.

``input.py``
   :class:`~miles.ui.canvas.input.InputController` (``canvas.input``):
   mouse, trackpad and keyboard, turned into view moves, tool actions and
   commands (see below).

``commands.py``
   :class:`~miles.ui.canvas.commands.Commands` (``canvas.commands``): what
   the menus do to the document and the view (undo, the clipboard, select a
   group, representations, focus mode, measurements, clean, add hydrogens).
   On structures of :data:`~miles.ui.canvas.commands.BACKGROUND_ATOMS`
   atoms or more, Clean Geometry and Add Hydrogens work on a copy on a
   worker thread and apply the result as one undoable step, unless the
   document changed meanwhile.

``surfaces.py``
   What draws the molecule, chosen once per window by
   :func:`~miles.ui.canvas.surfaces.make_surface`.
   :class:`~miles.ui.canvas.surfaces.PainterSurface` paints it with
   QPainter inside each canvas. :class:`~miles.ui.canvas.surfaces.GLView`
   is one ``QOpenGLWidget`` that the stage keeps under its canvases: it
   draws the scene of the canvas shown, so all tabs share one OpenGL
   context, one set of compiled shaders and one framebuffer; the canvases
   above it are transparent and paint only their overlay. Switching tabs
   sends the new canvas's arrays to the GPU once.

``view_state.py``
   What sessions keep of a canvas: :func:`~miles.ui.canvas.view_state.view_state`
   and :func:`~miles.ui.canvas.view_state.apply_view_state` (camera,
   display options, charge model and plugin state).

``tools.py``
   One class per tool (see below).

``overlay.py``
   Atom labels, the numbers of selected atoms, the rubber band while
   drawing a bond, the status line, hover text and toasts. Only atoms
   inside the window get a label, at most
   :data:`~miles.ui.canvas.overlay.LABEL_LIMIT` (the ones nearest the
   viewer); while the view turns, more than
   :data:`~miles.ui.canvas.overlay.MOVING_LABEL_LIMIT` labels are left out.
   The outline of each text is kept between paints; above
   :data:`~miles.ui.canvas.overlay.OUTLINED_LABELS` labels, cached
   pictures of the texts are copied instead of filling outlines, which
   keeps a repaint with 1,000 labels at a few milliseconds.

``context_menu.py``
   The right-click menu.

The API (:class:`~miles.api.DocumentView`) and plugins use the canvas's own
methods and attributes (``doc``, ``display``, ``selection``,
``set_style``, ``fit``, ``render_export`` and so on), not its parts.

Input
=====

A mouse press creates a :class:`~miles.ui.canvas.tools.Gesture`
recording where it started, which button was pressed and what dragging will
do, a :class:`~miles.ui.canvas.tools.Drag`. The choice is made in this
order:

1. With Option, drag the selection (or the molecule under the cursor):
   Shift moves it, otherwise it rotates.
2. With Shift, pan. With Cmd, spin about the line of sight.
3. Otherwise the active tool may claim the press: the build tool starts
   drawing a bond from an atom, the select tool starts moving the selection
   when a selected atom is pressed.
4. Otherwise the mouse mode decides: rotate, translate or scale.

The middle button always pans and the right button zooms. Moves smaller than
three pixels are not drags: releasing then counts as a click, which goes to
the tool (or opens the context menu for the right button or a
Control-click).

Tools
=====

Each tool is a small stateless object, a subclass of
:class:`~miles.ui.canvas.tools.CanvasTool`, kept in the ``TOOLS``
dictionary. The canvas asks the active tool four things:

``highlight(canvas, hit)``
   How to highlight the atom under the cursor. The build tool shows
   hydrogens in green, because clicking one attaches the fragment there;
   the delete tool shows everything in red.

``press(canvas, gesture, hit)``
   Whether a plain left press starts a tool-specific drag.

``click(canvas, hit, pos, additive)``
   What a click does: select, build or delete.

``double_click(canvas, hit)``
   What a double click does. The select tool selects the residue, then the
   chain, then everything; the build tool treats a quick second click as
   another build click.

This is the strategy pattern: the canvas contains no tool-specific branches,
and adding a tool means writing one class (see :doc:`extending`).

The scene cache
===============

Building a scene walks every atom and bond, so
:class:`~miles.ui.canvas.scene_cache.SceneCache` rebuilds it only when
the document revision, the display style or the theme changes. While only
positions change (a trajectory frame, a drag step), it keeps the scene
layout and fills in the new positions. Highlights and the cartoon tint
change on every mouse move: the highlights are an array of their own, one
value per sphere, worked out from the selection and the atom under the
cursor with array operations, and the cartoon tint is redone only when the
selection changes. Each frame the cache returns a
:class:`~miles.ui.canvas.scene_cache.Frame` with the arrays and flags
saying which of them changed, and the OpenGL canvas uploads only those to the
GPU.

The cache also holds the picker, rebuilt under the same rule, and the
bounding sphere of the structure, which gives the depth range for fading
and for the projection.

The status line
===============

The text at the bottom left shows the formula, the molar mass and, for
structures of up to 500 atoms, the point group. It is recomputed 120
milliseconds after the last change rather than on every change, so dragging
atoms does not run the point group search on every mouse move.

Menus
*****

The menu bar is declared as data in ``ui/menus.py``. Each item names what it
calls as a dotted path, ``slot``, which :func:`~miles.ui.menus.resolve`
looks up: a method of the main window (``"new_document"``) or of one of its
parts (``"files.save"``, ``"dialogs.show_settings"``,
``"builder.set_tool"``), or of the canvas in front, looked up each time the
item is chosen (``"canvas.commands.paste"``, ``"canvas.view.fit"``);
``"structure."`` is the canvas in front only when it has atoms, otherwise
the item does nothing. There are no forwarding methods in between.
``MENUS`` is a list of :class:`~miles.ui.menus.Menu` objects whose items are:

``Command(text, slot, shortcut, args)``
   calls ``slot(*args)``;

``ToggleItem(text, slot, name, shortcut, checked, option)``
   a checkable item that calls ``slot(checked)``, or
   ``slot(**{option: checked})`` when ``option`` is set (the display
   options, through ``"canvas.set_style"``);

``ChoiceItem(title, name, slot, options, option)``
   mutually exclusive items, in a submenu or in place, that call
   ``slot(value)`` (or ``slot(**{option: value})``);

``Submenu(title, items)``
   a nested menu;

``Dynamic(title, slot, name)``
   a submenu whose items change, such as **File > Open Recent**:
   ``slot(menu)`` clears and fills the ``QMenu`` when the menu is
   built and again each time it opens;

``SEPARATOR``.

:func:`~miles.ui.menus.build_menus` creates the Qt actions and returns
the ones that are updated later, such as the undo action and the checked item
of each choice. To keep a choice in step with the program, call
``menu_actions.check(name, value)``; dynamic submenus are kept in
``menu_actions.menus`` under their name, so they can be refilled when their
content changes (the recent files are refilled after every open and save,
which keeps the command palette up to date too).

The main window
***************

:class:`~miles.ui.mainwindow.window.MainWindow` builds its parts and
connects them; each is an object of its own:

* :class:`~miles.ui.mainwindow.tabs.DocumentTabs` (``window.tabs``): the
  open documents. The tab bar is the one record of them: each tab keeps its
  canvas in its data, so ``window.canvases`` is read from the bar, tabs can
  be dragged into another order, and the current document is announced
  (``currentDocumentChanged``) once per change.
* :class:`~miles.ui.mainwindow.files.FileController` (``window.files``):
  opening, saving, recent files, exported pictures and data, and fragment
  libraries.
* :class:`~miles.ui.mainwindow.dialogs.WindowDialogs` (``window.dialogs``):
  the popovers and dialogs, and the offer to restore documents.
* :class:`~miles.ui.stage.Stage` (``window.stage``): the canvases, the
  OpenGL view under them and the panels over them.

The window's own methods are what scripts, plugins and the command line use
(``current_canvas``, ``canvases``, ``new_document``, ``open_path``,
``open_structure``, ``save_document``, ``set_appearance``, ``run_script``,
``run_code`` and a few more); the menus call the parts directly.

The recent files live in the settings under ``recent_files``, a list of
absolute paths, most recent first, at most
:data:`~miles.ui.mainwindow.files.RECENT_LIMIT` long.

Plugins add their commands after ``build_menus``, as plain ``QAction``
objects inserted into these menus (or into menus they create before
Help), so nothing else needs to know about them; see :ref:`dg-plugins-window`.

Dialogs
*******

Dialogs are drawn in the app's theme, which may differ from the system's
appearance (a dark app on a light system, or the other way round). Each
dialog sets :func:`~miles.ui.controls.dialog_stylesheet`, which gives
the background, labels, plain buttons, check boxes, number and text fields,
tables and text browsers the theme's colours; text fields that should look
like fields get the object name ``Field``. Combo boxes keep their native
shape and only take the theme's colours through
:func:`~miles.ui.controls.theme_palette`, because a stylesheet would
replace the native popup. The atom list stays open while the appearance
changes, so the main window restyles open atom lists in ``apply_theme``.
Combo boxes inside floating panels also get their text colour from the
panel stylesheet (``#Panel QComboBox``): a panel that is hidden when the
theme changes is restyled when it is shown, and on macOS that restyling
brought back the palette of the old theme.

Unexpected errors
*****************

:func:`miles.ui.errors.install` makes an
:class:`~miles.ui.errors.ErrorReporter` the exception hook before the
main window is built. PySide passes exceptions raised in slots, timers and
event handlers (paint events included) to ``sys.excepthook``; the reporter
prints the traceback, as before, and shows it in a dialog with the details
folded away, so the user learns that something failed and can save. Only
one dialog is open at a time, and an error that repeats (one raised by
every repaint, for example) is shown once.

Starting the window
*******************

The ``miles`` command (:func:`miles.ui.app.main`) and
:func:`miles.launch` start the window the same way:
:func:`~miles.ui.app.application` asks for the OpenGL format and
creates the ``QApplication``, or takes the one already running (IPython's
Qt support may have made it), and gives it MILES's names, which Qt
uses for the plugin folder and similar paths;
:func:`~miles.ui.app.build_window` makes the main window with its
plugins and shows it. What differs is who runs Qt's event loop.

The command runs it itself with ``app.exec()``, after a ``SIGINT`` handler
and a timer that wakes Python every 200 ms so that :kbd:`Ctrl C` in the
terminal quits.

:func:`~miles.interactive.launch` never blocks. In IPython it calls
``enable_gui("qt")`` on the shell, so the terminal runs the loop while it
waits for a line and the Jupyter kernel while it waits for a request; if
the shell already runs another toolkit's loop, ``launch`` refuses before
anything is built, since only one can run. The application no longer quits
when its last window closes, because that would end the event loop the
kernel keeps entering; a closed window is deleted at the next ``launch``,
which makes a new one. :func:`~miles.interactive.running_app` tells
whether a window is open, so ``launch`` returns the open window's ``app``
instead, also when it was called from the console of a window started on
the command line. In plain Python nobody runs the loop between lines, so
``launch`` only warns at the interactive prompt, and scripts end with
:meth:`~miles.api.App.wait_until_closed`, which runs a nested
``QEventLoop`` until the window is hidden and turns :kbd:`Ctrl C` into a
``KeyboardInterrupt``.

IPython puts its own ``sys.excepthook`` in place while a cell runs and
restores the previous hook afterwards, so a hook set inside the cell that
called ``launch`` would be undone. ``launch`` therefore installs the error
reporter from IPython's ``post_run_cell`` event. From then on, errors
raised inside the window while IPython waits open the error dialog, and
the reporter prints their traceback to ``sys.stderr``, which the kernel
sends to the output of the last cell. Errors raised while a cell runs reach
IPython's own hook and show in that cell.

ipykernel runs cells on an asyncio loop and enters Qt's loop from it,
leaving Qt's loop only when a socket notifier reports a new request. Now
and then the asyncio side reads the next request just before the kernel
enters Qt's loop; that request waits for the previous cell's handler to
finish, which cannot happen while Qt's loop blocks asyncio, so the cell
hangs until another message arrives. This happens with a bare
``%gui qt`` as well; in a test of 150 quick cells in 30 new kernels, 6
hung without the workaround below and none with it. In a kernel,
``launch`` therefore starts a timer that
every 50 ms, while the notifier is enabled, takes a request the notifier
missed and leaves Qt's loop, as ipykernel does by itself on Windows; the
kernel enters it again a millisecond later. It relies on a few ipykernel
internals (``get_shell_stream``, ``kernel._qt_notifier``,
``kernel.app.qt_event_loop``) and turns itself off if they are missing.

Panels and popovers
*******************

The floating parts of the window are ordinary widgets placed over the
canvases by :class:`~miles.ui.stage.Stage`. The stage follows the
canvas of the current tab: :meth:`~miles.ui.stage.Stage.set_canvas`
points every panel at it, and the window forwards document and selection
changes to :meth:`~miles.ui.stage.Stage.doc_changed` and
:meth:`~miles.ui.stage.Stage.selection_changed`. Panels that appear and
disappear emit ``relayout``, and the stage moves the others out of the way:
the rail sits at the top left with the sequence strip beside it, the focus
bar at the top right, and the timeline and selection panel stack at the
bottom.

* :class:`~miles.ui.rail.Rail` shows the tools and the current fragment.
* :class:`~miles.ui.selection_panel.SelectionPanel` appears with the
  select tool. It shows atom details, or the bond length, angle or dihedral
  of 2-4 selected atoms with a :class:`~miles.ui.controls.ScrubField`
  and a :class:`~miles.ui.controls.Dial` to change it.
* :class:`~miles.ui.fragment_popover.FragmentPopover` is the fragment
  library with its search box and a rotatable preview.
* :class:`~miles.ui.symmetry_popover.SymmetryPopover` shows the point
  group and symmetrises.
* :class:`~miles.builtin_plugins.properties.panel.PropertiesPanel` lists the
  properties of the selection or the structure, stacked on the right under
  the vibrations and Ramachandran panels. The ``properties`` built-in
  plugin creates it and opens it from **Tools > Properties**
  (:ref:`dg-builtin-plugins`). It follows selection changes,
  structure edits and tab changes through a 200 ms single-shot timer, so a
  burst of changes is measured once; moves and frame changes are ignored,
  since they do not change these properties. The rows come from
  :func:`~miles.builtin_plugins.properties.panel.property_rows`, so the text that
  **Copy All** puts on the clipboard is the same as what is shown.
* :class:`~miles.ui.timeline.Timeline` appears when the document has
  frames (a trajectory or conformers) and steps or plays through them with
  :meth:`Document.set_frame <miles.ui.document.Document.set_frame>`.
  A :class:`~miles.ui.timeline.FramePlot` above the slider shows one
  :class:`~miles.ui.timeline.Quantity` per frame: the energy relative
  to the lowest (``Document.frame_energies``, set through
  :meth:`~miles.ui.document.Document.set_frames` or read from a file),
  a pinned measurement, or the RMSD from the first frame
  (:func:`~miles.ui.timeline.quantities` lists what a document offers).
  The plot is inset by half the slider handle, so its marker sits above the
  handle, and it follows the fractional playback position. Values are
  computed for all frames at once with ``core/trajectory.py`` and cached
  until the coordinates or frames change; showing another frame only moves
  the marker, so playback stays cheap. The plot shows by itself when there
  are energies or pinned measurements, and the plot button overrides that
  per document. The energies are dropped with the frames, and are only kept
  when there is one for every frame.
* :class:`~miles.ui.sequence_strip.SequenceStrip` shows the residues of
  proteins and nucleic acids; clicking and dragging selects them.
* :class:`~miles.ui.focus_bar.FocusBar` appears while focus mode hides
  atoms (the ``hidden`` field of each atom) and shows everything again.
* :class:`~miles.ui.floating.FloatingPanel` is the base of the grid,
  vibrations, Ramachandran and properties panels. Each puts its top row,
  with the base's chevron at the left, in its own layout and the rest in
  ``body``, which the chevron hides. Mouse presses that reach the panel
  itself (its background and labels; controls take their own) drag it; the
  position is kept in ``placed``, and
  :meth:`~miles.ui.stage.Stage.layout_overlays` places a moved panel
  there, inside the window, instead of in its column. A panel emits
  ``arranged`` when the user moves, folds or resets it; the main window then
  stores its place and fold under ``panels/<name>`` in the settings
  (:meth:`~miles.ui.stage.Stage.floating_panels` gives the names) and
  restores them at start-up.
* :class:`~miles.builtin_plugins.grids.panel.VolumePanel`, the panel of
  the Orbitals and Grids plugin, sits above the timeline and controls the
  isosurface of the document's volumetric data. It keeps the request of each
  document (:meth:`~miles.builtin_plugins.grids.panel.VolumePanel.request`)
  and emits ``volumeChanged``; the plugin stores the request in the
  document's plugin state, so sessions keep it, and draws the isosurface as
  a scene item, caching its mesh per grid, isovalue and colouring grid. Changes are applied after a short pause
  (40 ms), so dragging the isovalue does not recompute a surface for every
  pixel. For documents with molecular orbitals it lists the orbitals as
  well, and the electron and spin densities (menu entries such as
  ``("orbital", 7)`` or ``("density", 0)``); choosing one that has no grid
  yet evaluates it on a worker thread
  (:class:`~miles.ui.background.Background`) and appends the grid to
  ``JobInfo.volumes``, where later choices find it by name. The result is
  tagged with its document, so a grid that arrives after the user has moved
  to another tab or orbital is kept but not shown.
  :class:`~miles.builtin_plugins.grids.levels.OrbitalLevels` draws the orbital
  energies with QPainter and keeps its own energy window (zoom and pan);
  a click goes through the same path as the menu, so both stay in step.
  Which ticks share a place depends on the zoom, so the stacking is worked
  out again whenever the energy window or the width changes.
  For unsigned data the
  panel offers a second grid to colour the surface (the ``mapped`` field of
  :class:`~miles.render.surface.VolumeRequest`); **Add Grid...** asks
  the main window to read another file and append its grids, after checking
  that the atoms match, **Combine...** opens a
  :class:`~miles.builtin_plugins.grids.combine_dialog.CombineDialog` and adds the result
  (:meth:`~miles.builtin_plugins.grids.panel.VolumePanel.combine`), and
  **Save Grid...** asks it to write the chosen
  grid (:meth:`~miles.builtin_plugins.grids.panel.VolumePanel.selected_volume`) as
  a cube file; the plugin asks for the files through the scripting
  interface, which remembers the last folder.
* :class:`~miles.builtin_plugins.vibrations.panel.VibrationsPanel` shows the IR
  spectrum and the modes. Choosing a mode saves the document's frames and
  geometry, replaces the frames with the vibration and starts the timeline;
  Stop restores what was saved. The mode chosen last is remembered per
  document, and **Arrows** keeps its displacements, scaled by
  :func:`~miles.builtin_plugins.vibrations.modes.mode_arrows`, for each tab
  (``arrows_for``), which the Vibrations plugin draws as scene items.
* :class:`~miles.builtin_plugins.ramachandran.panel.RamachandranPanel` plots phi
  and psi of the current protein below the vibrations panel. The
  ``ramachandran`` built-in plugin creates it, feeds it the plugin events
  (current tab, document changes with their kind, selection and theme) and
  connects its ``exportRequested`` to the window's ``export_plot``;
  clicking a point selects that residue on the canvas.
* The plots of the vibrations panel, the timeline and the Ramachandran
  panel share a right-click menu (:mod:`miles.ui.plot_export`). Each
  of the three emits ``exportRequested`` with what was asked for (copy,
  data, spectrum curve or image) and offers ``export_name``,
  ``export_table`` and ``plot_widget``
  (:class:`~miles.ui.plot_export.PlotSource`);
  :meth:`~miles.ui.mainwindow.files.FileController.export_plot` copies the table
  or asks where to save it, so the file dialogs and the last folder are
  handled in one place. The tables come from :mod:`miles.core.tables`,
  which needs no window. Pictures are made by the plot's own ``paint``
  method, which ``paintEvent`` also calls, drawn in the light theme onto a
  white PNG at three times the size on screen or into an SVG file.

Playback is smooth because the timeline does not step through frames: it
advances a fractional position at the chosen speed about sixty times a
second, and :meth:`Document.set_position
<miles.ui.document.Document.set_position>` blends the two neighbouring
frames. Periodic frames (vibrations) also blend the last frame into the
first.

The command palette (:class:`~miles.ui.command_palette.CommandPalette`)
needs no list of its own: it walks the menu bar, so every command added to
``MENUS`` can be found by typing. A query matches when its characters appear
in order; matches at word starts, runs of consecutive characters and matches
in the command's own name (rather than only its menu path) rank higher.

Work that takes seconds runs on the application's one executor
(:func:`~miles.ui.tasks.executor`): a few worker threads that call a
function with plain arguments and deliver its result, or its error, on the
interface thread. A job is dropped when its
:class:`~miles.ui.tasks.CancelToken` was cancelled, when a newer job with
the same ``key`` replaced it (a panel recomputing after each change), or
when the Qt object that owns it was deleted. Owners are held only weakly and
the executor keeps each job until its result is delivered, so no Qt object
is ever finalised on a worker thread. A job that takes the token
(``pass_token``) checks ``token.cancelled`` between steps and stops early.
Downloads (:class:`~miles.ui.lookup.StructureLookup`) and
:class:`~miles.ui.background.Background`, which some panels use, run on it
too; only autosave keeps a writer thread of its own, because its writes must
happen in order and must all be on disk when the window closes.

The canvas keeps a few annotations of its own besides the selection:
pinned measurements and stereo labels, which are cached per document
revision.

Slow commands run on the executor through the stage's
:class:`~miles.ui.tasks.TaskRunner`; the
:class:`~miles.ui.busy_bar.BusyBar` at the top right, under the focus
bar, says what runs and cancels the task started last: its result is
dropped, and a job that takes the cancel token stops at its next check
(``App.run_in_background(..., cancellable=True)``). Help > Shortcuts
(:class:`~miles.ui.dialogs.ShortcutsDialog`) lists the gestures and
keys of :data:`~miles.ui.dialogs.SHORTCUTS` and then the shortcuts of
the menu commands, read from the menu bar by
:func:`~miles.ui.dialogs.menu_shortcuts`.

Everything else drawn over or with the molecule comes from plugins, through
the window's :class:`~miles.ui.canvas.extensions.Extensions`, which
every canvas receives when it is made. Each
:class:`~miles.ui.canvas.extensions.Layer` is a function of the canvas;
:meth:`~miles.ui.canvas.drawing.Drawing.layer_results` calls it and
keeps the result until the molecule, the document revision (which changes
with every edit and every frame shown), the display options, the theme or
the layer's ``generation`` change. Exported pictures and movie frames call
the layers again with ``export`` set, so a plugin that saves time on
screen (the Hydrogen Bonds plugin reuses a slow search for 120 ms while
frames play) is exact there.

* :meth:`~miles.ui.canvas.drawing.Drawing.annotations` gives the lines
  and labels (:mod:`miles.render.layers`).
  :meth:`~miles.ui.canvas.drawing.Drawing.resolve_lines` turns their
  anchors into points: atoms take their coordinates, lines with
  ``drawn_only`` to atoms that are not shown are dropped, and lines with
  ``ribbon`` move ends at backbone atoms that only a cartoon shows
  (:data:`~miles.ui.canvas.drawing.RIBBON_ATOMS`) to the nearest vertex of
  the cartoon mesh within 4 angstrom (found with
  :func:`~miles.core.neighbors.cross_pairs`), dropping lines to other
  undrawn atoms. The result is cached per layer result and cartoon mesh.
  :func:`~miles.ui.canvas.overlay.paint_annotations` draws the lines
  under the pinned measurements and the labels over everything, on screen
  and onto exported pictures and movie frames.
* :meth:`~miles.ui.canvas.drawing.Drawing.scene_items` merges the
  arrows, meshes and atom colours of the scene layers. The scene cache turns
  the arrow rows into a mesh only when they change, rebuilds the atoms only
  when their colours change (it compares a digest), and both renderers draw
  the meshes after the molecular surface.
* Status layers add text to the status line, worked out a moment after the
  last change or at once by
  :meth:`~miles.ui.canvas.canvas.Canvas.update_status`.
* The right-click menu (:mod:`miles.ui.canvas.context_menu`) adds the
  plugins' items after the atom's, the bond's and the selection's own.

The canvas also keeps ``plugin_state``, a dictionary per plugin that
:func:`~miles.ui.canvas.view_state.view_state` saves in session
files under ``plugins`` and that is filled from the session when the canvas
is made.

The source of partial charges is a property of the document
(``charge_model``, set with
:meth:`~miles.core.document.Document.set_charge_model`); the charges
themselves are the molecule's ``"charge"`` column, worked out when first
asked for and kept until atoms or bonds change. The canvas saves the model
in the view state.

The Python console
******************

:class:`~miles.ui.console.ConsolePanel` is the console of
**View > Python Console**. It is a
:class:`~miles.ui.floating.FloatingPanel` rather than a dock widget, so
it looks and behaves like the other panels: it folds to its title row, can
be dragged elsewhere, and its place and fold are stored with theirs under
``panels/python``. Unless it has been moved,
:meth:`~miles.ui.stage.Stage.layout_overlays` lays it across the whole
width at the bottom, and the timeline and selection panel stack above it, so
the console never covers them. A 7-pixel grip along its top edge changes its
height, which, like the history (``console/history``, the last 500 lines)
and whether it is open (``console/open``), is kept in the settings.

The interpreter is a :class:`code.InteractiveConsole` running on the
interface thread, so code can use the window directly and Qt calls are
safe; the price is that the window waits while a command runs. Each line
goes through :meth:`~miles.ui.console.ConsolePanel.execute`, which
redirects ``sys.stdout`` and ``sys.stderr`` to the output only while the
line runs. Tracebacks are written to the output by the interpreter itself
and never reach ``sys.excepthook``, so a typing mistake does not open the
error dialog. Scripts (:meth:`~miles.ui.console.ConsolePanel.run_file`,
used by **File > Run Script...** and ``--script``, and
:meth:`~miles.ui.console.ConsolePanel.run_source`, used by ``-c``) run
in the same namespace; their errors also go to the error dialog, since a
script may run without the user watching the console.

The namespace comes from :func:`miles.scripting.namespace`. It holds
``app`` (:class:`miles.api.App`), a few modules and the helpers of
:class:`~miles.scripting.Helpers`, which only call :mod:`miles.api`.
Before every command :func:`~miles.scripting.refresh` points ``doc`` and
``mol`` at the current tab; it remembers the objects it set, so a name the
user has bound to something else is left alone.

Code in the console may change a molecule without going through
``Document.edit``, for example ``mol.coords = mol.coords + 1``. The document's
``revision`` then stays the same and no view would notice. To catch this,
:func:`~miles.scripting.snapshot` records the revision and a hash of
atoms, bonds and coordinates of every document before a command, and
:func:`~miles.scripting.settle` compares them afterwards: a document
whose hash changed while its revision did not gets
``notify(Change.STRUCTURE)``, and every canvas is redrawn. Such a change is
drawn but is not a separate undo step.

.. _dg-plugins-window:

Plugins in the window
*********************

The main window creates a :class:`~miles.plugins.PluginManager` at the
end of its constructor, after ``api`` and the console, and asks it to load
the plugins (third-party ones only when ``--no-plugins`` was not given and
``MILES_NO_PLUGINS`` is not set). :doc:`plugins` describes the
interface plugins see; here is what the window provides for it.

Signals report events to plugins, the document ones carrying the canvas
concerned: ``documentOpened`` at the end of ``_add_document``,
``currentDocumentChanged`` from ``_tab_changed``,
``documentSelectionChanged`` from ``_selection_changed``,
``documentStructureChanged`` from ``_doc_changed`` for
``Change.STRUCTURE``, ``BONDS`` and ``COORDS``, ``documentChanged`` from
``_doc_changed`` for every change, with its kind, and ``themeChanged`` at
the end of ``apply_theme``, with the new theme.

A plugin's panel is a :class:`~miles.ui.plugin_panel.PluginPanel`, or
a :class:`~miles.ui.floating.FloatingPanel` of its own.
:meth:`~miles.ui.stage.Stage.add_panel` keeps it in
``stage.extra_panels``, gives it the theme's icon colour and lays it out
below the vibrations panel in the right-hand column, in the order plugins
added them; :meth:`~miles.ui.mainwindow.window.MainWindow.track_panel`
restores its place and fold from the settings and saves them when the user
moves or folds it, as ``_restore_panels`` and ``_save_panel`` do for the
core's panels.

The window knows nothing of the features that are plugins: their commands
are not in :mod:`miles.ui.menus`, the window has no slots for them and
the stage does not create their panels. Their dialogs and panels
(:mod:`miles.builtin_plugins.scan.dialog`, :mod:`miles.builtin_plugins.superpose.dialog`,
:mod:`miles.builtin_plugins.sequence.dialog`, :mod:`miles.builtin_plugins.prepare_protein.dialog`,
:mod:`miles.builtin_plugins.ramachandran.panel` and
:mod:`miles.builtin_plugins.properties.panel`) stay in ``ui`` as reusable widgets
that only the plugins import.

**Help > Plugins...** opens :class:`~miles.ui.plugins_dialog.PluginsDialog`:
a table of :class:`~miles.plugins.PluginInfo` rows with a check box
that writes ``plugins/disabled`` in the settings, the details of the
selected plugin (from :func:`~miles.ui.plugins_dialog.plugin_details`)
and a button that opens the plugin folder. The first failed plugin is
selected when the dialog opens, so its error is shown straight away. Its
**Settings...** button emits ``settingsRequested`` with the plugin's key;
the window then opens Settings at that plugin.

Settings
********

**Settings...** (:meth:`~miles.ui.mainwindow.dialogs.WindowDialogs.show_settings`)
opens :class:`~miles.ui.settings_dialog.SettingsDialog`, a single
non-modal window. Its menu item has the ``PreferencesRole``, so macOS moves
it to the application menu. The window is one form, built again each time
it opens (plugins may have been switched on or off), and every field writes
its setting as soon as it changes. The keys are:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Key
     - Value
   * - ``appearance``
     - ``system``, ``light`` or ``dark``
   * - ``autosave/interval``
     - minutes between autosaves; ``0`` switches autosave off
   * - ``defaults/mode``, ``defaults/labels``, ``defaults/cartoon_color``
     - style, labels and cartoon colours of new documents
   * - ``defaults/show_hydrogens``, ``defaults/depth_fade``
     - ``true`` or ``false`` for new documents
   * - ``plugins/settings/<plugin>/<name>``
     - a plugin's setting (:class:`~miles.plugins.Setting`)

:func:`~miles.ui.settings_dialog.new_document_style` reads the
``defaults`` group, and the window applies it to each new canvas before a
session's own view is restored. All fields, drop-down lists included, are
drawn by :func:`~miles.ui.settings_dialog.fields_stylesheet`, so they
have one width and line up in a column on every platform.

Autosave
********

:class:`~miles.ui.autosave.Autosave` belongs to the window. Each new
document is tracked; an announced change other than showing another frame
marks it for the next copy, and saving it (``renamed`` with the document
no longer modified) or closing its tab deletes its copy. Every
``autosave/interval`` minutes, :meth:`~miles.ui.autosave.Autosave.save_now`
takes a :func:`~miles.ui.autosave.snapshot` of each marked document
(the session data, with its frames and the canvas's view) on the interface
thread and hands it to a single writer thread, which encodes, compresses
and writes it, through a temporary file and a rename. For 4HHB with
hydrogens (10,000 atoms) the window is held up for about 4 ms and the copy
is on disk 16 ms later; with 50 frames, 11 ms and 180 ms (a 6 MB file).
Writing the same copy on the interface thread took half a second. Removals go
through the same thread, so a copy is never written after the request to
delete it.

Each window writes into its own folder, ``<pid>-<random>`` under
:func:`~miles.ui.autosave.autosave_folder` (``$MILES_AUTOSAVE_DIR``
or ``autosave`` in the application data folder), and holds a
:class:`~PySide6.QtCore.QLockFile` in it while it is open, created at the
first copy. A ``.miles`` file and a ``.json`` file with its tab name, file
and time make one copy. Closing the window normally deletes the folder.

When the window starts,
:meth:`~miles.ui.mainwindow.dialogs.WindowDialogs.offer_recovery` asks
:meth:`~miles.ui.autosave.Autosave.find_leftovers` for the copies of
windows that stopped: a folder whose lock can be taken belongs to a
process that is gone, since the lock is taken with no stale time and
``QLockFile`` only gives up the lock of a process that no longer runs.
Folders of the running process itself are skipped: a window deleted
without being closed releases its lock, but its process did not crash. The
folders found stay locked while the
:class:`~miles.ui.recovery_dialog.RecoveryDialog` is open, so a second
window starting at the same time does not offer them too. Restore opens
each copy with its original file as a modified tab and deletes the copies
(a copy that cannot be read is kept for another try); Discard deletes them;
closing the dialog unlocks them for the next start. The test suite points
``MILES_AUTOSAVE_DIR`` at a new folder for every run.

Sessions
********

Saving to a session file stores, besides the molecule, the document's
frames and the canvas's view:
:func:`~miles.ui.canvas.view_state.view_state` returns the camera
and display options as JSON-ready data, and
:func:`~miles.ui.canvas.view_state.apply_view_state` restores them
when the session is opened. Per-atom representations, hidden atoms and
pinned measurements travel with the molecule itself.

They are styled with one Qt stylesheet built from the theme by
:func:`~miles.ui.controls.stylesheet`, and widgets pick their look
through their object name (``Panel``, ``Popover``, ``Seg``, ``Tile`` and so
on).

Searches by name and PDB downloads can take seconds, so
:class:`~miles.ui.lookup.StructureLookup` runs them on a worker thread
and reports back with Qt signals, which Qt delivers on the interface thread.

Text
****

``ui/text.py`` collects the strings shared between widgets: units, atom
labels such as ``C12``, residue labels and the one-line descriptions of
atoms and bonds shown when hovering.
