.. _ug-python:

Python Console and Scripts
##########################

MILES can be driven from Python, by typing commands in the console at
the bottom of the window or by running a script file. Both see the same
names and act on the open tabs, as the menus do: a structure made in
Python gets a tab of its own, and anything changed is drawn at once. The
window can also be opened from a Jupyter notebook or IPython, and used
from there in the same way (:ref:`ug-jupyter`).

.. _ug-console:

The Console
***********

**View > Python Console** (:kbd:`Option Cmd P`) opens the console along the
bottom of the window, above the timeline and the status line. Drag its top
edge to make it taller or shorter, click the chevron to fold it to its title
row, or drag it by the title to put it somewhere else (a double click on its
background puts it back). The console stays open when you change tabs.
MILES remembers whether it was open, its height and its place.

Type a line and press :kbd:`Return`. The value of an expression is shown
below it, and so is anything the command prints. A line that opens a block
(``for``, ``if``, ``def``, ``with`` ...) or leaves a bracket open is
continued at the ``...`` prompt; an empty line runs the block. An error
shows its traceback in red and leaves the window as it was.

.. code-block:: pycon

   >>> m, _ = miles.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
   >>> app.new(m, "aspirin")
   <aspirin: 21 atoms>
   >>> select('smarts "[CX3](=O)[OX2H1]"')
   [10, 11, 12]
   >>> show("tube")
   >>> for i in doc.selection:
   ...     print(i, mol.atoms[i].symbol)
   ...
   10 C
   11 O
   12 O
   >>> print(f"C=O {dist(1, 2):.3f} A, C-O {dist(1, 3):.3f} A")
   C=O 1.224 A, C-O 1.391 A
   >>> round(angle(2, 1, 3), 1)
   126.6

Atoms are numbered from 0 in Python, one less than the numbers shown on the
canvas and in the atom list.

:kbd:`Up` and :kbd:`Down` step through earlier commands; the last 500 are
kept between sessions. :kbd:`Tab` completes the name before the cursor, or
lists the choices when there are several. :kbd:`Control L` (Ctrl+L on Linux
and Windows) or the **Clear** link empties the output, which can be selected
and copied like any text. Pasting several lines runs them one after the
other, as if they had been typed.

``help()`` lists what the console offers, and ``help(app)``, ``help(doc)``
or ``help(select)`` describe each part. ``exit()`` does not close MILES;
close the window to quit.

What the console knows
======================

These names are ready in the console and in scripts:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - What it is
   * - ``app``
     - The window (:class:`miles.plugin.App <miles.api.App>`): ``app.current``,
       ``app.documents``, ``app.open(path)``, ``app.new(molecule, title)``,
       ``app.run("Clean Geometry")``, ``app.commands()`` (every menu
       command), ``app.message(text)``.
   * - ``doc``
     - The current tab (:class:`miles.plugin.DocumentView <miles.api.DocumentView>`): its
       ``molecule``, ``selection``, ``select()``, ``style()``, ``edit()``,
       ``editing()``, ``undo()``, ``frames``, ``frame``, ``fit()``,
       ``turn()``, ``export_image()`` and ``save()``.
   * - ``mol``
     - The molecule of the current tab
       (:class:`miles.Molecule <miles.core.molecule.Molecule>`): ``mol.atoms``,
       ``mol.coords`` (an N x 3 NumPy array in angstrom), ``mol.bonds``.
   * - ``np``, ``miles``, ``Molecule``
     - NumPy, and the chemistry library (:doc:`../api/library`): SMILES,
       minimisation and conformers (``miles.from_smiles``,
       ``miles.minimize``, ``miles.conformers``), adding hydrogens
       (``miles.add_hydrogens``), measurements (``miles.distance``,
       ``miles.angle``), selections (``miles.find``) and file reading and
       writing.

``doc`` and ``mol`` follow the current tab: after you click another tab, or
after a command opens a new one, they name that tab. If you give either name
a value of your own, it keeps it until you delete it with ``del doc``.
Everything else in ``miles`` can be imported as usual.

Short commands
==============

A few helpers save typing for what is done most often. They act on the
current tab and are written with ``app`` and ``doc``, so ``help(select)``
and the :mod:`miles.scripting` reference show how to do the same in a
script of your own.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Command
     - Effect
   * - ``load("file.pdb")``
     - Open a file in a new tab
   * - ``select("protein and name CA")``, ``select([0, 1, 2])``, ``select()``
     - Select the atoms an expression picks (see :ref:`ug-selections`), or
       these atoms, or nothing; ``add=True`` adds to the selection. Returns
       the selected atoms.
   * - ``count("water within 3.5 of ligand")``
     - How many atoms an expression picks, without selecting them
   * - ``define("pocket", "byres within 5 of ligand")``
     - Name an expression; later ones can use ``pocket`` as a word
   * - ``show("tube")``
     - Style: ``ball``, ``tube``, ``wire``, ``spacefill`` or ``cartoon``
   * - ``fit()``, ``turn(y=90)``
     - Fit the structure to the window; turn the view, in degrees
   * - ``png("figure.png", scale=3)``
     - Save a picture of the view (``background`` can be ``"white"``,
       ``"transparent"`` or ``"theme"``)
   * - ``run("Add Hydrogens")``
     - Run a menu command by its name
   * - ``dist(i, j)``, ``angle(i, j, k)``, ``dihedral(i, j, k, l)``
     - Distance in angstrom, angles in degrees

Changing structures
===================

``doc.style()`` takes the display options of the **View** menu, for example
``doc.style(labels="symbol", show_hydrogens=False)``.

A change made through ``doc.edit`` or inside ``with doc.editing()`` is one
step that **Edit > Undo** (or ``doc.undo()``) takes back.
``doc.edit(func, ...)`` calls ``func(molecule, ...)``, so any function of
:mod:`miles` that changes a molecule in place can be used, as in
``doc.edit(miles.add_hydrogens)``. Other changes go in a ``with``
block:

.. code-block:: pycon

   >>> with doc.editing() as m:
   ...     m.coords = m.coords - m.coords.mean(axis=0)     # centre on the origin
   ...

Arrays read from a molecule (``mol.coords``, ``mol.column("resname")``) are
read-only: change a molecule by assigning a whole array, as above, or with
its methods, such as ``mol.set_positions(atoms, positions)`` or
``mol.atoms[0].symbol = "N"``. Reading whole columns is also the fast way
to work on large structures:

.. code-block:: pycon

   >>> water = mol.mask("resname", {"HOH"})
   >>> mol.coords[water].mean(axis=0)
   >>> mol.column("b_factor").max()          # PDB files give b_factor and occupancy

The structure can also be changed outside an edit, as in
``mol.coords = mol.coords + [0, 0, 1]``. The view still follows, but such a
change is not a step of its own, so undo cannot take it back by itself.

.. note::

   Commands typed in the console run in the window's own thread, and the
   window does not respond until they end. This is no problem for edits
   and measurements, but a long loop keeps it waiting. Slow work can run in
   the background with ``app.run_in_background(job, ..., done=...)``, as
   the conformer search and protein preparation do, or in plain Python
   without the window (see :ref:`ug-batch`), with the results opened
   afterwards.

.. _ug-scripts:

Scripts
*******

**File > Run Script...** runs a Python file in the console's namespace: the
script sees ``app``, ``doc``, the helpers and everything typed in the console
before, and what it defines is still there afterwards. ``__name__`` is
``"__main__"`` and ``__file__`` is the script's path. Its output appears in
the console, which opens if it was closed. An error shows its traceback in
the console and in the error dialog.

A script can also be given on the command line. It runs once the window is
up and the files named with it are open:

.. code-block:: console

   $ miles protein.pdb --script prepare.py
   $ miles ligand.sdf -c "show('tube'); select('element N'); png('ligand.png')"

``--script FILE`` runs a file and ``-c CODE`` runs a line of code; both can
be given, and the file runs first.

A module imported by a script, or a plugin, gets the same ``app`` with
:func:`miles.plugin.current_app <miles.api.current_app>`:

.. code-block:: python

   from miles.plugin import current_app

   app = current_app()
   for view in app.documents:
       print(view.name, len(view.molecule))

.. _ug-batch:

Without the window
==================

``import miles`` is the chemistry library of MILES, with the file formats
(:doc:`../api/library`). It needs neither the window nor Qt, so batch work
can be done in a plain Python session or a Jupyter notebook:

.. code-block:: python

   from pathlib import Path

   import miles

   for path in Path("ligands").glob("*.sdf"):
       mol, info = miles.read_file(str(path))
       coords = miles.minimize(mol)          # MMFF94, or UFF when MMFF94 has no parameters
       if coords is not None:
           mol.coords = coords
       miles.write_file(str(path.with_suffix(".xyz")), mol, info)

:mod:`miles` and :mod:`miles.plugin` are the public, stable part of MILES;
every other module is internal and may change between releases
(:ref:`dg-stability`).

The results open in MILES like any other file.

.. _ug-jupyter:

From Jupyter and IPython
************************

The window can also be opened from Python that is already running: a
Jupyter notebook (JupyterLab, the classic notebook or an editor's
notebooks), qtconsole, or IPython in a terminal. :func:`miles.launch`
opens the window and returns the same ``app`` the console has:

.. code-block:: python

   import miles

   app = miles.launch("ligand.sdf")       # files to open are optional
   doc = app.current
   doc.select("element N O")
   doc.style(mode="tube", labels="symbol")
   doc.image()                                  # a picture of the view below the cell

The window keeps working between cells: you can turn the structure, select
and build with the mouse, then go on in the notebook, and what a cell
changes shows in the window at once. ``launch`` turns on IPython's Qt
support for this, as ``%gui qt`` does. The kernel has to run on the
computer you are working at, since that is where the window opens.

Jupyter is not part of MILES's environment. Install it there, for
example with ``micromamba install -n miles jupyterlab`` (or only
``ipykernel`` to use the environment as a kernel from another Jupyter).

``launch`` takes a few options:

* ``appearance="dark"`` (or ``"light"``, ``"system"``) chooses the theme for
  this window only; **View > Appearance** is not changed.
* ``plugins=False`` starts without third-party plugins, like
  ``miles --no-plugins``.

Calling ``launch`` again while the window is open returns the same
``app``, opens the files given and brings the window to the front. After
the window has been closed, the next call opens a new one.

Pictures
========

``doc.image()`` returns a picture of the view, trimmed to the structure,
which the notebook shows below the cell. It is sharp on high-resolution
screens. ``scale=2`` shows it twice as large, ``background`` can be
``"white"`` (the default), ``"transparent"`` or ``"theme"``, ``crop=False``
keeps the whole view, and ``doc.image().save("view.png")`` writes it to a
file.

What works differently
======================

* The short commands of the console (``select``, ``show``, ``png`` and the
  others) are not defined in a notebook. Use the methods of ``doc`` and
  ``app``, or make them with
  ``from miles.scripting import Helpers; h = Helpers(app)`` and call
  ``h.select("element N")``, ``h.show("tube")`` and so on.
* An error in a command that a cell runs shows in that cell, as any error
  does. An error inside the window while the notebook waits, in a menu
  command or a panel, opens MILES's error dialog, and its traceback
  appears below the cell that ran last.
* Only one kind of window can run in a Python session. Matplotlib plots
  shown inline (the default in notebooks) are fine; plots in windows of
  their own need Qt too (``%matplotlib qt``). If another kind is already
  running, ``launch`` says so instead of opening the window.
* As in the console, the window waits while a cell runs, so a long
  computation keeps it from responding until the cell ends.

.. _ug-plain-python:

Plain Python
============

The plain ``python`` prompt cannot keep the window working between lines,
so ``launch`` prints a note when it is called there; use IPython instead
(``ipython``, then the same commands). A script run with ``python`` ends
with ``app.wait_until_closed()``, which hands over to the window until you
close it, and then returns:

.. code-block:: python

   import miles

   app = miles.launch("ligand.sdf")
   app.current.style(mode="tube")
   app.wait_until_closed()      # the window works until you close it
   print("done")

:kbd:`Ctrl C` in the terminal stops the waiting, with a
``KeyboardInterrupt``; the window stays open. A script can also run inside
the window with ``miles --script`` (:ref:`ug-scripts`), which needs no
``launch`` at all.

.. _ug-plugins:

Plugins
*******

A plugin is Python code that adds commands, file formats or panels to
MILES. Plugins come as single ``.py`` files or as packages installed
with ``pip``:

- **A file**: choose **Help > Plugins...**, click **Open Plugin Folder**,
  copy the file there and restart MILES.
- **A package**: install it into the environment MILES runs in, for
  example ``micromamba run -n miles pip install miles-charges``,
  and restart MILES.

A plugin's commands appear in the menus it chose, in the command palette and
in ``app.commands()``, and can be run from a script with ``app.run``.

Many of the built-in tools are plugins too, listed in Help > Plugins...
as *built in*: **File > Open from PDB...** and looking up names and PDB
codes, **File > Export Movie...**, **File > Export 2D Drawing...**,
**Edit > Copy 2D Drawing**, **Edit > Select...**, **View >
Hydrogen Bonds**, the **View > Charges** submenu, in the Tools menu **Build
from Sequence...**, **Generate Conformers**, **Scan Coordinate...**,
**Align Frames**, **Superpose...**, **Properties**, **Vibrations**,
**Prepare Protein...**, **Show Binding Pocket** and **Ramachandran Plot**,
the **View > Surface** submenu, the grid panel for grids and orbitals, and
**Mutate** in the right-click menu. They work exactly as described in this
guide, also when MILES starts with ``--no-plugins``, and like any
plugin they can be switched off in Help > Plugins..., which takes their
commands, panels and drawings out of the window at once.

Fragments of your own can be added without a plugin: put a ``.json`` file
in the format of the built-in fragment library (see :ref:`dg-plugins`) in
the ``fragments`` folder next to the plugin folder, and its categories
appear as tabs of the fragment window at the next start.

**Help > Plugins...** lists every plugin found, where it came from and
whether it loaded. Select one to see its description, the commands it added
and, if it failed, the error. A plugin that fails never stops MILES
from starting; an error in one of its commands is shown once in the error
dialog. Clear the check box in front of a plugin to switch it off, or tick
it to switch it back on; the change takes effect at once and is kept for
later starts. **Reload** loads the selected plugin again, reading its file
anew, which helps while you write one. Plugins with settings, such as the
cutoffs of hydrogen bonds or the radius of a binding pocket, have a section
in Settings (:ref:`ug-settings`); **Settings...** opens it at the selected
plugin. A change applies at once and is remembered.

If a plugin keeps MILES from working, start it without plugins from
packages and the plugin folder:

.. code-block:: console

   $ miles --no-plugins

then switch the plugin off in Help > Plugins... or remove it.

The repository has two examples. ``examples/plugins/principal_axes.py``
adds **View > Principal Axes**, which draws the axes of inertia as arrows
with their moments and names the kind of rotor in the status line.
``examples/plugins/selection_info.py`` has
a Selection Info panel (**View > Selection Info**, :kbd:`Option Cmd I`) with
the formula, mass, centre of mass and radius of gyration of the selection,
two commands under **Tools > Selection Info**, and reading and writing of
``.atoms`` files. :ref:`dg-plugins` in the developer guide explains how to
write a plugin.
