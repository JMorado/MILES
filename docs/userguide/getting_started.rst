.. _ug-getting-started:

Getting Started
###############

Installing MILES
*********************

MILES is a Python package. It needs Python 3.12 or newer and four
libraries: PySide6 (Qt), NumPy, PyOpenGL and Open Babel. The simplest way to
get all of them is a conda environment from conda-forge. The repository
contains an ``environment.yml`` file for this; with `micromamba
<https://mamba.readthedocs.io>`_ the steps are:

.. code-block:: bash

   micromamba env create -f environment.yml
   micromamba run -n miles pip install . --no-deps

The first command creates an environment called ``miles`` with Python
and the four libraries. The second installs MILES itself into that
environment; ``--no-deps`` stops pip from installing a second copy of the
libraries conda already provided. If you want to change MILES itself,
use the development environment described in the developer guide instead.

conda and mamba accept the same commands. If you prefer pip alone, ``pip
install .`` in the repository installs MILES with its dependencies from
PyPI (Open Babel comes from the ``openbabel-wheel`` package).

The miles Command
======================

The installation creates a ``miles`` command inside the environment.
To run it from any terminal without activating the environment, link it
into a directory on your ``PATH``:

.. code-block:: bash

   ln -s ~/.local/share/mamba/envs/miles/bin/miles ~/.local/bin/miles

The path to the environment depends on where micromamba keeps its
environments; ``micromamba env list`` shows it.

.. _ug-install-app:

In the Dock, the Start Menu or the Application Menu
===================================================

Started from a terminal, MILES runs inside Python, so the system names
it after Python: the Dock, the taskbar and the application switcher show
"python" and its icon. One command adds MILES to the system under its own
name and icon:

.. code-block:: bash

   miles --install-app

What it makes depends on the system:

macOS
   ``MILES.app`` in ``~/Applications``, so MILES can be started from
   Launchpad, Spotlight and the Dock like any other application. To have
   the ``miles`` command show the name too when you start it from a
   terminal, link the command to the launcher inside the bundle:

   .. code-block:: bash

      ln -sfn ~/Applications/MILES.app/Contents/MacOS/MILES ~/.local/bin/miles

Linux
   A ``miles.desktop`` entry and the icon in ``~/.local/share``, so MILES
   appears in the application menu and its windows get the MILES icon in
   the dock and the task switcher.

Windows
   A MILES shortcut with the icon in the Start menu. The taskbar shows
   MILES with its own icon however it was started.

The entry starts MILES with the Python of the environment the command
was run in. Run ``miles --install-app`` again after moving or recreating
that environment; the old entry is replaced. A folder after the option,
as in ``miles --install-app ~/Desktop``, puts the entry there instead.

Running MILES
******************

Start the program with ``miles``. Anything you write after the command
is opened in its own tab:

.. code-block:: bash

   miles                       # an empty canvas
   miles ligand.xyz 1abc.pdb   # open files
   miles 4hhb                  # download a PDB entry
   miles caffeine "CC(=O)O"    # build from a name or a SMILES string

Each item is treated as follows:

* If a file with that name exists, it is opened.
* If it ends in a structure extension (``.xyz``, ``.pdb``, ``.sdf`` and so
  on) but the file does not exist, MILES reports that the file was not
  found.
* Otherwise the text is looked up: four characters starting with a digit
  are a PDB code, valid SMILES is built directly, and anything else is
  treated as a compound name and fetched from PubChem. Lookups run in the
  background, so the window opens at once; PDB codes and names need an
  internet connection.

These options change how the program runs:

``--software``
   Draw with QPainter instead of OpenGL. Use it on machines whose graphics
   drivers cause trouble. MILES also falls back to QPainter by itself
   when no suitable OpenGL context is available. Setting the environment
   variable ``MILES_SOFTWARE`` has the same effect as the option.

``--no-plugins``
   Start without plugins from installed packages and from the plugin folder
   (see :ref:`ug-plugins`), for when a plugin keeps MILES from
   working. The environment variable ``MILES_NO_PLUGINS`` has the same
   effect.

``--install-app [FOLDER]``
   Add MILES to the Dock, the Start menu or the application menu under its
   own name and icon, then exit (see :ref:`ug-install-app`).

``--version``
   Print the version and exit.

Your First Molecule
*******************

This section builds ethanol and saves it. It takes about a minute and shows
the main ideas of building.

1. Start MILES. The canvas is empty, the build tool is active and the
   current fragment is a tetrahedral carbon, shown on the build button of
   the tool rail.

2. Click anywhere on the canvas. A methane molecule appears: the carbon with
   its four hydrogens.

3. Click one of the hydrogens. It turns green while the cursor is over it,
   which means the current fragment will attach there. The click replaces
   the hydrogen with a new carbon, and you have ethane.

4. Type :kbd:`o`. The current fragment becomes a bent (sp3) oxygen; a short
   message at the top of the canvas confirms the change. Click a hydrogen on
   either carbon. It is replaced by an OH group: ethanol.

5. To grow the chain instead of clicking, press :kbd:`c` and drag outwards
   from a carbon atom. A dashed line follows the cursor; when you let go, a
   new carbon is attached to the hydrogen that points most nearly in the
   direction you dragged. Press :kbd:`Cmd Z` to undo it if you want to keep
   ethanol.

6. Press :kbd:`Cmd L` (or the clean button in the rail) to relax the
   geometry. Bond lengths and angles settle to sensible values; for a small
   organic molecule like this one, the MMFF94 force field is used.

7. Press :kbd:`Cmd S`. Because the document is new, a save dialog opens;
   choose a name such as ``ethanol.xyz``. The extension you type selects the
   format.

The status line at the bottom left shows the formula, the molar mass and
the point group of what you built. :ref:`ug-building` explains every way of
building in detail.

Example Files
*************

The ``examples`` folder of the repository has files to try the viewers:

``alanine_dipeptide_md.xyz``
   A short molecular dynamics trajectory: the timeline appears; press
   :kbd:`Space` to play it.

``qm/benzene_homo.cube``, ``qm/chloromethane_density.cube``
   A molecular orbital and an electron density, shown as isosurfaces.

``qm/chloromethane_potential.cube``
   The electrostatic potential of the same molecule. Open the density, click
   **Add Grid...** in the grid panel and choose this file to colour the
   density surface by the potential.

``qm/formaldehyde.molden``
   Molecular orbitals and vibrations of formaldehyde: choose orbitals in
   the grid panel, and open **Tools > Vibrations** for the infrared spectrum
   and the animated modes.

``qm/nitrogen_dioxide.molden``
   Molecular orbitals of nitrogen dioxide, an open-shell molecule with
   alpha and beta orbitals.

The trajectory comes from ``make_trajectory.py`` (Open Babel's MMFF94) and
the files in ``qm`` from ``make_qm_examples.py`` (PySCF, B3LYP/def2-SVP);
``examples/README.md`` explains how to run them.
