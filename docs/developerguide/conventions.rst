Code Conventions
################

These are the rules the code base follows. Read them before changing
anything; a change that breaks one of them should come with a good reason.

What the program is
*******************

* A general molecule builder and viewer. It builds and edits structures; it
  does not set up or run calculations.
* Nothing is tied to one simulation or quantum chemistry package. Formats are
  supported generally, through one file layer.
* Small molecules and whole proteins are both first-class: a feature that
  works on benzene must also behave on a 10 000-atom PDB entry.

Interface
*********

* Minimal, modern and fast. The canvas fills the window; tools sit in a small
  floating rail and everything else lives in menus, popovers and keyboard
  shortcuts. No docks, no toolbars full of icons.
* Drawing happens on the GPU (OpenGL core profile, impostor spheres and
  cylinders), with a QPainter fallback for machines without OpenGL.
* Interaction must stay smooth with thousands of atoms. Anything that scales
  with the number of atoms is vectorised or uses a spatial grid.
* Light and dark themes follow the system.

Dependencies
************

* Keep them minimal. The runtime dependencies are PySide6, NumPy, PyOpenGL
  and Open Babel, and nothing else without a strong case.
* Open Babel is the layer for file formats, SMILES and force fields, and
  there is no second chemistry toolkit. A reader or writer of our own is
  allowed only where Open Babel cannot do the job and the gap matters to
  users: it cannot read the format at all (Gaussian input), drops what the
  file holds (Molden orbitals, the steps and modes of some program
  outputs, the B-factors and alternate locations of PDB records), or
  writes it wrongly (PDB files). Such code lives in ``miles/files``, never
  in the core, reads only the part Open Babel misses, and comes with
  tests built from small hand-made inputs.
* Prefer a short NumPy implementation to a new package (the neighbour grid
  replaces a KD-tree, for example).
* Environments are created with micromamba from conda-forge, not with venv
  and pip: ``environment.yml`` to run the program, ``environment-dev.yml`` to
  develop it.

Python
******

* Python 3.12 or newer.
* No ``from __future__ import annotations``. Use built-in generics and
  ``X | None``; when an import is only needed for typing, guard it with
  ``TYPE_CHECKING`` and quote the name.
* Type hints on every function; the package passes ``mypy``.
* The package passes ``ruff check`` with the configuration in
  ``pyproject.toml``.
* NumPy-style docstrings on modules, public classes and public functions,
  with Parameters, Returns and Attributes sections as needed. A one-line
  docstring is enough for a small helper.
* Fixed sets of options are enums (``StrEnum`` or ``IntEnum``), not loose
  string constants.
* Records are dataclasses, with ``slots=True`` where many are created.
* Matrices and point sets keep their mathematical names (``R``, ``P``,
  ``Vt``); this is the one exception to lower-case variable names. Qt's own
  camelCase method and signal names are the other.
* Qt event handlers carry ``@override``.
* Source files are plain ASCII; other characters are written as escapes.

Structure
*********

* Small modules with one job each, in four layers: ``core``, ``files``,
  ``render`` and ``ui`` (see :doc:`introduction`). A layer never imports a
  later one; ``core`` and ``files`` never import Qt.
* Use a design pattern where it removes duplication or branching, not for its
  own sake: a registry for file formats, one class per canvas tool, the menu
  bar declared as data.
* Pick data structures and algorithms that fit the problem size: adjacency
  lists for bonds, a grid for contacts, Tarjan's algorithm for ring bonds.
* Readable first: short functions, clear names, no clever one-liners.
* Keep the package lightweight: no dead code and no unused options. A
  clean-up keeps every feature exactly as it was.

Tests
*****

* Core logic is tested without a display. Interface tests drive the real
  window offscreen with the QPainter canvas and never touch the user's
  settings.
* ``ruff``, ``mypy`` and ``pytest`` all pass before a change is finished.

Documentation
*************

* The README stays short: installation, how to run, and a link to the
  documentation.
* The user guide, this developer guide and the API reference live in
  ``docs/`` and are published with GitHub Pages.
* A new feature comes with its section in the user guide, and a new module or
  algorithm with its description here.

Distribution
************

* The package installs a ``miles`` command that opens files, PDB codes,
  compound names or SMILES strings given on the command line.
* Releases are published as a conda package; the GitHub workflows test every
  change on Linux and macOS and build the package on every tag. Windows is
  not tested yet.
* Standalone executables for Linux, macOS and Windows are planned, built with
  PyInstaller by GitHub Actions.
