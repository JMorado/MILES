Setting Up a Development Environment
####################################

Getting the source
******************

Clone the repository and work from its top directory. Everything below
assumes that directory is the current one.

Creating the environment
************************

MILES uses conda environments from conda-forge, created with
micromamba. There are two: ``environment.yml`` holds only what the program
needs to run, and ``environment-dev.yml`` adds the tools for development:
pytest, ruff and mypy for checks, Sphinx and Furo for the documentation, and
rattler-build and python-build for packaging. Create the development one:

.. code-block:: bash

   micromamba env create -f environment-dev.yml
   micromamba run -n miles-dev pip install -e . --no-deps --no-build-isolation

The second command installs the package in editable mode, so changes to the
source take effect the next time the program starts. ``--no-deps`` keeps pip
from installing PyPI copies of packages that conda already provides; in
particular, conda's ``openbabel`` and the ``openbabel-wheel`` package named in
``pyproject.toml`` must not both be installed.

To run commands without typing ``micromamba run`` each time, activate the
environment in your shell with ``micromamba activate miles-dev``.

Add a package to ``environment.yml`` only if the program needs it to run,
and then to ``environment-dev.yml``, ``pyproject.toml`` and
``recipe/recipe.yaml`` as well; tools go only into ``environment-dev.yml``.

Running the program
*******************

The editable install creates the ``miles`` command:

.. code-block:: bash

   miles                     # empty window
   miles ligand.sdf 1abc.pdb # open files
   miles 4hhb caffeine       # download a PDB entry, build a compound by name
   miles --software          # draw with QPainter instead of OpenGL
   miles --no-plugins        # without plugins from packages or the plugin folder

``python -m miles`` does the same. To make the command available in every
shell, link it into a directory on your ``PATH``:

.. code-block:: bash

   ln -s ~/.local/share/mamba/envs/miles-dev/bin/miles ~/.local/bin/miles

Setting the environment variable ``MILES_SOFTWARE`` has the same effect
as ``--software``. It is useful for testing the fallback renderer and for
running the interface on machines without a GPU.

Checking your changes
*********************

Three tools must pass before a change is finished:

.. code-block:: bash

   ruff check miles tests
   mypy
   pytest

``ruff`` checks style, import order, naming, docstrings and a set of common
mistakes; its configuration is in ``pyproject.toml``. ``mypy`` checks types
across the whole package, including the calls into Qt. ``pytest`` runs the
tests, which take a few seconds; :doc:`testing` describes them.

Building the documentation
**************************

The documentation is written in reStructuredText and built with Sphinx and
the Furo theme, both in the development environment:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html
   open docs/_build/html/index.html        # xdg-open on Linux, start on Windows

The second command opens the site in your browser. The API reference is
generated from the docstrings, so a function documented in the source
appears there automatically.

The site is published with GitHub Pages; :doc:`releasing` explains how to
set it up.
