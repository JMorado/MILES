Continuous Integration and Releases
###################################

Three GitHub Actions workflows in ``.github/workflows`` test the code,
publish this documentation and release new versions.

Workflows
*********

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Workflow
     - Runs on
     - Does
   * - ``ci.yml``
     - every push to ``main`` and every pull request
     - Runs the tests on Linux and macOS, ruff and mypy on Linux,
       and builds and tests the conda package.
   * - ``docs.yml``
     - every push to ``main``
     - Builds this documentation, treating warnings as errors, and deploys
       it to GitHub Pages.
   * - ``release.yml``
     - every tag starting with ``v``
     - Checks the version, builds the Python packages and the conda package,
       uploads the conda package to anaconda.org if configured, and creates
       a GitHub release with the files attached.

All three use the development environment from ``environment-dev.yml``,
so they run with the same versions you use locally.

Setting up the repository
*************************

These steps are needed once, after the repository is created on GitHub.

The website
===========

1. Push the code to the ``main`` branch.
2. On GitHub, open **Settings > Pages** and set **Source** to
   **GitHub Actions**.
3. Open **Actions > Documentation** and choose **Run workflow**, or push
   any commit to ``main``.

When the workflow finishes, the site is live at
``https://<user>.github.io/<repository>/``, and every later push to
``main`` updates it. GitHub Pages is free for public repositories.

To use your own domain, add it under **Settings > Pages > Custom domain**
and create the DNS record GitHub shows there.

A conda channel on anaconda.org
===============================

The release workflow can upload the conda package to your own channel, so
people can install a release as soon as it is tagged:

1. Create an account on `anaconda.org <https://anaconda.org>`_.
2. Under your account settings, create an access token with permission to
   upload packages.
3. In the GitHub repository, open **Settings > Secrets and variables >
   Actions**. Add a secret ``ANACONDA_API_TOKEN`` with the token, and a
   variable ``ANACONDA_OWNER`` with your anaconda.org user or organisation
   name.

Users then install with:

.. code-block:: bash

   micromamba create -n miles -c <owner> -c conda-forge miles-studio

Without the secret and the variable the upload step is skipped and the
release still happens.

Making a release
****************

1. Set the new version in ``miles/__init__.py`` (``__version__``) and
   in ``recipe/recipe.yaml`` (``version`` under ``context``). The package
   version is read from ``__version__``, so these are the only two places.
2. Commit, push and wait for the CI workflow to pass.
3. Tag the commit and push the tag:

   .. code-block:: bash

      git tag v0.3.0
      git push origin v0.3.0

The release workflow checks that the tag matches both versions, builds the
source distribution and wheel with ``python -m build``, builds and tests
the conda package with rattler-build, uploads it to anaconda.org, and
creates the GitHub release.

You can build the conda package locally the same way:

.. code-block:: bash

   rattler-build build --recipe recipe/recipe.yaml -c conda-forge

The recipe
==========

The package is called ``miles-studio``, the name it will have on PyPI and
conda-forge too, because ``miles`` was already taken on PyPI. It is not
published on PyPI yet; releases attach the wheel to the GitHub release.
Everything else keeps the short name: the module (``import miles``), the
``miles`` command and the ``miles.plugins`` entry point group.

``recipe/recipe.yaml`` describes the conda package. It is a ``noarch:
python`` package, so one build serves every platform, and it declares the
``miles`` command as an entry point so conda creates it on Windows as
well. Its tests import the main modules in a fresh environment and run
``miles --version``.

The recipe turns off pip's dependency check: ``pyproject.toml`` names Open
Babel by its PyPI name, ``openbabel-wheel``, so that ``pip install`` works,
while conda provides it as the package ``openbabel``, which pip does not
recognise.

Publishing on conda-forge
*************************

conda-forge is the community channel most scientific Python users install
from. Packages there are built and updated by its own infrastructure, so
once MILES is accepted, users can run ``conda install -c conda-forge
miles-studio`` without adding any other channel.

Submitting the package is a one-time pull request:

1. Make a release as described above, so a source archive exists at
   ``https://github.com/<user>/<repository>/archive/refs/tags/v<version>.tar.gz``.
2. Compute its checksum:

   .. code-block:: bash

      curl -sL <archive URL> | shasum -a 256

3. Fork `conda-forge/staged-recipes
   <https://github.com/conda-forge/staged-recipes>`_ and add
   ``recipes/miles-studio/recipe.yaml``: a copy of ``recipe/recipe.yaml`` with
   the ``source`` section replaced by the archive

   .. code-block:: yaml

      source:
        url: https://github.com/<user>/<repository>/archive/refs/tags/v${{ version }}.tar.gz
        sha256: <checksum>

   and your GitHub user name under ``extra: recipe-maintainers``.
4. Open the pull request and answer the reviewers.

When it is merged, conda-forge creates a ``miles-studio-feedstock``
repository with you as maintainer and publishes the package. From then on,
a bot opens a pull request on the feedstock for every new release tag;
merging it publishes the new version.
