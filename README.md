<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/_static/images/miles-lockup-dark.svg">
    <img src="docs/_static/images/miles-lockup-light.svg" alt="MILES" width="440">
  </picture>
</p>

<p align="center">
  <b>Molecular Improvisation &amp; Ligand Exploration Studio</b><br>
  <i>Kind of cool. Kind of molecular. Improvise your structures.</i>
</p>

<p align="center">
  <a href="https://github.com/JMorado/MILES/actions/workflows/ci.yml"><img src="https://github.com/JMorado/MILES/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://jmorado.github.io/MILES/"><img src="https://github.com/JMorado/MILES/actions/workflows/docs.yml/badge.svg" alt="Documentation"></a>
  <a href="https://anaconda.org/jmorado/miles-studio"><img src="https://img.shields.io/conda/vn/jmorado/miles-studio.svg?label=conda" alt="Conda version"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-BSD--3--Clause-blue.svg" alt="License: BSD-3-Clause"></a>
  <img src="https://img.shields.io/badge/python-3.12%2B-blue.svg" alt="Python 3.12+">
  <img src="https://img.shields.io/badge/platform-macOS%20%7C%20Linux-lightgrey.svg" alt="Platforms: macOS and Linux">
</p>

MILES is a fast molecule builder and 3D viewer for proteins, ligands, binding
pockets and small molecules. Build molecules by clicking and dragging, open
about 150 chemistry file formats, download proteins by their PDB code and
export pictures for papers.

## Installation

With [micromamba](https://mamba.readthedocs.io) (conda and mamba work the
same way):

```bash
micromamba create -n miles -c jmorado -c conda-forge miles-studio
```

To install from a copy of this repository instead:

```bash
micromamba env create -f environment.yml
micromamba run -n miles pip install . --no-deps
```

## Running

```bash
micromamba run -n miles miles               # empty window
micromamba run -n miles miles ligand.sdf    # open a file
micromamba run -n miles miles 4hhb caffeine # a PDB entry, a compound by name
```

From a Jupyter notebook or IPython, the window opens in the running session
and stays usable between cells:

```python
import miles

app = miles.launch("ligand.sdf")
app.current.style(mode="tube")
app.current.image()          # a picture of the view below the cell
```

A script run with plain `python` ends with `app.wait_until_closed()`.

## Documentation

The [user guide, developer guide and API reference](https://jmorado.github.io/MILES/) are online.
They are built from `docs/`; to build the website locally:

```bash
micromamba env create -f environment-dev.yml
micromamba run -n miles-dev pip install -e . --no-deps --no-build-isolation
micromamba run -n miles-dev sphinx-build -b html docs docs/_build/html
```

then open `docs/_build/html/index.html`.

## Licence

MILES is released under the [BSD 3-Clause licence](LICENSE). It uses
[Open Babel](https://openbabel.org), which is licensed under the GPL version 2;
a package that bundles the two together must follow the GPL as a whole.
