# Examples

Files to try the viewers of MILES. Open them with File > Open, by dropping
them on the window, or from the command line, e.g.
`miles examples/alanine_dipeptide_md.xyz`.

| File | Shows |
| --- | --- |
| `alanine_dipeptide_md.xyz` | A trajectory: 200 frames of molecular dynamics (MMFF94, 300 K). The timeline appears; press Space to play. |
| `qm/benzene_homo.cube` | A molecular orbital on a grid: the HOMO of benzene, a pi orbital with positive and negative lobes. |
| `qm/chloromethane_density.cube` | The electron density of chloromethane. Click Add Grid... in the grid panel and choose `chloromethane_potential.cube` to colour the density surface by the electrostatic potential. |
| `qm/chloromethane_potential.cube` | The electrostatic potential of chloromethane, on the same grid as the density. |
| `qm/formaldehyde.molden` | Molecular orbitals and vibrations of formaldehyde. The HOMO is drawn on opening and other orbitals are chosen in the grid panel; Tools > Vibrations shows the infrared spectrum and animates the modes. |
| `qm/nitrogen_dioxide.molden` | Molecular orbitals of nitrogen dioxide, an open-shell molecule: alpha and beta orbitals are listed separately. |

## A plugin

`plugins/selection_info.py` is an example plugin. It adds a Selection Info
panel (View > Selection Info) with the formula, mass, centre of mass and
radius of gyration of the selected atoms, two commands under
Tools > Selection Info, and reading and writing of `.atoms` files (one
`symbol x y z` line per atom).

To try it, open Help > Plugins..., click Open Plugin Folder, copy the file
into that folder and restart MILES. Or point MILES at this folder
for one run:

```bash
MILES_PLUGIN_DIR=examples/plugins micromamba run -n miles-dev miles
```

`plugins/principal_axes.py` is a smaller example of drawing: arrows along
the axes of inertia, with the moments as labels.

To share a plugin, make it a package and declare it as an entry point in
the `miles.plugins` group of its `pyproject.toml`; installing the
package into the environment of MILES is then enough. `plugins/template`
is such a package to copy, with tests that run without a screen:

```bash
cd examples/plugins/template
micromamba run -n miles-dev python -m pytest          # its tests
micromamba run -n miles-dev pip install -e . --no-deps   # then start MILES
```

The developer guide chapter "Writing Plugins" (`docs/developerguide/plugins.rst`)
describes everything a plugin can do.

## How the files were made

All structure files here are made by the scripts in this folder and are covered by
the project's licence.

- `make_trajectory.py` runs the molecular dynamics with Open Babel's MMFF94
  and needs only the development environment:

  ```bash
  micromamba run -n miles-dev python examples/make_trajectory.py
  ```

- `make_qm_examples.py` makes the files in `qm/` with
  [PySCF](https://pyscf.org): each molecule is optimised with
  B3LYP/def2-SVP, then the grids, orbitals, harmonic frequencies and
  infrared intensities are computed. PySCF is not a dependency of
  MILES, so the script runs in an environment of its own:

  ```bash
  micromamba create -n miles-examples -c conda-forge python=3.12 pyscf geometric numpy
  micromamba run -n miles-examples python examples/make_qm_examples.py
  ```
