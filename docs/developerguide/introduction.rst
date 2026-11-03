Introduction
############

What MILES is
******************

MILES is a program for building, editing and looking at molecular
structures. It covers small organic and inorganic molecules, metal complexes
and whole proteins. It does not set up or run calculations: its output is a
structure file or a picture, which other programs then use.

The program is written in Python. All numerical work is done with NumPy
arrays, the interface uses Qt through PySide6, molecules are drawn with
OpenGL shaders, and Open Babel reads and writes the file formats. There is no
compiled code of its own.

This guide describes how the program is put together. It assumes you know
Python and are comfortable reading NumPy code; it does not assume you know Qt
or OpenGL. The :doc:`../api/index` documents every public function, and the
source itself is meant to be read alongside this guide.

Design principles
*****************

A few principles explain most of the decisions you will meet in the code.

**The interface stays out of the way.** The canvas fills the window. Tools
sit in a small floating rail, everything else lives in the menu bar,
popovers that appear next to the control that opened them, and keyboard
shortcuts. New features follow the same pattern instead of adding toolbars or
docked panels.

**It is fast at every size.** Building a methane and rotating a 10 000-atom
protein should both feel immediate. Anything whose cost grows with the number
of atoms is vectorised with NumPy, and spatial queries use a uniform grid
instead of comparing all pairs. Drawing is done on the GPU, where the cost of
a frame hardly depends on the molecule.

**It is general.** No feature is tied to one simulation or quantum chemistry
package. Formats are handled in one place, through Open Babel, and every
format is treated the same way.

**It has few dependencies.** The runtime needs PySide6, NumPy, PyOpenGL and
Open Babel and nothing else. A short NumPy implementation is preferred to a
new package: the neighbour search, the Clean force field, the point group
detection and the cartoon meshes are all written here.

**The code is meant to be read.** Modules are small and have one job, public
functions have NumPy-style docstrings, and the package passes ``ruff`` and
``mypy``. :doc:`conventions` lists the rules in full.

Choice of language and libraries
********************************

Python makes the program easy to change and to extend, and with NumPy the
numerical parts run at compiled speed as long as they work on whole arrays.
The few loops that remain in hot paths run over bonds or residues, not over
atom pairs.

**PySide6** is the official Python binding of Qt. It provides the windows,
menus and widgets, native file dialogs, the clipboard, high-resolution
displays and trackpad gestures on every platform, and it ships type stubs that
``mypy`` uses.

**PyOpenGL** gives access to OpenGL 3.3 core shaders. Atoms and bonds are
drawn as impostors (see :doc:`rendering`), so the renderer needs nothing from
OpenGL beyond instanced quads and a fragment shader. When OpenGL is not
available, a QPainter renderer draws the same scene more slowly.

**Open Babel** reads and writes about 150 chemistry formats, parses and writes
SMILES, generates 3D coordinates from a connection table and provides the
MMFF94 and UFF force fields. Using one library for all of this keeps the file
layer to a few hundred lines and means every format behaves the same way.
Other toolkits were considered and left out to keep the dependency list
short: RDKit, for instance, would duplicate most of what Open Babel already
provides.

A tour of the source tree
*************************

.. code-block:: text

   miles/
       core/
           molecule.py      the molecule as columns: atoms, coordinates, bonds, properties
           elements.py      periodic table data and typical bond lengths
           geometry.py      vectors, rotations, internal coordinates, VSEPR templates
           building.py      placing and attaching fragments, hydrogens, deleting, joining
           internal.py      setting bond lengths, angles and dihedrals; rigid moves
           perception.py    bonds and bond orders from coordinates; atom types
           forcefield.py    the small force field behind Clean
           chem.py          Open Babel: conversion, SMILES, 3D building, MMFF94
           fragments.py     element fragments, rings, groups, biomolecules, ligands
           protein.py       residue classes, atom groups, secondary structure
           symmetry.py      point groups and symmetrisation
           stereo.py        R/S and E/Z labels, inverting and flipping
           surface.py       molecular surfaces and the values that colour them
           interactions.py  binding pockets and non-covalent interactions
           mutation.py      replacing amino acid side chains
           biopolymer.py    peptides, DNA and RNA from sequences
           volume.py        volumetric data (cube files) and isosurfaces
           vibrations.py    vibrational modes, animation frames, IR spectra
           neighbors.py     grid-based neighbour search
           document.py      documents: molecule, display, frames, undo, change events
           display.py       per-atom representations and focus-mode hiding
           jobinfo.py       what travels with a structure besides its atoms
           data/fragments.json
       files/
           common.py        JobInfo and FileFormatError
           formats.py       the format registry, pasted text
           openbabel.py     every format Open Babel knows
           session.py       MILES session files (.miles)
       render/
           scene.py         display style and the scene arrays
           picking.py       atom and bond under the cursor
           camera.py        orthographic camera and smooth camera moves
           cartoon.py       protein and nucleic acid ribbons
           surface.py       colouring molecular surfaces for drawing
           theme.py         light and dark colours
           shaders.py       GLSL sources
           gl.py            the OpenGL renderer and picture export
           painter.py       the QPainter renderer and thumbnails
       ui/
           app.py           the miles command; making the application and window
           mainwindow/      the window: its tabs, files and dialogs
           menus.py         the menu bar, declared as data
           document.py      documents with Qt signals
           builder.py       current tool and fragment
           stage.py         the canvases and the panels floating over them
           canvas/          the molecule view (see the Interface chapter)
           selection_panel.py, rail.py, timeline.py, sequence_strip.py,
           focus_bar.py, command_palette.py, fragment_popover.py,
           symmetry_popover.py, dialogs.py, export_dialog.py, console.py,
           plugin_panel.py, plugins_dialog.py, controls.py, icons.py,
           clipboard.py, tasks.py, background.py, lookup.py, text.py
       __init__.py          the version and launch()
       interactive.py       the window from IPython, Jupyter and scripts
       api.py               the scripting interface
       scripting.py         the console's namespace and short commands
       plugin.py            what plugins import: the stable surface for them
       plugins/             finding and loading plugins, what they can add:
                            discovery.py, context.py, settings.py,
                            layers.py, manager.py
       testing.py           pytest fixtures for testing plugins
       builtin_plugins/     features written as plugins, one module or
                            package each, with the code only they use:
           ramachandran/    __init__.py (registration), panel.py, angles.py
           grids/           __init__.py, panel.py, combine_dialog.py
           ...              (the Writing Plugins chapter lists them all)

:doc:`architecture` explains how these modules fit together. The chapters
after it describe the data model, building, chemistry, files, rendering and
the interface, then how to extend the program and how it is tested.
