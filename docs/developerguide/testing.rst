Testing
#######

The tests live in ``tests/`` and run with ``pytest`` in a few seconds. They
fall into two groups.

Core tests
**********

Most tests exercise ``core`` and ``files`` without a display. They build
molecules with the same functions the interface uses and check geometric and
chemical facts rather than exact coordinates:

* building: attaching fragments, ring closure, hydrogens, bonding types and
  element changes give the expected formula, bond orders and bond lengths;
* internal coordinates: every move mode changes only the side it should,
  including inside rings;
* perception and Open Babel: SMILES round trips, conversion, rebonding of
  known structures;
* symmetry: point groups of reference molecules and symmetrisation;
* proteins: secondary structure and the residue groups;
* stereochemistry: R/S and E/Z labels of reference molecules, inversion and
  flipping, and the speed of labelling a large structure;
* analysis: pocket residues and interactions, superposition, mutations;
* surfaces: closed, outward-facing meshes of the right area, colouring
  values, and the time for a large structure;
* the neighbour grid: the same pairs as a brute-force search;
* files: reading and writing several formats, pasted text and sessions.

Prefer testing a property ("every carbon has four bonds", "the dihedral is
now 60 degrees") to comparing coordinates, which change whenever an
algorithm is improved.

Interface tests
***************

``tests/test_gui.py`` drives the real main window. ``tests/conftest.py`` sets
``QT_QPA_PLATFORM=offscreen``, so no window appears, and
``MILES_SOFTWARE``, so the QPainter canvas is used and no GPU is needed.
Each test gets a fresh window with its own settings file in a temporary
directory, so tests never change the user's preferences.

The tests click and type with ``QTest`` and drag by sending mouse events
directly to the canvas. A helper converts an atom's position to screen
coordinates, so a test can click a specific atom:

.. code-block:: python

   def test_build_by_clicking(window):
       canvas = build_ethane(window)
       assert canvas.mol.formula() == "C2H6"
       window.undo()
       assert canvas.mol.formula() == "CH4"

These tests cover building, typing elements, mouse modes, drawing bonds,
the selection panel, deleting, copy and paste, the menus, the fragment
popover, the context menu, picture export, opening a protein, focus mode,
pinned measurements, conformers and the timeline, sessions with their view,
surfaces, stereo labels, mutations with the sequence strip, superposition,
movie frames and the pocket view. They are smoke tests: they check that
each path works end to end, not every detail.

The OpenGL renderer
*******************

The OpenGL code cannot run offscreen. When you change ``render/gl.py`` or the
shaders, start the program and look at a small molecule, a protein cartoon
and an exported picture with ambient occlusion, in both themes. Running with
``--software`` shows the QPainter renderer for comparison.

Before finishing a change
*************************

.. code-block:: bash

   ruff check miles tests
   mypy
   pytest

All three must pass.
