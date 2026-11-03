Extending MILES
####################

This chapter walks through the most common extensions. Each section lists
the files to change and shows the code involved.

Adding a ring or a functional group
***********************************

Rings and functional groups are defined in ``core/fragments.py`` and built
the first time they are used.

A ring is an entry in ``_RINGS``: a key, the name shown in the popover and a
function that returns the molecule and its attachment hydrogen. The helpers
``_six`` (an aromatic six-membered ring from six atom specifications),
``_five``, ``_polycyclic`` (fused benzene rings from hexagon centres) and
``_ring`` (any ring from coordinates, atom specifications and bond orders)
cover most cases:

.. code-block:: python

   _RINGS = {
       ...
       "triazine": ("1,3,5-Triazine", lambda: _six([_N, _C, _N, _C, _N, _C])),
   }

A functional group is a recipe in ``_RECIPES``: the name, the root atom's
element and bonding type, and a list of growth steps. Each step is
``(parent, element, bonding type)``, where ``parent`` indexes the atoms
grown so far and 0 is the root. The group is grown atom by atom with the
same code the builder uses, then cleaned:

.. code-block:: python

   _RECIPES = {
       ...
       "ethoxy": ("Ethoxy", "O", "Bent", [(0, "C", "Tetrahedral"), (1, "C", "Tetrahedral")]),
   }

Add the key to ``_GROUP_ORDER`` to place it in the popover.

Amino acids, nucleobases and ligands come from ``core/data/fragments.json``.
Each entry has a key, a name, the SMILES it was made from, the attachment
atom, symbols, coordinates, bonds and, for residues, PDB atom names and the
three-letter code.

Adding a file format
********************

Open Babel already reads and writes most formats, so a new format usually
needs nothing. To handle one yourself, write a reader and a writer and
register them for the extension in ``files/formats.py``:

.. code-block:: python

   def read_frames(path: str) -> tuple[Molecule, JobInfo]:
       ...

   def write_frames(path: str, mol: Molecule, info: JobInfo) -> None:
       ...

   REGISTRY.register("frm", read_frames, write_frames)

A plugin can register a format without changing ``files/formats.py``
(see :doc:`plugins`). Registered formats appear in the file dialogs
automatically. Raise
:class:`~miles.files.common.FileFormatError` with a readable message for
files you cannot handle.

Adding a command
****************

A command is a ``MainWindow`` method and an entry in ``MENUS``. If it
changes the molecule, do the work in a ``core`` function and call it through
the document, so it can be undone:

.. code-block:: python

   # core/building.py
   def remove_hydrogens(mol: Molecule) -> int:
       """Delete every terminal hydrogen; returns how many were removed."""
       ...

   # ui/mainwindow.py
   def remove_hydrogens(self) -> None:
       """Delete every hydrogen."""
       canvas = self.current_canvas()
       if canvas is not None:
           n = canvas.doc.edit(building.remove_hydrogens)
           canvas.toast(f"Removed {plural(n, 'hydrogen')}")

   # ui/menus.py, in the Edit menu
   Command("Remove Hydrogens", "remove_hydrogens"),

Write a test for the ``core`` function, and add the command to the user
guide.

A self-contained feature can instead be a plugin in
``miles/builtin_plugins``, which keeps it out of ``MainWindow`` and
shows other people how to write one (see :doc:`plugins`). Add the module to
:data:`~miles.plugins.discovery.BUILTIN` and to the entry points in
``pyproject.toml``.

Adding a tool
*************

A tool is a subclass of :class:`~miles.ui.canvas.tools.CanvasTool`.
Override only what differs from the default, which highlights hovered atoms
and does nothing on click:

.. code-block:: python

   class MeasureTool(CanvasTool):
       """Click two atoms to pin their distance."""

       def click(self, canvas, hit, pos, additive):
           ...

Then:

1. add a member to :class:`~miles.ui.builder.Tool`;
2. add the instance to ``TOOLS`` in ``ui/canvas/tools.py``;
3. add a button to the rail in ``ui/rail.py``;
4. add the tool to the ``tool`` choice in ``ui/menus.py`` with its shortcut.

Adding a representation
***********************

Representations are members of
:class:`~miles.core.molecule.Representation`. To add one:

1. add the member and its menu name in ``REPRESENTATION_NAMES``
   (``render/scene.py``);
2. give it a sphere radius in :func:`~miles.render.scene.atom_radius`
   and decide how its bonds are drawn in ``_bond_cylinders``;
3. if picking should behave differently, adjust ``Picker._pick_radii``.

The style menu, the context menu and the selection panel list
representations from the enum and the name table, so they pick up the new
one.

Adding a colour scheme
**********************

Themes are instances of :class:`~miles.render.theme.Theme`. A new theme
needs only new colours; :func:`~miles.render.theme.theme_for` and the
``Appearance`` choice in the View menu decide which one is used.
