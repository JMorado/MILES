.. _ug-shortcuts:

Shortcut Reference
##################

This chapter lists every mouse action, key and menu command. On Linux and
Windows, read Ctrl for Cmd and Alt for Option. A summary is available in the
program under **Help > Shortcuts** (:kbd:`Cmd /`): the mouse actions and
keys, then the shortcut of every menu command, taken from the menus
themselves, so commands of plugins appear there and those of plugins that
are switched off do not.

Mouse
*****

On the canvas
=============

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Drag empty space
     - Rotate the view (or translate or zoom, depending on the mouse mode)
   * - Shift-drag, middle-drag
     - Move the view
   * - Cmd-drag, two-finger rotate
     - Spin the view about the line of sight
   * - Scroll, pinch, right-drag
     - Zoom
   * - Right-click, Control-click
     - Context menu
   * - Double-click empty space
     - Fit to window (select and delete tools)

With the select tool
====================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Click an atom or bond
     - Select it
   * - Shift-click, Cmd-click
     - Add to or remove from the selection
   * - Click empty space
     - Clear the selection
   * - Double-click
     - Select the residue, then the chain, then everything
   * - Drag a selected atom
     - Move the selection
   * - Option-drag
     - Rotate the selection
   * - Option-Shift-drag
     - Move the selection
   * - Option-scroll (2 to 4 atoms selected)
     - Step the bond length, angle or dihedral

With the build tool
===================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Click empty space
     - Place the current fragment
   * - Click a hydrogen
     - Replace it with the fragment
   * - Click a heavy atom
     - Change its element (element fragments) or attach to it (other
       fragments, terminal atoms only)
   * - Click a bond
     - Cycle single, double, triple
   * - Drag from an atom
     - Grow the fragment in that direction
   * - Drag onto another atom
     - Bond the two atoms (closes rings)
   * - Option-drag, Option-Shift-drag
     - Rotate or move the molecule under the cursor

With the delete tool
====================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Click an atom
     - Delete it with its hydrogens
   * - Click a bond
     - Delete the bond

In the selection panel
======================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Drag the value sideways
     - Scrub the value
   * - Click the value, type, Return
     - Set the value
   * - Up and Down arrows (Shift: ten times more)
     - Step the value
   * - Drag or scroll the dial
     - Change the value
   * - Option while dragging
     - Five times finer
   * - Pin button (2 to 4 atoms)
     - Keep the measurement on screen, or remove it

In the sequence strip
=====================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Click a residue
     - Select it
   * - Drag across residues
     - Select a stretch of residues
   * - Shift-click, Cmd-click
     - Add residues to the selection
   * - Scroll
     - Move along the sequence

Keyboard
********

These keys work when the cursor is over the canvas.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Key
     - Effect
   * - Element symbol: :kbd:`c`, :kbd:`n`, :kbd:`o`, :kbd:`c` :kbd:`l` ...
     - Choose that element (build tool)
   * - :kbd:`1` to :kbd:`9`
     - Choose a bonding type of the current element
   * - :kbd:`r`, :kbd:`t`, :kbd:`s`
     - Rotate, translate or scale mouse mode (with the build tool, :kbd:`s`
       is sulfur)
   * - :kbd:`Esc`
     - Clear the selection, or switch to the select tool
   * - :kbd:`Delete`, :kbd:`Backspace`
     - Delete the selection, or the atom or bond under the cursor
   * - Arrow keys
     - Turn the view by 6 degrees
   * - :kbd:`+`, :kbd:`=`, :kbd:`-`
     - Zoom in and out
   * - :kbd:`.`, :kbd:`,`
     - Next and previous frame (structures with several frames)
   * - :kbd:`Space`
     - Play or pause the frames

In the Python console
=====================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Key
     - Effect
   * - :kbd:`Return`
     - Run the line (a block runs after an empty line)
   * - :kbd:`Up`, :kbd:`Down`
     - Earlier and later commands
   * - :kbd:`Tab`
     - Complete the name before the cursor
   * - :kbd:`Control L`
     - Clear the output (Ctrl+L on Linux and Windows)

Menus
*****

Shortcuts marked "standard" follow the conventions of the operating system;
the macOS keys are shown.

File
====

.. list-table::
   :widths: 40 60

   * - New
     - :kbd:`Cmd N` (standard)
   * - Open
     - :kbd:`Cmd O` (standard)
   * - Open from PDB
     - :kbd:`Shift Cmd O`
   * - Save
     - :kbd:`Cmd S` (standard)
   * - Save As
     - :kbd:`Shift Cmd S` (standard)
   * - Export Image
     - :kbd:`Shift Cmd E`
   * - Export Movie
     - no shortcut
   * - Export 2D Drawing
     - no shortcut
   * - Run Script
     - no shortcut
   * - Close
     - :kbd:`Cmd W` (standard)

Edit
====

.. list-table::
   :widths: 40 60

   * - Undo, Redo
     - :kbd:`Cmd Z`, :kbd:`Shift Cmd Z` (standard)
   * - Cut, Copy, Paste
     - :kbd:`Cmd X`, :kbd:`Cmd C`, :kbd:`Cmd V` (standard)
   * - Duplicate
     - :kbd:`Cmd D`
   * - Copy SMILES
     - :kbd:`Shift Cmd C`
   * - Copy 2D Drawing
     - no shortcut
   * - Select All
     - :kbd:`Cmd A` (standard)
   * - Select...
     - :kbd:`Cmd F`
   * - Clean Geometry
     - :kbd:`Cmd L`
   * - Point Group
     - :kbd:`Cmd G`
   * - Add Hydrogens
     - :kbd:`Option Cmd H`
   * - Rebond
     - no shortcut
   * - Move to Origin
     - no shortcut
   * - Clear Measurements
     - no shortcut
   * - Atom List
     - :kbd:`Shift Cmd A`
   * - Settings
     - :kbd:`Cmd ,` (standard; in the MILES menu on macOS)

Tools
=====

.. list-table::
   :widths: 40 60

   * - Command Palette
     - :kbd:`Shift Cmd P`
   * - Fragments and SMILES
     - :kbd:`Cmd K`
   * - Build from Sequence, Generate Conformers, Scan Coordinate, Superpose
     - no shortcut
   * - Vibrations
     - :kbd:`Shift Cmd V`
   * - Properties
     - :kbd:`Cmd I`
   * - Prepare Protein, Show Binding Pocket
     - no shortcut
   * - Ramachandran Plot
     - :kbd:`Option Cmd R`
   * - Select, Build, Delete
     - :kbd:`Cmd 1`, :kbd:`Cmd 2`, :kbd:`Cmd 3`

View
====

.. list-table::
   :widths: 40 60

   * - Style: Ball and Stick, Tube, Wire, Space Fill, Cartoon
     - :kbd:`Option Cmd 1` to :kbd:`Option Cmd 5`
   * - Labels, Cartoon Colors
     - no shortcut
   * - Mouse Mode: Rotate, Translate, Scale
     - :kbd:`r`, :kbd:`t`, :kbd:`s` on the canvas
   * - Hydrogens
     - :kbd:`Shift Cmd H`
   * - Hydrogen Bonds
     - :kbd:`Shift Cmd B`
   * - Stereo Labels, Depth Fade
     - no shortcut
   * - Surface: Show, Hide, colours, See-Through
     - no shortcut
   * - Show Only Selection
     - :kbd:`Shift Cmd F`
   * - Hide Selection
     - no shortcut
   * - Show Everything
     - :kbd:`Shift Cmd U`
   * - Python Console
     - :kbd:`Option Cmd P`
   * - Appearance
     - no shortcut
   * - Fit to Window
     - :kbd:`Cmd 0`
   * - Center
     - :kbd:`Cmd E`
   * - Reset Orientation
     - :kbd:`Shift Cmd R`
   * - Look Along X, Y, Z
     - no shortcut

Help
====

.. list-table::
   :widths: 40 60

   * - Shortcuts
     - :kbd:`Cmd /`
   * - Plugins
     - no shortcut
