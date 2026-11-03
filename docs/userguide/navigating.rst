.. _ug-navigating:

The Window and the View
#######################

Moving the View
***************

The view is an orthographic camera that looks at the structure. Moving it
never changes the coordinates of the atoms. You can turn, move and zoom the
view with the mouse, a trackpad or the keyboard.

With the mouse
==============

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Action
     - Effect
   * - Drag
     - Rotate about the centre of the view.
   * - Cmd-drag
     - Spin about the line of sight, following the cursor around the centre.
   * - Shift-drag, or drag with the middle button
     - Move the view.
   * - Drag with the right button
     - Zoom: up zooms in, down zooms out.
   * - Scroll
     - Zoom, keeping the point under the cursor fixed.
   * - Double-click empty space
     - Fit the structure to the window (select and delete tools).

A drag only starts once the mouse has moved a few pixels, so a slightly
shaky click still counts as a click.

Mouse modes
===========

A plain drag rotates the view by default. Like VMD, MILES lets you
change this with a single key, so you can move or zoom with the plain
mouse button too. The mode decides what a plain left drag does on empty
space, or anywhere a tool does not start a drag of its own:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Key
     - Mode
     - A plain drag
   * - :kbd:`r`
     - Rotate
     - turns the view (the default).
   * - :kbd:`t`
     - Translate
     - moves the view.
   * - :kbd:`s`
     - Scale
     - zooms: dragging right or up zooms in.

The modes are also in **View > Mouse Mode**, and the mode is shared by all
tabs. The cursor changes shape in the translate and scale modes (a move
cursor and a vertical resize cursor), and a message names the new mode.

Shift-drag, Cmd-drag, middle-drag and right-drag keep their meaning in every
mode. So do the drags that tools start themselves: dragging from an atom with
the build tool still grows the structure, and dragging a selected atom with
the select tool still moves the selection.

The mode keys share the keyboard with element symbols
(:ref:`ug-element-keys`):

* :kbd:`r` and :kbd:`t` always switch the mode, since no element has these
  symbols. If the next letter, typed quickly, completes a two-letter symbol
  such as Br, Ru or Ti, that element is chosen and the previous mode comes
  back.
* :kbd:`s` switches to the scale mode, except with the build tool, where it
  chooses sulfur.
* A mouse click ends a symbol, so typing :kbd:`s`, dragging and then typing
  :kbd:`r` gives the rotate mode, not strontium.

With a trackpad
===============

Pinching zooms around the cursor and the two-finger rotate gesture spins the
view about the line of sight. Two-finger scrolling zooms, like a mouse
wheel.

With the keyboard
=================

The arrow keys turn the view in steps of 6 degrees, and :kbd:`+` (or
:kbd:`=`) and :kbd:`-` zoom in and out. The canvas has the keyboard focus
whenever the cursor is over it, unless you are typing in a text field.

.. _ug-panels:

Panels
******

Some commands open a panel over the view: the grid and orbital panel, the
vibrations, the Ramachandran plot and the properties. Each can be folded
and moved, so it never hides the part of the structure you are looking at.

* The arrow at the left of a panel's top row folds the panel to that row
  (for the grid panel: the orbital or grid, the isovalue and the buttons)
  and opens it again.
* Drag a panel by its background or its title to move it anywhere in the
  window; it stays where you leave it, also when you switch tabs, and is
  kept inside the window when the window gets smaller. Double-click its
  background to put it back in its usual place.

Where each panel is and whether it is folded are remembered the next time
MILES starts.

If something goes wrong inside MILES (a bug), a message says what
happened, with details you can copy into a bug report, and the window
stays open so you can save your work.

Fitting and Orienting
=====================

The **View** menu has commands that move the camera smoothly to a new
position:

**Fit to Window** (:kbd:`Cmd 0`)
   Zoom so the whole structure fills most of the window. Panels along the
   top and bottom edges (the sequence, the timeline, the grid panel, the
   Python console) are left clear: the structure fits in the space between
   them.

**Center** (:kbd:`Cmd E`)
   Put the centre of the structure in the middle of the free space,
   keeping the zoom and the orientation.

**Reset Orientation** (:kbd:`Shift Cmd R`)
   Return to the default orientation and fit.

**Look Along X**, **Look Along Y**, **Look Along Z**
   Look down one of the coordinate axes.

To rotate about a particular molecule instead of the centre of the whole
structure, right-click one of its atoms and choose **Rotate Around This
Fragment**.

When you build, the view zooms out by itself if new atoms would end up
outside the window. It never zooms in on its own.

.. _ug-focus:

Showing Part of a Structure
***************************

In a large structure it often helps to see only the part you are working
on. Focus mode hides atoms without deleting them:

**View > Show Only Selection** (:kbd:`Shift Cmd F`)
   Hide every atom that is not selected, and zoom to what remains.

**View > Hide Selection**
   Hide the selected atoms and clear the selection.

**View > Show Everything** (:kbd:`Shift Cmd U`)
   Show all hidden atoms again and fit the view.

The same commands are in the context menu when they apply. While atoms are
hidden, a small bar at the top right of the canvas says how many atoms are
shown, for example "Showing 120 of 3,400 atoms", and has a **Show
Everything** button.

Hidden atoms are not drawn, cannot be picked and are left out of surfaces.
Fitting, centring and resetting the view use the visible atoms only. Hiding
and showing are edits, so they can be undone, and session files remember
which atoms are hidden (:ref:`ug-sessions`). This is separate from the
**Hidden** representation (:ref:`ug-representations`): an atom hidden in
focus mode comes back with its own representation when you show everything.

Documents and Tabs
******************

Each structure is a document with its own tab. **File > New**
(:kbd:`Cmd N`) opens an empty document, and opening a file or downloading a
PDB entry adds a tab for it. When you open a structure while the current
tab is empty and unchanged, the new structure takes its place.

The tab bar appears only when more than one document is open. A dot after
a tab's name means the document has unsaved changes. **File > Close**
(:kbd:`Cmd W`) closes the current tab; if it has unsaved changes you are
asked whether to save them first. Quitting asks the same question for every
modified document.

Each document has its own undo history (see :ref:`ug-undo`) and its own view
and display style. The current tool and fragment are shared by all tabs.

Appearance
**********

MILES has a light and a dark theme. By default it follows the
appearance of the operating system and switches when the system does.
**View > Appearance** lets you choose **Match System**, **Light** or
**Dark**; the choice is remembered.

.. _ug-settings:

Settings
********

**Settings** (:kbd:`Cmd ,`) is in the MILES menu on macOS and at the
bottom of the **Edit** menu on Linux and Windows. Every change applies at
once and is remembered, so the window has no OK button.

General
   **Appearance**, the same choice as **View > Appearance**, and
   **Autosave**: every 1, 2, 5 or 10 minutes, or off
   (:ref:`ug-autosave`).

New Documents
   How documents opened from now on are drawn: the style of atoms without
   a representation of their own (Ball and Stick, Tube, Wire or Space
   Fill), the labels, the cartoon colours, whether hydrogens are shown and
   whether distant atoms fade. The **View** menu still changes the
   document in front, and a session opens with the view it was saved with.

Plugins
   Plugins with settings have a section each: the distance and angle
   cutoffs of **Hydrogen Bonds**, the radius of **Binding Pocket** and the
   length of the dipole arrow of **Charges and Dipole**, and those of any
   plugin you add (:ref:`ug-python`).

Finding Commands
****************

Press :kbd:`Shift Cmd P` (or Tools > Command Palette) and type part of a
command's name: ``surf`` finds the surface commands, ``add h`` finds Add
Hydrogens. The palette searches every menu; use the arrow keys to choose and
:kbd:`Return` to run the command, or :kbd:`Esc` to close it. Commands that
cannot be used at the moment are left out.

What the Window Tells You
*************************

Status line
   The bottom left corner shows the formula of the structure, its molar mass
   and, for structures of up to 500 atoms, its point group (found with a
   tolerance of 0.1 angstrom).

Hover text
   The bottom right corner describes whatever is under the cursor. For an
   atom it shows its label (element symbol and number, such as ``C12``), the
   element name, its bonding geometry and its coordinates; for an atom in a
   residue it shows the residue name, number and chain and the atom name. For
   a bond it shows the two atoms, the bond order and the length.

Messages
   Short messages appear at the top of the canvas for a moment, for example
   after saving a file, changing the element or adding hydrogens.

Bars
   Some bars appear only when they have something to show: the sequence
   strip at the top for proteins and nucleic acids (:ref:`ug-sequence`), the
   focus bar at the top right while atoms are hidden (:ref:`ug-focus`), and
   the timeline at the bottom for structures with several frames
   (:ref:`ug-timeline`).

Atom numbers
   Atoms are numbered from 1 in the order they appear in the structure. The
   same numbers are used in labels, in the selection panel and in the atom
   list.
