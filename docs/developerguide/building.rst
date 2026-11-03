Building Algorithms
###################

This chapter explains how MILES decides where new atoms go. The code is
in ``core/geometry.py`` (coordination templates), ``core/building.py``
(fragments, hydrogens, deleting and joining) and ``core/internal.py`` (bond
lengths, angles and dihedrals). All of these functions edit the molecule
they are given; the interface wraps each call in an undo step.

Atom types and coordination templates
*************************************

Every atom the builder places has an intended bonding pattern: a geometry, a
number of sites and a valence (see :doc:`molecule_model`). The geometries are
the usual VSEPR shapes, defined in
:class:`~miles.core.geometry.Geometry`, and each has a template of unit
vectors in ``TEMPLATES``:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Geometry
     - Sites
     - Ideal angles (degrees)
   * - linear
     - 2
     - 180
   * - trigonal planar
     - 3
     - 120
   * - tetrahedral
     - 4
     - 109.47
   * - square planar
     - 4
     - 90, 180
   * - trigonal bipyramidal
     - 5
     - 90, 120, 180
   * - octahedral
     - 6
     - 90, 180

The builder offers each element as one or more bonding types, the
:class:`~miles.core.building.AtomType` records listed by
:func:`~miles.core.fragments.element_types`. Carbon, for example, can be
tetrahedral (four sites, valence 4), trigonal (three sites, valence 4, so one
double bond), linear or a bare atom without hydrogens. A bent oxygen has the
tetrahedral geometry with two sites, so its two bonds sit at the tetrahedral
angle and the missing sites stand for lone pairs.

Free sites
**********

To add bonds to an atom that already has some, the existing bond directions
are fitted onto the template, and the template vectors that were not used are
the free sites. :func:`~miles.core.geometry.free_sites` does this:

* With no existing bonds, all template vectors are free.
* With one bond, the template is rotated so its first vector lies along the
  bond. The remaining vectors are free sites, but their rotation about the
  bond is still undetermined; the caller chooses it (see below).
* With two or more bonds, every assignment of template vectors to the
  existing bonds is tried. For each, the Kabsch algorithm gives the rotation
  that best maps the template vectors onto the bonds, and the assignment with
  the smallest error wins. Templates have at most six vectors, so the number
  of assignments stays small.

Avoiding crowding
*****************

Whenever a free choice remains, such as the rotation of a methyl group or the
sites a new hydrogen takes, MILES chooses the least crowded option. The
crowding of a set of new positions :math:`\mathbf{p}_a` with respect to the
existing atoms :math:`\mathbf{r}_b` is the soft repulsion computed by
:func:`~miles.core.building.clash_score`:

.. math::

   E = \sum_{a} \sum_{b} \max(d_{ab}, 0.5)^{-6}

The steep power means the closest contact dominates, and the floor at 0.5
angstrom keeps the sum finite. Options are tried in a fixed order and an
option replaces the current best only when it is lower by more than
:math:`10^{-12}`, so the result does not depend on rounding noise.

Placing a fragment
******************

:func:`~miles.core.building.place_fragment` adds a free-standing copy of
a fragment. The fragment is rotated about a chosen origin (its centroid, or
the central atom of an element fragment) and moved so that the origin lands
on the clicked point. The rotation is the one shown in the fragment preview,
expressed relative to the current view, so a fragment appears on the canvas
the way it looked in the preview.

Attaching a fragment
********************

Most building replaces a hydrogen with a fragment.
:func:`~miles.core.building.attach_fragment` takes the hydrogen ``h``
(bonded to atom ``x``), a fragment, and the fragment's attachment atom: a
terminal atom of the fragment that stands for the bond and is discarded. Its
neighbour, the anchor, is the atom that bonds to ``x``.

1. **Bond order.** The spare valence of each side is its valence minus its
   current bond-order sum minus one (for the new single bond), not counting
   the bond being replaced. The new bond's order is one plus the smaller
   spare valence, rounded down and capped at a triple bond. Attaching a
   trigonal carbon to another trigonal carbon therefore makes a double bond,
   while attaching it to a tetrahedral carbon makes a single bond.

2. **Bond length.** The anchor is placed along the direction from ``x`` to
   ``h`` at the typical length for the two elements and the bond order, from
   :func:`~miles.core.elements.bond_length`.

3. **Orientation.** The fragment is rotated so that the direction from its
   anchor to its attachment atom points back at ``x``.

4. **Rotation about the new bond.** The fragment may still turn about the
   bond. If both ``x`` and the anchor are trigonal, only the two rotations
   that make their planes coincide are tried, which keeps conjugated systems
   flat. Otherwise rotations from 0 to 350 degrees in 10-degree steps are
   tried. The least crowded rotation wins.

5. **Indices.** The anchor takes the index of the replaced hydrogen; the
   other fragment atoms are appended. Existing indices do not change.

A hydrogen without a parent atom (a lone H) is simply replaced by a
free-standing fragment.

Growing towards the cursor
**************************

Dragging from an atom grows the fragment in the direction of the drag.
:func:`~miles.core.building.best_hydrogen_towards` projects the atom and
each of its hydrogens to the screen and returns the hydrogen whose bond
points most nearly along the drag; the fragment is then attached there as
above. If the drag ends on another atom instead, the two atoms are joined
(see below).

Hydrogens
*********

:func:`~miles.core.building.missing_hydrogens` gives the number of
hydrogens an atom lacks: the smaller of its free sites and its spare valence
rounded down. :func:`~miles.core.building.site_directions` chooses where
they go:

1. The existing bond directions are fitted onto the atom's geometry, giving
   the free sites. If the atom has more bonds than its geometry allows, the
   geometry for the total number of bonds is used instead.
2. With exactly one existing bond, the free sites can turn about it. If the
   atom and its neighbour are both trigonal, the two rotations that keep them
   coplanar are tried; otherwise 36 rotations are tried.
3. For every rotation and every choice of the required number of free sites,
   the hydrogens are placed at the X-H bond length and the crowding is
   computed. The least crowded combination wins.

:func:`~miles.core.building.add_hydrogens` applies this to every atom,
or to the selected ones. Atoms read from files do not have a reliable
intended valence, so their type is first inferred again with standard
valences (see :doc:`chemistry`); a carbon with three single bonds then gets
its fourth hydrogen. New hydrogens inherit the residue data and the
representation of their parent, so hydrogens added to a protein belong to the
right residue and stay hidden inside a cartoon.

:func:`~miles.core.building.add_valence` adds one hydrogen even to a
saturated atom by raising its number of sites and choosing the geometry for
the new coordination number.

Changing an element
*******************

:func:`~miles.core.building.mutate_atom` gives an atom another element
or bonding type. Bonds to heavy neighbours are kept, but their orders are
reduced until they fit the new type's valence, highest orders first. Then the
atom's hydrogens are rebuilt: existing hydrogens are moved onto the new free
sites, extra ones are added and surplus ones deleted, so most indices
survive.

Deleting and joining
********************

:func:`~miles.core.building.delete_atoms` removes the given atoms
together with the terminal hydrogens of every deleted heavy atom, so
deleting a carbon does not leave its hydrogens floating.

:func:`~miles.core.building.join_atoms` bonds two atoms, for example to
close a ring. A hydrogen given as one end stands for its parent and is
removed. A saturated end gives up the hydrogen that points most nearly
towards the other end. After the bond is made, the molecule is usually far
from a sensible shape, so the interface relaxes the joined molecule with a
fold: its 3D coordinates are rebuilt from the bonds by Open Babel and aligned
back onto the old position, with a soft force-field relaxation as fallback
(see :doc:`chemistry`).

Bond lengths, angles and dihedrals
**********************************

``core/internal.py`` changes one internal coordinate by moving part of the
molecule rigidly:

* :func:`~miles.core.internal.set_bond_length` slides one side along the
  bond axis;
* :func:`~miles.core.internal.set_angle` turns one side about the normal
  of the angle's plane, through the vertex;
* :func:`~miles.core.internal.set_dihedral` turns one side about the
  central bond. The change is taken modulo 360 degrees, so the shortest turn
  is used.

Which atoms move is chosen by a :class:`~miles.core.internal.MoveMode`:
the side of the last atom (the default), the side of the first atom, both
sides by half the change each, or only the last atom with its terminal
neighbours. A side is found with a graph search: starting from the end atom
and blocking the pivot atom, :meth:`~miles.core.molecule.Molecule.component`
returns everything on that side of the bond. If that search reaches another
neighbour of the pivot, the two ends are connected through a ring and there
is no separate side; then only the end atom and its terminal neighbours move.

The selection panel calls these functions as you drag the dial, recording a
single undo step for the whole drag.
