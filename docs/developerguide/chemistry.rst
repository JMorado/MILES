Chemistry
#########

This chapter covers the chemical perception of structures read from files,
the force fields behind Clean, the bridge to Open Babel, stereochemistry,
point groups, the neighbour grid, proteins, binding pockets, superposition,
mutations and molecular surfaces.

Perception
**********

Files such as XYZ carry no bonds, and many others carry bonds without
reliable orders. Open Babel perceives connectivity and aromaticity when it
reads a file, and MILES converts the result (see :doc:`files`). The
functions in ``core/perception.py`` provide the program's own perception,
used by Rebond and by the Clean force field.

**Bonds.** :func:`~miles.core.perception.perceive_bonds` connects two
atoms when their distance is less than 1.2 times the sum of their covalent
radii (0.85 angstrom for two hydrogens) and more than 0.4 angstrom. Candidate
pairs come from the neighbour grid, so the cost grows linearly with the
number of atoms. Orders of bonds that already existed are kept.

**Bond orders.** :func:`~miles.core.perception.assign_bond_orders`
compares each bond length with reference single, aromatic, double and triple
lengths for the element pair and takes the closest. Two corrections follow.
An "aromatic" length outside a ring, such as the partial double bond of a
peptide, becomes a single bond. Then any atom whose bond-order sum exceeds its
maximum valence has its least certain multiple bond reduced, preferring
aromatic bonds and bonds most stretched relative to their reference, until
the valence fits. This keeps a carbonyl's C=O while reducing a neighbouring
bond.

**Atom types.** :func:`~miles.core.perception.infer_atom_type` fills in
an atom's geometry, number of sites and valence from its bonds:

* the number of sites is the number of bonds, or, when completing an atom for
  Add Hydrogens, the number needed to reach the element's standard valence;
* two sites are linear for a triple bond, two double bonds or a straight
  angle, trigonal next to a pi system and tetrahedral otherwise;
* three sites are trigonal planar next to a pi system or when the bonds lie in
  a plane, and tetrahedral otherwise; nitrogens next to a pi system (amides,
  anilines, guanidines) count as conjugated;
* four sites on a metal with two straight angles are square planar,
  otherwise tetrahedral; five and six sites are trigonal bipyramidal and
  octahedral.

Clean
*****

Clean gives a structure a sensible geometry. :func:`~miles.core.forcefield.clean_geometry`
chooses between two force fields:

* **MMFF94**, through Open Babel, for organic molecules: structures of up to
  400 atoms whose elements are all covered by MMFF94 (H, C, N, O, F, P, S,
  Cl, Br, I), without residue data, and without delocalised bonds outside
  rings. A nitro group drawn with two order-1.5 bonds, for example, has no
  Kekule form that MMFF94 types correctly, so such structures use the
  built-in force field. MMFF94 gives realistic bond lengths, including the
  shortening of single bonds next to conjugated atoms.
* The **built-in force field** for everything else: metals and unusual
  coordination, large structures, and experimental structures.

Cleaning a selection uses the same choice for the region around it, with
the surrounding atoms held fixed.

The built-in force field
************************

The built-in force field (``core/forcefield.py``) is written for MILES.
It handles any element and coordination the builder can produce, including
metals, and it only needs to reproduce the intended geometry, not accurate
energies.

Terms
=====

The energy is a sum of five terms, each evaluated with vectorised NumPy over
index arrays built once per structure:

.. math::

   E = \sum_\text{bonds} k_b (r - r_0)^2
     + \sum_\text{angles} k_\theta (\cos\theta - \cos\theta_0)^2
     + \sum_\text{planar} k_p \sin^2\phi
     + \sum_\text{sp3-sp3} k_s (1 + \cos 3\phi)
     + \sum_\text{contacts} k_r \max(0, d_0 - d)^2

* **Bonds** pull towards the typical length for the elements and bond order.
* **Angles** pull towards the ideal angle of the central atom's geometry; for
  atoms whose missing sites are lone pairs (water, ammonia, hydrogen sulfide,
  phosphines) the experimental angle is used instead. For square planar,
  trigonal bipyramidal and octahedral centres, each pair of bonds targets the
  nearest of the geometry's ideal angles.
* **Planarity** torsions about bonds between two trigonal atoms keep
  conjugated systems flat.
* **Staggering** torsions on single bonds between tetrahedral atoms prefer
  staggered conformations.
* **Contacts** push apart non-bonded atoms of the same molecule that are
  closer than the sum of their van der Waals radii scaled down. For
  structures over 150 atoms, contacts are found with the neighbour grid
  within 6 angstrom.

The force constants are module-level constants (``K_BOND`` and so on).

Minimisation
============

:func:`~miles.core.forcefield.minimize` uses L-BFGS with a backtracking
line search and stops when the largest gradient component falls below a
tolerance or after a step limit. Frozen atoms have their gradient set to
zero.

Clean modes
===========

:func:`~miles.core.forcefield.clean_geometry` chooses between two modes:

* **Full.** A structure built in MILES is relaxed freely.
* **Gentle.** If any atom has residue data, the structure probably came from
  an experiment. Every such atom is tied to its starting position with a
  spring, and fewer steps are taken. Clashes and bad bond lengths are
  repaired without rebuilding the protein.

:func:`~miles.core.forcefield.clean_selection` relaxes only the selected
atoms. Every atom within 6 angstrom of the selection is included in the
calculation but held fixed, so a rebuilt side chain or ligand settles into
its surroundings.

:func:`~miles.core.forcefield.clean_atoms` with ``fold=True`` is used
after closing a ring. A chain whose ends were just bonded is far from any
ring shape, and a local minimiser cannot get there. Open Babel rebuilds 3D
coordinates from the bonds, the result is aligned onto the old position with
the Kabsch algorithm, and only if that fails a soft version of the force
field (no torsions or contacts) is minimised first.

The Open Babel bridge
*********************

``core/chem.py`` is the only module in ``core`` that imports Open Babel.

* :func:`~miles.core.chem.from_openbabel` and
  :func:`~miles.core.chem.to_openbabel` convert molecules. Open Babel
  numbers atoms from 1; everything outside this module uses 0-based indices.
  Aromatic bonds become order 1.5, and on the way back they are kekulised so
  that formats without aromatic bonds receive alternating single and double
  bonds. Residue data is kept only for formats where it is meaningful (PDB,
  mmCIF, GRO and a few more).
* :func:`~miles.core.chem.from_smiles` builds a 3D structure from SMILES,
  adds hydrogens and optimises it briefly. A ``*`` in the SMILES marks an
  attachment point; it becomes a hydrogen whose index is returned, so a
  SMILES fragment can be attached like any other.
* :func:`~miles.core.chem.to_smiles` writes canonical SMILES with
  stereochemistry taken from the 3D geometry.
* :func:`~miles.core.chem.minimize` minimises with MMFF94, or with UFF
  when MMFF94 lacks parameters, optionally with fixed atoms or with a
  chosen list of force fields.
* :func:`~miles.core.chem.reembed` rebuilds coordinates from the bonds
  alone.
* :func:`~miles.core.chem.conformers` runs Open Babel's genetic
  conformer search, minimises each conformer, drops those within 0.25
  angstrom RMSD of a lower-energy one and superposes the rest on the input.
  Open Babel prints its progress from C++; the search runs with the
  standard output and error streams redirected at the file-descriptor
  level so the terminal stays quiet.

Quantities along a trajectory
*****************************

``core/trajectory.py`` measures a whole stack of frames, shape
``(n_frames, n_atoms, 3)``, in one NumPy pass:
:func:`~miles.core.trajectory.measurement_series` gives a distance,
angle or dihedral in every frame with the same formulas as the single
measurements in ``core/geometry.py``, and
:func:`~miles.core.trajectory.unwrap_degrees` removes the 360-degree
jumps of a dihedral that turns past 180 degrees, for plotting.
:func:`~miles.core.trajectory.rmsd_series` superposes every frame on a
reference with the Kabsch method, batched: the singular values of the
3 x 3 covariance matrices of all frames give the RMSD directly, without
rotating any coordinates, and a reflection is excluded by flipping the
sign of the smallest singular value when the best orthogonal matrix has
determinant -1.

:func:`~miles.core.trajectory.align_frames` uses the same batched
singular value decomposition to rotate the frames themselves: each frame is
centred on its fitted atoms, rotated by ``V D U^T`` (with ``D`` flipping the
last axis when needed, so the matrix is a proper rotation) and moved to the
reference's centre. All atoms move with the fitted ones, so the internal
geometry of each frame is unchanged. **Tools > Align Frames** records an
undo step with ``Document.checkpoint(frames=True)``, which keeps a copy of
the frames, their labels and energies, so undo restores the trajectory and
not only the geometry on screen.

Substructure patterns
*********************

``core/patterns.py`` matches SMARTS patterns with Open Babel's
``OBSmartsPattern``. :func:`~miles.core.patterns.match_pattern`
converts the molecule with :func:`~miles.core.chem.to_openbabel`, which
keeps the atom order (Open Babel atom ``i + 1`` is MILES atom ``i``),
and returns the unique matches (``GetUMapList``) as
:class:`~miles.core.patterns.PatternMatch`. Open Babel's SMARTS
primitives for charge and hydrogen count need formal charges and hydrogen
counts that the molecule model does not store, so they are filled in per
connected fragment: fragments that contain hydrogens get the formal charges
of :func:`~miles.core.charges.formal_charges`, worked out from the
bonding, and fragments without any hydrogen get Open Babel's typical
implicit hydrogens (``OBAtomAssignTypicalImplicitHydrogens``), so a protein
without hydrogens bound to a ligand with hydrogens is handled correctly.
Matching a 2600-atom protein takes a few tens of milliseconds, so the
dialog matches again on every pause in typing.
:data:`~miles.builtin_plugins.select_atoms.examples.FUNCTIONAL_GROUPS` lists the
patterns the Select Atoms dialog offers as examples.

Relaxed scans
*************

:func:`~miles.builtin_plugins.scan.relaxed_scan.scan` steps a distance (two atoms), angle
(three) or dihedral (four) from a start value to a stop value. At each
point:

1. The coordinate is set by moving one side of the molecule with the
   functions of ``core/internal.py``, as when editing a measurement.
2. Open Babel minimises the whole structure with MMFF94 (UFF when MMFF94
   has no parameters) while an ``OBFFConstraints`` restraint holds the
   coordinate at its target: 300 steepest-descent steps, which cope with a
   geometry far from the minimum, then 500 conjugate-gradient steps. Open
   Babel's default restraint constant lets a dihedral drift by several
   degrees, so a much stiffer one (``RESTRAINT``) is used; the coordinate
   then stays within a few tenths of a degree.
3. The energy is evaluated again without the restraint, so it is the plain
   force-field energy of the relaxed geometry.

Each point starts from the geometry of the one before, so the scan follows
one continuous path; scanning back the other way can end on a different
path where the structure has several minima, as with any relaxed scan. The
result, a :class:`~miles.builtin_plugins.scan.relaxed_scan.ScanResult`, holds the measured
values, energies and frames.

Stereochemistry
***************

``core/stereo.py`` labels stereocentres R or S and double bonds E or Z from
the 3D coordinates, and can invert or flip them.

Priorities follow rule 1a of the Cahn-Ingold-Prelog rules: higher atomic
number ranks higher, and ties are broken by exploring the hierarchical
digraph sphere by sphere. A multiple bond adds duplicate atoms, and a ring
closure ends in a duplicate of the atom it returns to; an atom with aromatic
bonds gets one duplicate, which matches a Kekule structure. Isotopes and the
stereo rules 3 to 5 are not used, so pseudoasymmetric centres are not
labelled.

A centre is labelled when it has four single bonds, at most one hydrogen
and four branches of different priority. Its label comes from the sign of
the signed volume of its substituents ordered by priority, which also works
for distorted geometries. A double bond is labelled when it is not in a ring
of fewer than eight atoms and each end has two substituents of different
priority (or one substituent and a lone pair on N or P).

:func:`~miles.core.stereo.invert_center` turns the two smallest
substituents that are not in a ring through the centre by 180 degrees about
the bisector of their bonds, which swaps them without mirroring them.
:func:`~miles.core.stereo.flip_double_bond` turns the smaller side of a
double bond by 180 degrees about it. Labelling a 5 000-atom structure takes
a few tens of milliseconds, so the labels are simply recomputed after each
edit.

Errors from Open Babel are raised as :class:`~miles.core.chem.ChemError`.

Point groups
************

``core/symmetry.py`` finds the point group of a structure and can move atoms
onto exact symmetry positions.

1. The structure is centred on its centre of mass and its principal axes of
   inertia are computed. The shape of the inertia tensor (all moments
   different, two equal or all equal) decides which directions can be
   symmetry axes.
2. Candidate rotation axes and mirror normals are built from the principal
   axes and from shells of equivalent atoms: sets of atoms of the same
   element at the same distance from the centre. Midpoints and cross
   products of shell members give the candidate directions.
3. Each candidate operation is applied to the structure, and it is a
   symmetry operation if every atom lands within the tolerance of an atom of
   the same element.
4. The elements found (highest axis order, perpendicular two-fold axes,
   horizontal, vertical and dihedral mirrors, inversion centre) are
   classified into a Schoenflies symbol, and the group's operations are
   generated from a few exact generators.

:func:`~miles.core.symmetry.symmetrize` averages every atom over the
images of its equivalent atoms under all operations, which places the
structure exactly on the symmetry elements.

The neighbour grid
******************

``core/neighbors.py`` finds all pairs of points closer than a cutoff. Points
are binned into cubic cells as wide as the cutoff, so every close pair lies
in the same or adjacent cells. For each of the 27 cell offsets, all points
are matched against the points of the offset cell at once with array
operations, and the candidate pairs are filtered by distance. The cost grows
linearly with the number of points for a fixed density. Bond perception,
force-field contacts and the environment of a cleaned selection use it.

Binding pockets
***************

``builtin_plugins/binding_pocket/interactions.py`` finds the binding pocket of a ligand and the
non-covalent interactions between them.
:func:`~miles.builtin_plugins.binding_pocket.interactions.ligand_atoms` takes the ligand's hetero
residue (or its molecule), and
:func:`~miles.builtin_plugins.binding_pocket.interactions.pocket_atoms` every residue, apart from
water, with an atom within 5 angstrom of it.
:func:`~miles.builtin_plugins.binding_pocket.interactions.find_interactions` then applies the
usual geometric criteria:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Interaction
     - Criterion
   * - Hydrogen bond
     - N, O or S donor and N or O acceptor 2.4 to 3.5 angstrom apart; when
       the donor carries hydrogens, one of them must make a D-H...A angle
       above 120 degrees.
   * - Salt bridge
     - Centres of oppositely charged groups within 4.0 angstrom. Charged
       groups are found by atom name in amino acids (Lys, Arg, Asp, Glu,
       C-terminus) and by bonds elsewhere (carboxylates, phosphates,
       sulfates, quaternary nitrogens, amidines).
   * - Pi stacking
     - Flat five- or six-membered aromatic rings: parallel within 30 degrees
       and 5.5 angstrom, or T-shaped (60 to 90 degrees) within 6.0 angstrom.

Each interaction stores the atoms at its two ends, so the lines on screen
follow the atoms when they move.

Hydrogen bonds within a structure
*********************************

:class:`~miles.builtin_plugins.shared.hbonds.HydrogenBondFinder` applies the
hydrogen bond criterion to a whole structure. Its constructor does the work
that depends only on the bonds: which atoms donate (with their hydrogens)
and accept, the heavy neighbours of each, and the keys of polar pairs at
most :data:`~miles.builtin_plugins.shared.hbonds.HBOND_SEPARATION` bonds apart,
found by a short breadth-first search from each polar atom.
:meth:`~miles.builtin_plugins.shared.hbonds.HydrogenBondFinder.find` then takes any
coordinates: candidate pairs come from the neighbour grid, the excluded
pairs are dropped with ``numpy.isin`` on integer keys, and the D-H...A angle
(over up to four hydrogens per donor) and the 90 degree heavy-atom angles
are computed for all pairs at once. The structure is set up once and the
search is repeated for every frame of a trajectory; for a 4600-atom
helical protein setup takes about 5 ms and each search about 5 ms.
:func:`~miles.builtin_plugins.shared.hbonds.hydrogen_bonds` does both for one
geometry.

Roles come from atom names where bonds are ambiguous: amino acids (backbone
N except proline, the polar side chains, backbone O and OXT) and nucleobases
(ring NH and amino donors; carbonyl and ring nitrogen acceptors). Elsewhere
donors without modelled hydrogens are atoms with unfilled valence (a
hydroxyl, water, an NH or NH2 group) and nitrogens accept unless they carry
four substituents, or three and are conjugated with a neighbour (amides,
anilines, pyrrole-type nitrogens). A nitrogen donor with two heavy
neighbours and no modelled hydrogen is planar, so its hydrogen is placed on
the external bisector at 1.01 angstrom, as DSSP does; this keeps the
i to i+3 and i to i+2 contacts of a helix without hydrogens from counting.
Each :class:`~miles.builtin_plugins.shared.hbonds.HydrogenBond` records the donor,
the hydrogen (None when there is none in the model) and the acceptor.

Superposition
*************

``builtin_plugins/superpose/matching.py`` (the Superpose plugin) overlays
one structure on another. Atoms are paired
first, by the first method that applies:

1. **Proteins**: the one-letter sequences are aligned with
   ``difflib.SequenceMatcher`` and the CA atoms of identical aligned
   residues are paired.
2. **The same atoms in the same order**: atoms are paired by index.
3. **The same molecule numbered differently**: both bond graphs are refined
   together, Morgan-style, until atoms with the same label are equivalent;
   symmetric atoms are then separated one choice at a time, trying the
   alternatives for the first few and keeping the pairing with the lowest
   RMSD.

The Kabsch algorithm gives the rotation that best overlays the paired atoms
after centring, and :func:`~miles.builtin_plugins.superpose.matching.apply` moves the whole
mobile structure.

Mutations
*********

:func:`~miles.builtin_plugins.shared.mutation.mutate_residue` replaces an amino acid's
side chain. The new residue's template from the fragment library is
superposed on the residue's N, CA and C atoms, the old side chain is removed
and the template's side chain added, named with standard PDB atom names.
The chi1 torsion (and chi2 where the side chain allows) is then scanned in
30-degree steps for the conformation least crowded by atoms within 8
angstrom. Proline also gets its ring bond to N. Hydrogens are added only if
the residue had them.

Building from sequences
***********************

``core/biopolymer.py`` builds peptides, DNA and RNA from one-letter
sequences.

* **Peptides** are built by MILES itself: a poly-glycine backbone is
  laid out from internal coordinates (standard bond lengths and angles, the
  chosen phi and psi, omega 180 degrees), with NH2 and COOH termini and
  hydrogens. Every residue is then turned into its amino acid with
  :func:`~miles.builtin_plugins.shared.mutation.mutate_residue`, and side chains that
  still clash, such as a proline inside a helix, are relaxed with the
  backbone held fixed. Open Babel's own peptide builder was tried and left
  out because its helices place each carbonyl oxygen within 1 angstrom of
  the next nitrogen.
* **DNA** comes from Open Babel's FASTA reader, which builds a B-DNA double
  helix (or one strand) with residue data.
* **RNA** also comes from the FASTA reader, whose geometry has clashes and
  short O3'-P bonds; it is relaxed with MMFF94 four nucleotides at a time
  with the neighbouring nucleotides fixed, which keeps the helix and scales
  linearly with the length.

Protein preparation
*******************

``builtin_plugins/prepare_protein/preparation.py`` (the Prepare Protein
plugin) completes structures from the PDB.
:func:`~miles.builtin_plugins.prepare_protein.preparation.incomplete_residues` compares each amino
acid's side-chain atoms with its template and
:func:`~miles.builtin_plugins.prepare_protein.preparation.complete_side_chains` rebuilds the
incomplete ones. :func:`~miles.builtin_plugins.prepare_protein.preparation.add_hydrogens_at_ph`
removes all hydrogens and adds them again. For the standard amino acids the
number of hydrogens on each heavy atom follows from its name
(:data:`~miles.builtin_plugins.prepare_protein.preparation.SIDE_CHAIN_HYDROGENS`), adjusted for
the pH with the pKa values in :data:`~miles.builtin_plugins.prepare_protein.preparation.PKA`, for
the termini, for disulfide bridges and for bonds to metals or ligands; the
hydrogens are then placed by the builder's own geometry code
(:func:`~miles.core.building.add_hydrogens`). Open Babel's pH model is
used only for the remaining atoms. It used to protonate everything, but on
crystal structures it depends on perceived bond orders and implicit
hydrogen counts that are often inconsistent, and it left carboxylates and
histidine rings with two or more extra protons. Heavy atoms keep their
records and order; the new hydrogens inherit their parent's residue data
and get PDB-style names.

Molecular surfaces
******************

``core/surface.py`` computes smooth molecular surfaces.
:func:`~miles.core.surface.molecular_surface` follows the approach of
VMD's QuickSurf: each atom contributes a Gaussian density that equals the
iso-value exactly at its van der Waals radius, and the surface is the level
set of the summed density. An isolated atom is a sphere of its radius; where
atoms are close the densities add up and the surface bridges crevices, much
as a probe sphere would. Hydrogens contribute half a Gaussian, so the surface
looks similar with or without them.

The density is evaluated on a grid only near atoms (0.5 angstrom spacing
below 2 000 atoms, coarser above), each atom adding to the grid points
within its cutoff with vectorised NumPy. The isosurface is extracted with
Surface Nets: one vertex per grid cell the surface crosses, placed at the
mean of the edge crossings and moved onto the level set along the density
gradient; each crossed grid edge gives a quad between the four cells around
it. Distant groups of atoms are meshed on separate grids. Each vertex
remembers its nearest atom, which colours it.

:func:`~miles.core.surface.atom_property` gives the per-atom values for
colouring, between -1 and 1: hydrophobicity from the Kyte-Doolittle scale
for amino acids and from element and polarity otherwise, and charge from
charged residues, termini and common ligand groups.

Electrostatic potential on surfaces
***********************************

``core/potential.py`` computes the Coulomb potential of point charges,
:math:`\phi = 332.06 \sum_i q_i / (\varepsilon d_i)` in kcal/(mol e), for the
``POTENTIAL`` surface colouring. Two screening models are offered
(:class:`~miles.core.potential.Dielectric`): vacuum for small
molecules, where the map should resemble an electrostatic potential from a
quantum chemistry calculation, and a distance-dependent dielectric
:math:`\varepsilon = 4d` for proteins and nucleic acids
(:func:`~miles.core.potential.default_dielectric`). In vacuum the
potential of a protein is dominated by its net charge and by groups far
away; the screened potential reflects the groups near each part of the
surface, which is what the colouring is for. It is the default of the
Coulombic surface colouring in UCSF Chimera and ChimeraX, and like them the
potential is evaluated 1.4 angstrom outside the surface along its normals
(:data:`~miles.core.potential.PROBE_OFFSET`).

The charges of polymers are those of their ionisable groups
(:func:`~miles.core.surface.ionisable_charges`): crystal structures
have no hydrogens, and partial charges computed without them, or with
hydrogens placed by a simple protonation model, put large spurious charges
on hydroxyl, amide and backbone atoms. Small molecules use the charge model
chosen by the user (``surface_charges`` of the Molecular Surface plugin).

:func:`~miles.core.potential.electrostatic_potential` sums over all
charges with matrix products in single precision: squared distances come
from :math:`|p|^2 + |c|^2 - 2 p \cdot c` (one matrix product per chunk of
points), then the reciprocal (and square root, in vacuum) and a
matrix-vector product with the charges; coordinates are centred first so
single precision keeps its accuracy. 130 000 points and 4 500 charges take
about 0.8 s. Zero charges are dropped, so the ionisable-group charges of a
protein (a few hundred) take a few hundredths of a second. When the number
of point-charge pairs exceeds :data:`~miles.core.potential.EXACT_LIMIT`,
:func:`~miles.core.potential.surface_potential` evaluates the potential
exactly at one representative point (the mean) per small cell of the
surface, the cell size chosen to stay under the limit, and gives every
vertex the Gaussian-weighted average of the representatives within 1.5
cells. For hydrogenated hemoglobin with partial charges on all 9 000 atoms
this takes under a second, and the interpolated values correlate with the
exact ones at 0.99.

Volumetric data and isosurfaces
*******************************

:class:`~miles.core.volume.Volume` holds a grid read from a cube file:
an origin, three step vectors (which need not be orthogonal) and the values.
Open Babel reads the grids, but the number of points along each axis is not
reachable from its Python bindings; it is recovered from the far corner of
the grid.

:func:`~miles.core.volume.isosurface` extracts the surface where the
values equal an isovalue with the same Surface Nets code as the molecular
surface, generalised to any grid (:func:`~miles.core.surface.surface_nets`).
The grid is padded by one layer, so lobes cut by the grid's edge still come
out closed. Orbitals have both signs;
:func:`~miles.core.volume.signed_isosurfaces` returns the positive and
the negative lobes, each with normals pointing out of the lobe.
:meth:`~miles.core.volume.Volume.suggested_level` gives the first
isovalue: 0.02 for signed data and 0.002 otherwise, or a quarter of the
largest value when the data are smaller.

A second grid can colour an isosurface. :meth:`~miles.core.volume.Volume.sample`
interpolates a grid trilinearly at any points (clamped to the grid's edge),
so the two grids need not share spacing or origin;
:func:`~miles.core.volume.mapped_isosurface` returns the isosurface of
one volume with the other's values at its vertices, and
:func:`~miles.core.volume.symmetric_range` suggests the colour range
(the 95th percentile of the magnitudes, so a few extreme values near the
nuclei do not wash out the rest). The colours themselves are chosen in
``render/surface.py`` (:data:`~miles.render.surface.MAPPED_COLORS`,
the same red-white-blue scale as the charge colouring of molecular
surfaces), and :class:`~miles.render.surface.VolumeRequest` names the
mapped volume and the range.

:func:`~miles.builtin_plugins.grids.combine.combine_volumes` makes the difference or sum
of two grids (:class:`~miles.builtin_plugins.grids.combine.Combination`). When
:func:`~miles.builtin_plugins.grids.combine.same_points` finds the grids identical the
values are combined point by point; otherwise the second grid is sampled at
the points of the first, one plane at a time so large grids need little
memory. A grid's sign test, :attr:`~miles.core.volume.Volume.signed`,
normally looks at the values (both signs must reach 5% of the largest
magnitude), but :attr:`~miles.core.volume.Volume.both_signs` can
declare it: differences and spin densities take both signs by nature even
when one of them is small, and are drawn with two lobes at the density
isovalue. Sessions keep the marker.

Molecular orbitals
******************

``files/molden.py`` reads the basis set and the orbital coefficients of
Molden files, and ``core/orbitals.py`` evaluates the orbitals on grids,
which then go through the same isosurface code as cube files. :class:`~miles.core.orbitals.OrbitalSet`
holds the orbitals (energy, occupation, spin, symmetry label), the shells
of contracted Gaussians and the coefficient matrix;
:meth:`~miles.core.orbitals.OrbitalSet.label` names orbitals relative to
the frontier orbitals (HOMO-1, LUMO+2, with the spin for unrestricted
calculations) and :meth:`~miles.core.orbitals.OrbitalSet.evaluate`
returns a :class:`~miles.core.volume.Volume` named by
:meth:`~miles.core.orbitals.OrbitalSet.name`, which the interface uses
to find orbitals it has already computed.
:meth:`~miles.core.orbitals.OrbitalSet.density` sums the occupied
orbitals squared, weighted by their occupations, and
:meth:`~miles.core.orbitals.OrbitalSet.spin_density` the alpha minus
the beta orbitals of an unrestricted calculation, on the same grids
(:data:`~miles.core.orbitals.DENSITY_NAME`,
:data:`~miles.core.orbitals.SPIN_NAME`). Each occupied orbital is
evaluated in turn, which takes under a second for 40 occupied orbitals over
20 atoms. The spin density of the nitrogen dioxide example integrates to
1.01 electrons on the default grid; the total density only roughly gives
the number of electrons there (33.9 for 23), because the sharp core
densities are sampled at grid points that happen to lie near the nuclei,
and 22.1 on a 0.15 angstrom grid. The valence density that isosurfaces
show is unaffected.

Gaussians factorise along x, y and z, so each primitive is evaluated as an
outer product of three one-dimensional factors, and only inside the box
where it exceeds :data:`~miles.core.orbitals.CUTOFF` of its peak. A
grid of 60 points per side for a 20-atom molecule takes a few hundredths of
a second, fast enough to compute an orbital when it is chosen (on a worker
thread, so the window never waits).

Programs write Molden files in their own ways, and the reader
(:func:`~miles.files.molden.read_molden_orbitals`) works out which from the
orbitals themselves: every orbital must have norm 1 with the
exact overlap of the basis (:func:`~miles.core.orbitals.overlap_matrix`,
by the Obara-Saika recurrence, vectorised over the pairs of primitives of
each pair of shell kinds). The Molden format specifies contraction
coefficients for normalised primitives, but ORCA and xtb write them for
unnormalised ones, and cartesian d, f and g functions are normalised one
by one in the format but differently by CFOUR, by Psi4 up to 1.3.2 and by
Turbomole. The readings are tried in turn, starting with the format
itself, and the first that makes every orbital normalised (to
:data:`~miles.files.molden.NORM_TOLERANCE`) is kept; if none comes
within :data:`~miles.files.molden.NORM_FAILURE`, the orbitals are not
read at all rather than drawn wrong. On the Molden files of IOData's and
cclib's test data, from ORCA, Psi4, Molpro, CFOUR, Turbomole and Molden
itself, all orbitals come out orthonormal, to 1e-9 or better where the file
gives enough digits. ORCA also flips the sign of
the spherical functions with ``|m|`` of 3 or 4 in f, g and h shells (not
``|m| = 5``); files whose title says ``orca_2mkl`` are corrected. This was
checked on real ORCA output (CuH in the cc-pVQZ basis, from IOData's test
data): with the exact overlap of the basis, the orbitals are orthonormal to
1e-10 when corrected this way, as Psi4's orbitals for the same molecule are,
against 1e-3 uncorrected and 3e-6 when the h functions with ``|m| = 5`` are
flipped too. The angular tables come from
:func:`~miles.core.orbitals.solid_harmonics`, the usual recurrences for
real solid harmonics, which reproduce the hand-written d, f and g tables
used before exactly. Open Babel reads the structure; when its atoms
match the ``[Atoms]`` section, the molecule takes those coordinates so the
orbitals line up exactly.

Vibrations
**********

``core/vibrations.py`` holds the vibrational modes read from frequency
calculations (:class:`~miles.core.vibrations.VibrationalMode`: the
wavenumber, the infrared intensity and the displacement of every atom).
Open Babel provides them; its list of displacements can carry a spurious
extra vector, so only the first n are kept.
:func:`~miles.builtin_plugins.vibrations.modes.mode_frames` turns a mode into one period
of sinusoidal motion, scaled so the atom that moves most moves 0.4
angstrom, and :func:`~miles.builtin_plugins.vibrations.modes.ir_spectrum` sums a
Lorentzian band for each mode, weighted by its intensity.

Partial charges
***************

``core/charges.py`` provides partial charges and the dipole moment.
:func:`~miles.core.charges.partial_charges` runs one of Open Babel's
charge models (:class:`~miles.core.charges.ChargeModel`: Gasteiger or
MMFF94) on the structure. The molecule model has no formal charges, so
:func:`~miles.core.charges.formal_charges` works them out from the
bonds after kekulizing aromatic rings: four-bonded N and P are +1, two-bonded
N -1, one-bonded O, S and Se -1, three-bonded O, S and Se +1, four-bonded B
-1, and atoms without bonds take the charge of the common monatomic ion.
They are set on the Open Babel atoms before the model runs. MMFF94 perceives
charged groups itself; Open Babel's Gasteiger code leaves formal charges
out, so they are added to the atoms that carry them afterwards and both
models sum to the net charge. A model that fails (no parameters for some
element, or non-finite results) raises
:class:`~miles.core.chem.ChemError`. Open Babel's EEM, QEq and QTPIE
models were left out: on simple organic molecules they gave charges of the
wrong sign or size, or NaN for ions.

:func:`~miles.core.charges.dipole_moment` sums ``q r`` about the centre
of mass and converts e*angstrom to debye. For an ion the result depends on
the origin.

Charges and the dipole stored in calculation files are read by
:func:`~miles.files.openbabel.file_charges` into ``JobInfo.charges``
and ``JobInfo.dipole``. Open Babel marks a molecule whose charges came from
the file as having perceived partial charges; that flag must be read before
anything asks for a charge, because asking computes Gasteiger charges and
sets the same flag. The canvas caches the charges of each model until
atoms or bonds change (coordinates do not matter to either model) and uses
the reported dipole only while the coordinates are those it was read with.

Molecular properties
********************

:func:`~miles.builtin_plugins.properties.descriptors.molecular_properties` returns a
:class:`~miles.builtin_plugins.properties.descriptors.Properties` record. The formula and the
molecular weight come from the molecule itself (the same numbers as the
status line), the exact mass from Open Babel's table of isotope masses.
Everything else goes through an Open Babel molecule with the formal charges
of :func:`~miles.core.charges.formal_charges` set and stereochemistry
perceived from the 3D coordinates: canonical SMILES and the InChI and
InChIKey formats for the identifiers, and Open Babel's descriptors
``logP``, ``TPSA``, ``MR``, ``HBD``, ``HBA2`` and ``rotors``. Open Babel
has two acceptor counts; ``HBA2`` (all N and O atoms) is the one Lipinski
used. :attr:`~miles.builtin_plugins.properties.descriptors.Properties.lipinski_violations`
names the limits exceeded.

Two checks keep the numbers meaningful. Formal charges, SMILES and the
descriptors assume complete hydrogens, so
:func:`~miles.builtin_plugins.shared.molecules.hydrogens_missing` looks for carbon atoms
whose bond orders add up to less than 3.5 (aromatic bonds count 1.5); such
structures get their composition only. Identifiers of proteins are long and
useless and drug-likeness means nothing for them, so
:func:`~miles.builtin_plugins.shared.molecules.is_small_molecule` limits those to
structures with at most
:data:`~miles.builtin_plugins.shared.molecules.MAX_HEAVY_ATOMS` heavy atoms and fewer
than three amino acids or nucleotides. Small molecules take a few
milliseconds and a 3000-atom protein about 20 ms; the panel
(``builtin_plugins/properties/panel.py``) moves structures of more than 5000 atoms to
a worker thread.

2D drawings
***********

:func:`~miles.builtin_plugins.drawing.depiction.depiction_svg` draws a skeletal formula
with Open Babel on a copy of the structure. Formal charges are set and
stereochemistry is perceived from the 3D coordinates, then hydrogens on
carbon are made implicit (deleting a hydrogen in Open Babel does not
increment its neighbour's implicit hydrogen count, so that is done by
hand) and ``gen2D`` lays the molecule out. The writer's ``C`` option
leaves terminal carbons unlabelled, as in a skeletal formula; a molecule
with a single carbon, such as methane, still shows its label. One hydrogen
is kept on each stereocentre: Open Babel draws the wedge or hash bond only towards an
explicit atom, and without it the configuration would be lost from the
picture. The SVG writer nests the drawing in a fixed-size outer element,
which Qt's renderer skips (SVG Tiny allows no nested SVG), so the inner
element is lifted out with ``xml.etree`` and becomes a standalone document
sized to the drawing. ``builtin_plugins/drawing/files.py`` renders it to PNG with
``QSvgRenderer`` for files and the clipboard.

Proteins and nucleic acids
**************************

``core/protein.py`` knows the residue names of amino acids (including common
protonation variants), nucleotides and water. When a structure with a
polymer is opened, :func:`~miles.core.protein.default_representations`
gives the document's display the polymer residues as a cartoon and water
hidden, and leaves ligands and everything else in the view's style.

The cartoon needs secondary structure. Many PDB files carry it in ``HELIX``
and ``SHEET`` records; for those that do not,
:func:`~miles.core.protein.estimate_secondary_structure` assigns it from
the CA trace in the spirit of P-SEA. For each residue ``i`` it compares the
CA distances to residues ``i + 2``, ``i + 3`` and ``i + 4`` with the values
typical of each structure:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * -
     - i to i+2
     - i to i+3
     - i to i+4
   * - Helix
     - 5.5 +- 0.5
     - 5.3 +- 0.5
     - 6.4 +- 0.6
   * - Strand
     - 6.7 +- 0.6
     - 9.9 +- 0.9
     - 12.4 +- 1.1

Matching residues are marked, helices shorter than five residues and strands
shorter than three are dropped, and strands never override helices.

:func:`~miles.core.protein.group_atoms` returns the atoms of a residue,
a chain, the polymer, the ligands (hetero atoms that are not water or amino
acids) or the water, for the selection commands.

Backbone dihedrals
==================

``builtin_plugins/ramachandran/angles.py`` computes phi and psi for the Ramachandran plot in
two steps, so that playing a trajectory only repeats the cheap one.
:func:`~miles.builtin_plugins.ramachandran.angles.backbone_residues` finds the amino acid
residues with N, CA and C and links them through their peptide bonds: the
residue before or after in file order when its C and this N are bonded or
closer than 2 angstrom, otherwise any residue bonded to the N or C.
Residue numbers are never compared, so numbering gaps are harmless while
real chain breaks leave the angle undefined. Each
:class:`~miles.builtin_plugins.ramachandran.angles.BackboneResidue` keeps the four atoms
of its phi and psi and its class (general, glycine, proline or the residue
before a proline). :func:`~miles.builtin_plugins.ramachandran.angles.backbone_angles`
then evaluates all dihedrals for a set of coordinates with NumPy at once,
with NaN where an angle is undefined.

:func:`~miles.builtin_plugins.ramachandran.angles.classify` places a (phi, psi) pair in a
favoured, allowed or outlier region. The regions in
:data:`~miles.builtin_plugins.ramachandran.angles.REGIONS` are unions of ellipses per
residue class, compared with differences wrapped into -180..180 so they
continue across the edges of the plot. They are deliberately simple
approximations of the usual shapes, not contours fitted to a reference set
of structures.

The plot panel (:class:`~miles.builtin_plugins.ramachandran.panel.RamachandranPanel`)
keeps the result of ``backbone_residues`` until the structure, the bonds or
the frames change, and recomputes the angles from the coordinates shown on
every other change, including each step of trajectory playback.
