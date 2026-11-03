Files
#####

The ``files`` layer reads, writes, pastes and downloads structures. It does
not use Qt, so it can be used from scripts:

.. code-block:: python

   from miles.files.formats import read_file, write_file

   mol, info = read_file("ligand.sdf")
   print(mol.formula(), info.title)
   write_file("ligand.pdb", mol, info)

The format registry
*******************

``files/formats.py`` keeps a :class:`~miles.files.formats.FormatRegistry`
that maps file extensions to a reader and a writer. A reader takes a path and
returns a :class:`~miles.core.molecule.Molecule` with a
:class:`~miles.files.common.JobInfo`; a writer takes a path, a molecule
and a ``JobInfo``. Extensions without an entry fall through to Open Babel, so
the registry only lists the formats MILES handles itself, which today
are session files, Gaussian input (``.gjf``, ``.com``), which Open Babel
can write but not read (:mod:`miles.files.gaussian`), and PDB files for
writing (:mod:`miles.files.pdb`, see below).
:func:`~miles.files.formats.read_file` and
:func:`~miles.files.formats.write_file` look up the extension and
dispatch.

``JobInfo`` carries what travels with a structure besides its atoms: a
title, the total charge, the spin multiplicity when the file states one,
trajectory frames with their energies, grids of volumetric data,
vibrational modes, molecular orbitals, partial charges and the dipole
moment from a calculation and, for sessions, the saved view. Its ``notes``
say what the reader found in the file but could not use, such as orbitals
written in a way it does not know, in sentences meant for the user.

Errors are raised as :class:`~miles.files.common.FileFormatError`, a
subclass of ``ValueError``, with a message suitable for showing to the user.

Reading through Open Babel
**************************

:func:`miles.files.openbabel.read` reads any of the roughly 150 formats
Open Babel knows.

1. The format is chosen from the extension. Files ending in ``.log``,
   ``.out`` or ``.output`` (also gzipped) are output from many different
   programs, so the first part of the file is searched for the signature of
   each program Open Babel reads (Gaussian, ORCA, Q-Chem, NWChem, GAMESS-UK,
   GAMESS and Firefly, Jaguar, Molpro and MOPAC, in that order, since
   GAMESS-UK outputs also say "GAMESS"). Outputs of programs Open Babel
   cannot read, or reads wrongly (Psi4, Turbomole, DALTON, Molcas, xtb,
   Serenity and CP2K), are recognised the same way, so the error can say so
   (:func:`~miles.files.openbabel.unreadable_program`). Formats Open
   Babel can only write are not used for reading.
2. The first structure is read and converted with
   :func:`~miles.core.chem.from_openbabel`, which builds the columns and
   sets all bonds in one call
   (:meth:`~miles.core.molecule.Molecule.set_bond_arrays`). Open Babel
   perceives bonds for formats without them, and aromatic bonds become
   order 1.5. Isotopes are kept, and so are residue insertion codes. Formal
   charges are kept only from formats that state them
   (:data:`~miles.files.openbabel.CHARGE_FORMATS`: SDF, MOL, CML, SMILES),
   since Open Babel reports 0 both for neutral and for not given; elsewhere
   the ``formal_charge`` column stays unset. Two inputs are
   prepared first: Open Babel's GAMESS-UK reader aborts the whole program on
   the normal modes of a frequency calculation, so it only gets the text
   before them, and the atom labels of Molden files are written as element
   symbols, since Psi4's capitals (``CU``) become dummy atoms otherwise.
   When nothing is read, the error names the program or the format the
   file was read as ("read as GAMESS Input" for an ORCA input, which shares
   the ``.inp`` extension).
3. Further structures in the same file, or further conformers of the first,
   become trajectory frames when they are the same molecule: the same
   elements in the same order and, for formats whose bonds come from the
   file (:data:`~miles.files.openbabel.BOND_FORMATS`: SDF, MOL2, CML), the
   same bonds. That covers models in a PDB file, steps in an XYZ file, the
   conformers of an SDF file and the optimisation steps of a quantum
   chemistry output. Structures of other molecules, such as the rest of an
   SDF library, are counted but not read, and ``JobInfo.notes`` says how
   many there were. For output files
   the last geometry is shown, for everything else the first. Reading stops
   at the first empty structure, since the MOPAC and GAMESS input readers
   report empty structures forever at the end of a file.
   :func:`~miles.files.openbabel.frame_energies` attaches an energy to
   each frame (``JobInfo.energies``, kcal/mol) when the file has one for
   every frame. Open Babel's Gaussian, NWChem and GAMESS readers keep the
   energy of each optimisation step in ``OBConformerData``, the first two
   in kcal/mol and GAMESS in hartree, so the unit is looked up per format.
   For XYZ files, :func:`~miles.files.openbabel.xyz_comment_energies`
   reads the comment lines written by xtb (``energy: -5.07``), ORCA
   trajectories (``... E -76.3``), extended XYZ (``energy=-12.3`` with
   ``Properties=``, in eV) and CREST (the energy alone, in hartree). Every
   comment line must follow the same convention; other comments are left
   alone rather than guessed at.
4. :mod:`miles.files.outputs` reads from the text of quantum chemistry
   outputs what Open Babel leaves out. :func:`~miles.files.outputs.job_steps`
   finds the geometry of every step: Open Babel keeps only the first
   geometry of an ORCA scan, the last of a Q-Chem or GAMESS optimisation and
   all but the last point of a Gaussian unrelaxed scan. ORCA and Q-Chem
   print an energy with each step, which is kept too; Gaussian steps are
   only used when every geometry has its own energy, because an
   optimisation prints its last geometry twice. For GAMESS and Firefly a
   step starts with an ``NSERCH=`` line.
   :func:`~miles.files.outputs.charge_and_multiplicity` reads the total
   charge and multiplicity, which Open Babel takes from Gaussian outputs
   only (its value for the others is a guess, so it is not used).
   :func:`~miles.files.outputs.orca_modes` reads ORCA's normal modes,
   because Open Babel takes the infrared intensity from the wrong column of
   ORCA 5 and later and finds no modes at all when ORCA 6 labels them by
   symmetry. For formatted checkpoint files,
   :func:`~miles.files.outputs.checkpoint` reads the Mulliken charges,
   the dipole moment (in atomic units in the file, which Open Babel passes
   on unconverted), the normal modes and, for the files Q-Chem writes for an
   optimisation, every step.
5. Grids attached by Open Babel (cube files) become
   :class:`~miles.core.volume.Volume` objects. A cube file's second
   comment line says what the grid holds ("Electron density from Total SCF
   Density"), so a file with one grid is named after it; the first line is
   the job title. Vibrational modes from frequency calculations become
   :class:`~miles.core.vibrations.VibrationalMode` objects; when Open
   Babel's list is the same modes several times over (Q-Chem prints the
   harmonic analysis again in an anharmonic calculation), one copy is kept.
6. Open Babel reads the geometry and the vibrations of Molden files, but not
   their orbitals, so the file is read a second time by
   :func:`~miles.files.molden.read_molden_orbitals` for the basis set
   and orbitals (the basis functions themselves are
   :mod:`miles.core.orbitals`). When the ``[Atoms]`` section matches the structure Open
   Babel read, the molecule takes its coordinates, so the orbitals line up
   exactly. When it does not, or the orbitals cannot be read, the file
   still opens, without orbitals, and ``JobInfo.notes`` says why. Open Babel
   also reads infrared intensities (``[INT]``) only when they come before
   the normal modes; files that write them after, as many programs do, have
   them read separately.

7. For PDB files, :func:`~miles.files.openbabel.add_pdb_columns` reads
   what the ATOM and HETATM records hold beyond what Open Babel's bindings
   give back: occupancy and B-factor (the properties ``occupancy`` and
   ``b_factor``), the alternate location (``altloc``) and a formal charge
   written in columns 79-80. Open Babel makes one atom of each record, in
   file order, so the records are matched to the atoms by residue key and
   atom name, taking repeated keys (atoms in alternate locations) in order.

Writing goes the other way through :func:`~miles.core.chem.to_openbabel`,
which carries formal charges (as stated, or worked out from the bonds; see
:func:`~miles.core.charges.formal_charges`), isotopes, residues with their
insertion codes, and tells Open Babel the chains are already perceived so
its writers keep the chains and residue numbers.
The format is again chosen from the extension. PDB files are the
exception: Open Babel's PDB writer leaves out alternate locations,
occupancies and B-factors, so :mod:`miles.files.pdb` writes the records
from the columns itself, with TER after each polymer chain and CONECT
records for the bonds residue templates do not cover (ligands, disulfides
and other links between residues), a double bond written by repeating the
partner. It refuses structures of more than 99,999 atoms, which need mmCIF.
:func:`~miles.files.openbabel.write_text` writes to a string; the
clipboard uses it to offer XYZ text to other programs.

Gaussian input
**************

Open Babel writes Gaussian input but cannot read it, so
:mod:`miles.files.gaussian` does. It skips the Link 0 lines, takes the
route section, the title, and the charge and multiplicity, then reads the
molecule: Cartesian lines (with an optional freeze code, and the extra
fields of ONIOM inputs ignored) or Z-matrix lines whose references are atom
numbers or labels and whose values are numbers or variables defined after a
blank line. Z-matrix atoms are placed with
:func:`~miles.core.geometry.place_atom`. Dummy atoms (``X``) and ghost
atoms (``Bq``, ``C-Bq``) place other atoms but are left out; ``units=bohr``
and radians are understood, and only the first job of a ``--Link1--`` file
is read. The structure then goes through Open Babel as XYZ text, so bonds
are perceived as for any XYZ file. An input that takes its geometry from a
checkpoint file (``geom=check``) has no molecule, and the error says so.
Pasted text that starts with ``%`` or ``#`` goes through the same reader.

Pasted text
***********

:func:`~miles.files.formats.parse_text` turns clipboard text into a
structure. :func:`~miles.files.formats.sniff_format` guesses the format
from the content, in this order:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - The text
     - is read as
   * - starts with an atom count on its own line
     - XYZ
   * - contains ``@<TRIPOS>``
     - MOL2
   * - contains ``V2000`` or ``V3000``
     - MOL/SDF
   * - has ``ATOM`` or ``HETATM`` records
     - PDB
   * - starts with ``data_``
     - mmCIF
   * - starts with ``<`` and mentions ``molecule``
     - CML
   * - starts with ``%`` or ``#``
     - Gaussian input
   * - has lines of a symbol followed by numbers
     - plain coordinates
   * - is a single word on one line
     - SMILES

Plain coordinates are lines such as ``C 0.0 0.0 0.0``. Atom labels may carry
numbers (``C12``), extra columns are ignored and the last three numbers of
each line are taken as the position, so tables copied from papers or
spreadsheets usually work.

Session files
*************

A session (``.miles``) stores everything MILES knows about a structure:
the atoms with all their fields, including representation and residue data,
the bonds, pinned measurements, trajectory frames and their energies, grids (including
orbitals computed on a grid), vibrational modes, molecular orbitals and the
view, including the isosurface shown (its grid, isovalue, see-through and
colouring grid, under ``view["volume"]``). The canvas reads that entry when
it is created, before the grid panel first looks at it, so a session with
orbitals opens on its own orbital rather than on a newly computed HOMO.
:mod:`miles.files.session` writes it as gzip-compressed JSON. The
coordinates, frames, mode displacements and orbital coefficients are
stored as base64-encoded float64 arrays, so a session gives back exactly
what was saved; grids are float32, as they are in memory. Every column of
the atom table is saved by name, so a column added to
:data:`~miles.core.molecule.ATOM_COLUMNS` is saved without changing the
session code, and columns an older version does not know are skipped when
it reads the file. The selection and the named selections are saved with
the view; results that no longer apply to the structure are left out. The
file starts with a format name and version, so later versions can read
older files (version 1 stored float32 coordinates).

Saving writes a temporary file next to the session and renames it into
place only once it is complete, so a save that fails (a full disk, or a
value in plugin state that JSON cannot hold) leaves the previous file as
it was; NumPy numbers in plugin state are saved as plain numbers. A damaged
file raises :class:`~miles.files.common.FileFormatError` rather than
another error.

Cube files
**********

Grids are read through Open Babel like any other file, and written by
:mod:`miles.files.cube`, which needs no Open Babel:
:func:`~miles.files.cube.write_cube` follows the layout of Gaussian's
``cubegen`` (two comment lines, the origin and step vectors in bohr, which
need not be orthogonal, one line per atom, then six values per line with
the last axis fastest). The second comment line holds the grid's name,
which is what the reader takes as the name, so a grid keeps its name
through a write and a read. The tests check this round trip on orthogonal
and skewed grids and that Open Babel reads the files.

Downloads
*********

``builtin_plugins/structure_lookup/remote.py``, part of the Structure
Lookup plugin, fetches structures from two services.

* :func:`~miles.builtin_plugins.structure_lookup.remote.fetch_pdb` downloads an entry from the RCSB
  Protein Data Bank by its four-character code, trying the PDB format first
  and mmCIF for entries too large for it. Only the first model is kept.
* :func:`~miles.builtin_plugins.structure_lookup.remote.fetch_compound` asks PubChem for a
  compound by name and downloads its 3D structure, falling back to its SMILES
  string when there is no 3D record.

:func:`~miles.builtin_plugins.structure_lookup.remote.resolve` is what the search box calls: valid
SMILES is built locally with Open Babel, anything else is treated as a name
and looked up on PubChem. Both functions use only the standard library for
the network requests. The interface runs them on a worker thread (see
:doc:`interface`) so a slow connection never freezes the window.

Checked against real outputs
****************************

The readers were run over the test data of cclib (BSD licence) and IOData
(LGPL), about 1,200 files from 20 programs, each file in its own process
with a time limit, and compared with what cclib reads from the same files:
atoms and their order, the last geometry, the number of steps, frequencies,
infrared intensities, orbital energies, charges, charge and multiplicity.
The orbitals of every Molden file were checked to be orthonormal with the
exact overlap of their basis. The files themselves are not part of the
project; the cases they showed are tested with short hand-written inputs in
``tests/test_outputs.py``.

What was wrong and is fixed:

* GAMESS and MOPAC inputs hung the reader; GAMESS-UK outputs were read as
  GAMESS (no structure), and its frequency outputs, once read as GAMESS-UK,
  aborted the program.
* An ORCA scan showed its first geometry, a Gaussian unrelaxed scan lost its
  last point and with it every energy, Firefly optimisations had no steps,
  and Q-Chem ones only the last geometry.
* ORCA 5 and 6 infrared intensities came from the wrong column, and ORCA 6
  frequency outputs with symmetry had no modes at all.
* Formatted checkpoint files had their dipole in atomic units shown as
  debye, and no charges or modes.
* Total charge and multiplicity were 0 and 1 for ions and radicals from
  most programs, and a TDDFT triplet job took the multiplicity of its
  excited states.
* Molden orbitals from CFOUR, Turbomole and Psi4 up to 1.3.2 were not
  normalised, because these programs normalise cartesian d, f and g
  functions in their own way (:func:`~miles.core.orbitals.overlap_matrix`
  now picks the reading that makes every orbital normalised), and Psi4's
  capital atom labels became dummy atoms.
* Gaussian input files, listed as readable, could not be read.

What Open Babel cannot read, so MILES cannot either: outputs of Psi4,
Turbomole, DALTON, Molcas, xtb, Serenity and CP2K (the error names the
program), most ADF outputs, and ORCA and NWChem input files. From the
outputs it reads, it misses the steps of ADF, Jaguar and GAMESS-UK
optimisations (the optimised geometry is shown) and of Gaussian and NWChem
dynamics runs (a single geometry is shown), the normal modes of ADF, Jaguar
and GAMESS-UK, the Mulliken charges of NWChem, Jaguar and GAMESS-UK, and the
dipole moment of ORCA, Q-Chem and NWChem. Where a file reports several
dipole moments (MP2 or CIS after Hartree-Fock, or a solvated one), the first
is the one read. ORCA before version 5 prints no infrared intensities in
km/mol, so its modes have none. Molecular orbitals are read from Molden
files only.
