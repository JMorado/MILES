.. _ug-files:

Reading and Writing Files
#########################

MILES reads and writes structures with `Open Babel
<https://openbabel.org>`_. The format is chosen from the file extension, so
``ligand.sdf`` is read as an SD file and ``model.pdb`` written as PDB.

Opening Files
*************

There are three ways to open files:

* **File > Open** (:kbd:`Cmd O`) lets you choose one or more files; each
  opens in its own tab.
* **File > Open Recent** lists the last ten files you opened or saved, most
  recent first. A file that has been moved or deleted is taken off the
  list when you choose it; **Clear Menu** empties the list.
* **Drag** files from the Finder or a file manager onto the window.
* **Give them on the command line**: ``miles a.xyz b.pdb``.

The open dialog lists every format Open Babel can read, about 150 of them,
starting in the folder you used last. The most common are:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Extension
     - Format
   * - ``.xyz``
     - XYZ coordinates
   * - ``.mol``, ``.sdf``
     - MDL molfile and SD file
   * - ``.mol2``
     - Tripos MOL2
   * - ``.pdb``, ``.cif``
     - Protein Data Bank and mmCIF
   * - ``.cml``
     - Chemical Markup Language
   * - ``.gjf``, ``.com``
     - Gaussian input (see below)
   * - ``.log``, ``.out``
     - Output of quantum chemistry programs (see below)
   * - ``.fchk``
     - Gaussian formatted checkpoint file, also written by Q-Chem
   * - ``.molden``, ``.cube``
     - Molecular orbitals and grids (see below)
   * - ``.miles``
     - MILES session (see :ref:`ug-sessions`)

What is read
============

Bonds
   Formats that store bonds, such as MOL, SDF and MOL2, keep their bonds and
   bond orders. For formats with coordinates only, such as XYZ, bonds and
   bond orders are worked out from the distances between atoms.

Residues
   PDB, mmCIF and the other formats that carry residue information keep it:
   residue names and numbers, insertion codes (residue ``52A``), chains,
   atom names and the heteroatom flag. PDB files also give each atom's
   alternate location, occupancy and B-factor; atoms in two alternate
   locations are both kept. This is what the protein features of
   :ref:`ug-proteins` use.

Charges and isotopes
   SDF, MOL, CML and SMILES state formal charges, and PDB files may give
   them in columns 79-80; MILES keeps them. Where no charge is given it is
   worked out from the bonds when needed: an ammonium nitrogen with four
   bonds is +1. Isotopes, such as the carbon-13 of ``C[13CH3]``, are kept
   too.

Several structures
   When a file holds several structures, for example a multi-model PDB file
   or an XYZ trajectory, the first one is shown. Further structures of the
   same molecule (the same elements in the same order, and for SDF, MOL2
   and CML also the same bonds) become frames of the document, which you can
   step through and play with the timeline (:ref:`ug-timeline`). Other
   molecules, such as the rest of an SDF library, are not read: a message
   says how many more structures the file holds. XYZ files written by xtb,
   ORCA, CREST or ASE (extended XYZ) carry an energy on each comment line;
   it is read too and plotted on the timeline.

Quantum chemistry output
   Files ending in ``.log``, ``.out`` or ``.output`` (also gzipped, as
   ``.out.gz``) are recognised by their content as output of Gaussian,
   ORCA, Q-Chem, NWChem, GAMESS or Firefly, GAMESS-UK, Jaguar, Molpro or
   MOPAC; ADF output is read from ``.adfout`` files. MILES shows the
   last geometry in the file, which for an optimisation is the optimised
   structure; the earlier steps, and the points of a scan, are frames of
   the document, with the energy of each step where the output has it
   (Gaussian, ORCA, Q-Chem, NWChem, GAMESS and Firefly), which the timeline
   plots. Vibrational modes from frequency calculations are read as well
   (:ref:`ug-vibrations`), and so are Mulliken charges and the dipole
   moment (:ref:`ug-charges`), and the total charge and spin multiplicity
   of the calculation.

   Output of Psi4, Turbomole, DALTON, Molcas, xtb, Serenity and CP2K
   cannot be read; the message says which program wrote the file. Save the
   structure from that program as an XYZ or Molden file instead, which
   most of them can write.

Formatted checkpoint files
   ``.fchk`` files from Gaussian and Q-Chem give the structure, the
   Mulliken charges, the dipole moment and the normal modes of a frequency
   calculation. Q-Chem writes one for every step of an optimisation into
   the same file; the steps become frames. Their molecular orbitals are not
   read; open a Molden file for those.

Gaussian input
   ``.gjf`` and ``.com`` files give the structure, the title, and the
   charge and multiplicity. The molecule may be given by Cartesian
   coordinates or as a Z-matrix, with numbers or with variables defined
   after it; dummy atoms (``X``) and ghost atoms (``Bq``) are left out. An
   input that takes its geometry from a checkpoint file (``geom=check``)
   holds no structure. Pasting the text of an input works too.

Volumetric data and orbitals
   Cube files keep their grids, which are shown as isosurfaces, and Molden
   files keep their molecular orbitals, which are computed on a grid when
   chosen (:ref:`ug-volumes`). Programs write the orbitals of Molden files
   in slightly different ways (ORCA, CFOUR, Turbomole and older Psi4
   versions each have their own); MILES works out which from the
   orbitals themselves, which must come out normalised. If they do not in
   any known way, the file opens without its orbitals rather than with
   wrong ones. **Add Grid...** in the grid panel adds the
   grids of another file of the same structure, such as an electrostatic
   potential to colour a density surface.

If a file cannot be read, a dialog explains why and nothing else changes.

Saving Files
************

**File > Save** (:kbd:`Cmd S`) writes the current document back to its file,
in the file's format. A new document has no file yet, so Save asks for a
name; so does **File > Save As** (:kbd:`Shift Cmd S`), which always asks.

The save dialog offers the common formats first: XYZ, SDF, MOL, MOL2, PDB,
CIF, CML and MILES sessions. Choose **Any other format** and type an
extension to use any other format Open Babel can write, such as ``.gjf`` for
a Gaussian input file. A name without an extension is saved as XYZ.

Files ending in ``.log`` or ``.out`` are program output, which MILES
never overwrites. Saving a document opened from such a file asks for a new
name.

Grids are saved from the grid panel rather than with Save As:
**Save Grid...** writes the grid shown as a cube file (:ref:`ug-volumes`).

What a file keeps depends on its format: XYZ stores elements and
coordinates only; MOL and SDF add bonds, bond orders, formal charges and
isotopes; and PDB and mmCIF keep residue information. PDB files are
written by MILES itself, so they keep everything a PDB file can hold:
chains, residue numbers and insertion codes as they were read, alternate
locations, occupancies, B-factors, formal charges, and the bonds of ligands
and between residues (with their orders) in CONECT records. Atoms without
residue data are written as residue ``UNL`` with names such as ``C1``. A
PDB file holds at most 99,999 atoms; save larger structures as mmCIF
(``.cif``). A document with several frames is saved with the frame on
screen. None of the standard formats remember how each atom is drawn, the
frames or the view. For that, use a session.

.. _ug-sessions:

Sessions
********

A session file (``.miles``) stores a document exactly as MILES holds it:

* elements, coordinates, bonds and bond orders;
* the bonding type of every atom, so building continues as before;
* residue and chain information, secondary structure and atom names;
* the representation of every atom (cartoon, balls and sticks, hidden and so
  on) and which atoms are hidden in focus mode;
* pinned measurements;
* all frames, with their labels and energies;
* grids and molecular orbitals, including orbitals already computed, and
  which one is shown, at what isovalue, see-through or not and coloured by
  which other grid;
* the view: the camera position, orientation and zoom, and the display
  options of the **View** menu (default style, labels, hydrogens, depth
  fade, cartoon colours, stereo labels and hydrogen bonds);
* the title and the total charge.

Save a session with **File > Save As** and the **Session with styles** file
type, and open it like any other file: the structure comes back looking
exactly as it did. Surfaces and the interaction lines of a pocket view are
not stored; show them again after opening.

Session files are compressed JSON, so they stay small even for large
proteins.

.. _ug-autosave:

Autosave and Recovery
*********************

Every minute, MILES keeps a copy of each document with unsaved
changes, in the session format, so the frames, the view and the atoms
hidden in focus mode are kept too. Only documents changed since their last
copy are written, and the writing happens in the background, so even a
large protein does not interrupt your work. The copies are kept in the
application data folder, apart from your own files, and are deleted when
you save or close the document and when you quit MILES normally.

If MILES stops without quitting normally (a crash, a forced quit, or
the computer shutting down), it offers the documents with unsaved changes
the next time it starts, with the time each copy was made:

**Restore**
   Open them in tabs. Each comes back with its file, if it had one, and
   marked as changed: nothing is written to your files until you save, and
   Save writes to the file the document came from.

**Discard**
   Delete the copies. Your files stay as they were last saved.

Closing the window without choosing keeps the copies, and they are offered
again at the next start. While several MILES windows run at the same
time, none of them offers the copies of another that is still open.

The interval, or switching autosave off, is chosen in Settings
(:ref:`ug-settings`). Autosave does not replace saving: it only protects
the changes made since the last save.
