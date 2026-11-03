The miles library
#################

What scripts and plugins compute with: ``import miles``. Together with :mod:`miles.plugin` this is the public, stable surface of MILES (see :ref:`dg-stability`).

.. automodule:: miles
   :no-members:

Opening the window
******************

Start MILES from IPython, Jupyter or a script.

.. autofunction:: miles.launch

Structures and files
********************

Molecules, their atoms, and reading and writing structure files.

.. autoclass:: miles.Molecule
   :no-index:

.. autoclass:: miles.Atom
   :no-index:

.. data:: miles.AtomLike
   :no-index:

   Either kind of atom record: an :class:`Atom` or a view of one atom of a molecule.

.. autoclass:: miles.Representation
   :no-index:

.. autoclass:: miles.Geometry
   :no-index:

.. autoclass:: miles.JobInfo
   :no-index:

.. autoexception:: miles.FileFormatError
   :no-index:

.. autofunction:: miles.read_file
   :no-index:

.. autofunction:: miles.write_file
   :no-index:

.. autofunction:: miles.parse_text
   :no-index:

.. autofunction:: miles.read_text
   :no-index:

.. autofunction:: miles.write_cube
   :no-index:

Selections
**********

The selection language and its helpers.

.. autofunction:: miles.find
   :no-index:

.. autoexception:: miles.SelectionError
   :no-index:

.. autofunction:: miles.residue_ids
   :no-index:

Geometry, neighbours and trajectories
*************************************

Measuring and setting geometry, finding neighbours, and aligning frames.

.. autofunction:: miles.distance
   :no-index:

.. autofunction:: miles.angle
   :no-index:

.. autofunction:: miles.dihedral
   :no-index:

.. autofunction:: miles.kabsch
   :no-index:

.. autofunction:: miles.rotation_matrix
   :no-index:

.. autofunction:: miles.unit_vector
   :no-index:

.. autofunction:: miles.place_atom
   :no-index:

.. autofunction:: miles.set_bond_length
   :no-index:

.. autofunction:: miles.set_angle
   :no-index:

.. autofunction:: miles.set_dihedral
   :no-index:

.. autofunction:: miles.pairs_within
   :no-index:

.. autofunction:: miles.cross_pairs
   :no-index:

.. autofunction:: miles.near_mask
   :no-index:

.. autofunction:: miles.align_frames
   :no-index:

.. autofunction:: miles.rmsd_atoms
   :no-index:

.. autofunction:: miles.rmsd_series
   :no-index:

Elements
********

Element data.

.. autofunction:: miles.element
   :no-index:

.. autofunction:: miles.normalize_symbol
   :no-index:

Open Babel, SMILES and force fields
***********************************

Converting to and from Open Babel, SMILES, minimisation and conformers.

.. autofunction:: miles.to_openbabel
   :no-index:

.. autofunction:: miles.from_openbabel
   :no-index:

.. autofunction:: miles.openbabel_coordinates
   :no-index:

.. autofunction:: miles.silence_openbabel
   :no-index:

.. autofunction:: miles.from_smiles
   :no-index:

.. autofunction:: miles.to_smiles
   :no-index:

.. autofunction:: miles.is_smiles
   :no-index:

.. autoexception:: miles.ChemError
   :no-index:

.. data:: miles.FORCE_FIELDS
   :no-index:

   Force fields Open Babel offers for :func:`minimize` and :func:`conformers`, in order of preference.

.. autofunction:: miles.minimize
   :no-index:

.. autofunction:: miles.conformers
   :no-index:

.. autofunction:: miles.clean_geometry
   :no-index:

.. autofunction:: miles.clean_selection
   :no-index:

Building and perception
***********************

Adding hydrogens, residue templates and working out bonds and atom types.

.. autofunction:: miles.add_hydrogens
   :no-index:

.. autofunction:: miles.clash_score
   :no-index:

.. autofunction:: miles.amino_acid_codes
   :no-index:

.. autofunction:: miles.amino_acid_template
   :no-index:

.. autofunction:: miles.rebond
   :no-index:

.. autofunction:: miles.infer_atom_types
   :no-index:

.. autofunction:: miles.infer_geometry
   :no-index:

Proteins and nucleic acids
**************************

Recognising residues, polymers and water.

.. data:: miles.AMINO_ACIDS
   :no-index:

   Residue names of the standard amino acids, including protonation variants such as HID.

.. data:: miles.NUCLEOTIDES
   :no-index:

   Residue names of DNA and RNA nucleotides.

.. data:: miles.ONE_LETTER
   :no-index:

   One-letter codes of amino acid and nucleotide residue names, including protonation variants.

.. autofunction:: miles.residues
   :no-index:

.. autofunction:: miles.is_polymer
   :no-index:

.. autofunction:: miles.has_polymer
   :no-index:

.. autofunction:: miles.is_water
   :no-index:

.. autofunction:: miles.water_mask
   :no-index:

Spectra, quantum chemistry and charges
**************************************

Vibrations, orbitals, grids and partial charges.

.. autoclass:: miles.VibrationalMode
   :no-index:

.. autoclass:: miles.OrbitalSet
   :no-index:

.. autoclass:: miles.Spin
   :no-index:

.. autoclass:: miles.Volume
   :no-index:

.. autofunction:: miles.mapped_isosurface
   :no-index:

.. autoclass:: miles.ChargeModel
   :no-index:

.. autofunction:: miles.dipole_moment
   :no-index:

.. autofunction:: miles.formal_charges
   :no-index:

Units
*****

Factors that convert the units of files into the units of MILES: angstrom, kcal/mol, elementary charges and debye.

.. data:: miles.BOHR
   :no-index:

   One bohr in angstrom.

.. data:: miles.HARTREE
   :no-index:

   One hartree in kcal/mol.

.. data:: miles.HARTREE_EV
   :no-index:

   One hartree in electronvolts.

.. data:: miles.ELECTRONVOLT
   :no-index:

   One electronvolt in kcal/mol.

.. data:: miles.DEBYE_PER_AU
   :no-index:

   One atomic unit of dipole moment (e times bohr) in debye.
