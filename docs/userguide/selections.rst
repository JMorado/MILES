.. _ug-selections:

Selecting with Expressions
##########################

Clicking selects the atoms you can point at. An expression selects the ones
you can describe: every alpha carbon, the water near a ligand, the residues
that line a binding pocket. MILES reads expressions in a small language
close to those of PyMOL and VMD::

   protein and name CA
   water within 3.5 of ligand
   byres within 5 of resname HEM
   sidechain and chain A and resid 40-60
   b_factor > 60 and not water
   smarts "[CX3](=O)[OX2H1]"

The same language works in the **Select** dialog, in the Python console and
in scripts and plugins. It is quick on large structures: picking the
residues within 5 angstrom of the ligands of a structure of 100,000 atoms
takes a few hundredths of a second.

The Select Dialog
*****************

**Edit > Select...** (:kbd:`Cmd F`) opens a field for an expression. While
you type, the atoms it picks are selected in the window and the line below
the field counts them ("212 atoms in 18 residues"), or says what is wrong
and where ("Expected 'of' after 'within 5' (at character 10)"). An unknown
word gets a suggestion: ``protien`` asks "Did you mean 'protein'?".

* **Select** keeps those atoms selected; **Cancel** brings the old
  selection back.
* With **Add to the current selection** checked, the atoms are added to
  what was selected before.
* **Examples** fills in common expressions: the protein, the ligands, the
  pocket around them, nearby water, the backbone, charged residues, polar
  hydrogens, metal ions, and the functional groups of
  :ref:`ug-select-pattern`. Named selections of the document (see
  :ref:`ug-named-selections`) come first.
* The :kbd:`Up` and :kbd:`Down` keys go through the expressions you selected
  before, as in a terminal.

The selection is an ordinary one: show it in another style, hide everything
else with **View > Show Only Selection**, measure it, or copy it.

Words
*****

Keywords can be typed in any case (``Protein``, ``NOT``). They fall into
four kinds: groups, values, numbers and distances.

Groups
======

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Word
     - Atoms
   * - ``all``, ``none``
     - Every atom, no atom
   * - ``protein``
     - Atoms of amino acid residues, including the usual protonation
       variants (HID, HIE, HIP, CYX, ASH ...)
   * - ``nucleic``
     - Atoms of DNA and RNA nucleotides
   * - ``water``
     - Water molecules (HOH, WAT, TIP3, SOL ...)
   * - ``ion``
     - Metal and halide atoms that stand alone: without bonds, or alone in
       their residue (ZN, NA, CL ...)
   * - ``ligand``
     - Everything that is not protein, nucleic acid, water or ion: bound
       ligands, cofactors such as heme, sugars, and in a file without
       residue data the molecules themselves
   * - ``hetero``
     - Atoms read from HETATM records
   * - ``hydrogen``, ``heavy``
     - Hydrogen atoms, and all the others
   * - ``backbone``
     - N, CA, C, O and OXT of amino acids, and the sugar-phosphate atoms of
       nucleotides
   * - ``sidechain``
     - Atoms of amino acids that are not backbone, nor hydrogens on it
   * - ``polar``, ``apolar``
     - Nitrogen, oxygen and the hydrogens bonded to them; and the rest
   * - ``aromatic``
     - Atoms of aromatic bonds
   * - ``metal``
     - Atoms of metallic elements, bonded or not
   * - ``helix``, ``sheet`` (or ``strand``), ``loop`` (or ``coil``)
     - Residues of each kind of secondary structure, as the cartoon shows it
   * - ``selected`` (or ``sele``)
     - The atoms selected when the expression is worked out
   * - ``visible``, ``hidden``
     - Atoms that are drawn, and those hidden by focus mode or by the
       **Hidden** representation

Values
======

These words are followed by one value or more; an atom matches if it has
any of them.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Word
     - Example
   * - ``name``
     - ``name CA CB``, ``name C*`` (``*`` and ``?`` are wildcards)
   * - ``resname`` (or ``resn``)
     - ``resname HIS HID HIE``
   * - ``resid`` (or ``resi``)
     - ``resid 10-20 42``, ``resid 10 to 20``
   * - ``chain``
     - ``chain A B``
   * - ``element`` (or ``elem``)
     - ``element N O``, ``element Fe``, ``element 26``
   * - ``index``
     - ``index 1-50``: atom numbers, counted from 1 as in the window

Atom and residue names also match the capitals of what you type, so
``resname his`` finds HIS. Chain identifiers are matched exactly, because
large structures use both ``A`` and ``a``. A value in quotes is matched
literally, without wildcards: ``name "O5*"`` finds atoms named O5* in old
nucleic acid files. Atom names with a prime need no quotes: ``name O5'``.

Numbers
=======

A name, a comparison (``<``, ``<=``, ``>``, ``>=``, ``==`` or ``=``,
``!=``) and a number pick atoms by a value:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Name
     - Value
   * - ``x``, ``y``, ``z``
     - Coordinates in angstrom: ``z > 10``
   * - ``resid``, ``index``
     - Residue number and atom number: ``resid >= 100``
   * - ``mass``
     - Atomic mass
   * - ``charge``
     - Partial charge of the model chosen in **View > Charges**:
       ``charge < -0.4``
   * - ``b_factor`` (or ``beta``), ``occupancy`` (or ``occ``)
     - From PDB files: ``b_factor > 60``
   * - any other number column
     - Values that a reader or a plugin added to the atoms

Distances and structure
=======================

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Words
     - Atoms
   * - ``within 5 of S``
     - Atoms at most 5 angstrom from any atom of S, S included
   * - ``around 5 of S``
     - The same without S itself
   * - ``byres S``
     - Whole residues that have an atom in S (a whole molecule for atoms
       without residue data)
   * - ``bychain S``
     - Whole chains that have an atom in S
   * - ``bonded to S``
     - Atoms bonded to an atom of S, not counting S itself
   * - ``smarts "..."``
     - Atoms of every match of a SMARTS pattern (see :ref:`ug-select-pattern`)

Combining
*********

``and``, ``or`` and ``not`` combine words, with parentheses where needed.
``&``, ``|`` and ``!`` mean the same. ``and`` binds more tightly than
``or``: ``resname LIG or water and hydrogen`` is the ligand and the water
hydrogens.

``not``, ``within ... of``, ``around ... of``, ``byres``, ``bychain`` and
``bonded to`` apply to the word or bracket right after them, so
``byres within 5 of ligand and chain A`` is
``(byres (within 5 of ligand)) and chain A``. Put brackets around a longer
selection: ``within 5 of (resname LIG and chain A)``.

As in PyMOL, ``water within 3.5 of ligand`` is read as
``water and within 3.5 of ligand``.

Examples
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Expression
     - Selects
   * - ``byres within 5 of ligand``
     - The binding pocket: residues with an atom within 5 angstrom of a
       ligand (the ligands included)
   * - ``protein and byres around 5 of resname ATP``
     - Protein residues that line the ATP site
   * - ``water within 3.5 of ligand``
     - Waters close enough to bridge the ligand and the protein
   * - ``sidechain and byres within 4 of ligand``
     - Side chains in contact with ligands
   * - ``protein and polar and within 3.5 of (ligand and polar)``
     - Protein atoms that could form hydrogen bonds with a ligand
   * - ``bonded to element Zn``
     - The atoms that coordinate a zinc ion (when the file has the bonds)
   * - ``ion within 3 of (protein and element S)``
     - Ions close to cysteine or methionine sulfurs
   * - ``chain A and resid 25-40 and backbone``
     - A stretch of backbone
   * - ``protein and name CA and helix``
     - Alpha carbons of helices
   * - ``b_factor > 60 and not water``
     - Flexible or poorly resolved parts
   * - ``not (protein or water or ion)``
     - The same as ``ligand``
   * - ``hydrogen and bonded to element N O``
     - Polar hydrogens
   * - ``byres selected``
     - Grows the current selection to whole residues
   * - ``bychain resname HEM``
     - The chains that bind a heme

.. _ug-named-selections:

Named Selections
****************

A name can stand for an expression. In the console::

   >>> define("pocket", "byres within 5 of ligand")
   >>> select("pocket and not hydrogen")

The name is then a word like the others, in the Select dialog too, where
named selections are listed first under **Examples**. Names are plain words
of letters, digits and underscores that are not keywords. A name keeps its
expression, not a list of atoms, so it follows edits: after mutating a
residue or moving the ligand, ``pocket`` is the new pocket. Session files
keep the names of each document (:ref:`ug-sessions`). ``define("pocket",
None)`` removes one.

.. _ug-select-pattern:

SMARTS Patterns
***************

``smarts "..."`` selects every group of atoms that matches a SMARTS pattern,
the substructure language of cheminformatics: ``[CX3]=[OX1]`` is a carbonyl
group, ``[OX2H]`` an OH group, ``a`` any aromatic atom. **Examples** in the
Select dialog offers common ones: carbonyl groups, carboxylic acids and
carboxylates, amides, hydroxyl groups of alcohols and phenols, amine
nitrogens, aromatic ring atoms, halogens, hydrogen-bond donors and
acceptors, and sulfur. Patterns combine with the rest of the language:
``ligand and smarts "[OX2H]"``.

Hydrogens are counted the way the structure suggests. A molecule that has
hydrogens is taken as complete, so an oxygen with one bond and no hydrogen
is a charged oxygen, as in a carboxylate, and ``[O-]`` finds it. A protein
read from a PDB file usually has no hydrogens at all; then each atom is
given its usual number, so ``[OX2H]`` still finds the OH groups of serine,
threonine and tyrosine. The definitions of donors and acceptors are simple
ones based on elements and bonding, meant for picking atoms rather than for
counting them precisely.

From Python
***********

In the console, ``select``, ``count`` and ``define`` take expressions
(see :ref:`ug-python`):

.. code-block:: pycon

   >>> select("protein and name CA")
   [1, 9, 18, ...]
   >>> count("water within 3.5 of ligand")
   7

Scripts and plugins use the methods of the document,
:meth:`~miles.api.DocumentView.select`, :meth:`~miles.api.DocumentView.find`
(the atoms, without selecting them), :meth:`~miles.api.DocumentView.count`
and :meth:`~miles.api.DocumentView.define`. A wrong expression raises an
error that says what is wrong and where.
