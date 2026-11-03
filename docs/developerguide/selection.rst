.. _dg-selection:

The Selection Language
######################

Expressions such as ``byres within 5 of ligand`` are read and worked out by
:mod:`miles.core.selection`, a module without Qt that every part of MILES
shares: the Select dialog (the ``select_atoms`` built-in plugin), the
console helpers in :mod:`miles.scripting`, the API
(:meth:`~miles.api.DocumentView.select`, ``find``, ``count``, ``define``)
and the document core (:meth:`miles.core.document.Document.find`). The user
guide (:ref:`ug-selections`) describes the language; this chapter describes
how it works and how to use and extend it.

From Text to Atoms
******************

Working out an expression has two steps.

:func:`~miles.core.selection.parse` turns the text into a tree of small
immutable nodes (:class:`~miles.core.selection.Flag`,
:class:`~miles.core.selection.Values`,
:class:`~miles.core.selection.Compare`,
:class:`~miles.core.selection.Within`,
:class:`~miles.core.selection.Expand`,
:class:`~miles.core.selection.Smarts`,
:class:`~miles.core.selection.Named`, and
:class:`~miles.core.selection.Not`, :class:`~miles.core.selection.Both`,
:class:`~miles.core.selection.Either` for the logic). A regular expression
splits the text into words, quoted strings, brackets and operators, and a
recursive-descent parser reads them with this grammar, loosest first::

   either   := both ("or" both)*
   both     := prefixed (("and" prefixed) | within-term)*
   prefixed := "not" prefixed
             | ("within" | "around") NUMBER "of" prefixed
             | ("byres" | "bychain") prefixed
             | "bonded" "to" prefixed
             | primary
   primary  := "(" either ")" | GROUP | FIELD value+ | NAME OP NUMBER
             | "smarts" PATTERN | NAME

A ``within`` or ``around`` term right after another term is joined to it
with ``and``, which is how ``water within 3.5 of ligand`` reads. Values
after a field word run until a keyword, a bracket or an operator. Trees are
cached by text (:func:`functools.lru_cache`), so the Select dialog, which
works the field out after every pause in typing, parses each text once.

:func:`~miles.core.selection.evaluate` walks the tree and returns one
boolean mask over the atoms; :func:`~miles.core.selection.find` returns
the indices of the true ones. Errors of both steps are
:class:`~miles.core.selection.SelectionError`, a :class:`ValueError` with
the character where the problem is (``position``) and a ``pointer()`` that
draws a caret under it. Unknown words and numbers get a suggestion from
:mod:`difflib`.

Working on Columns
******************

Every node becomes array arithmetic on the columns of the molecule (see
:doc:`molecule_model`); nothing loops over atoms in Python. A few ideas
keep it fast:

Runs of string columns
   Residue names, chains and secondary structure are object columns of
   interned strings. Atoms of a residue are stored together, so these
   columns are long runs of equal values. The evaluator compares each
   column with itself shifted by one atom to find the runs, decides each
   distinct value once in Python, and spreads the answer back over the atoms.
   A protein of 100,000 atoms has about 12,000 residue runs.

Direct comparisons for atom names
   Atom names change from one atom to the next, so runs do not help; a short
   list of exact names is compared with the whole column instead.

Group numbers
   ``byres`` and ``bychain`` need a number per residue or chain.
   :func:`~miles.core.selection.residue_ids` finds runs of equal chain,
   residue number and residue name, then gives equal keys one number, so a
   residue split by later edits (hydrogens added at the end) is still one
   residue. Atoms without residue data are grouped by molecule, using the
   connected components of the bond graph.

A grid for distances
   ``within`` first keeps the atoms inside the box around the inner atoms
   widened by the distance, then sorts both sets into cubic cells as wide as
   the distance and compares each atom with the 27 cells around it. The
   smaller set is walked and the larger one looked up in sorted cell order;
   candidate pairs are looked at in batches of about two million, so memory
   stays bounded on large structures. Next to ``and``, the search only tests
   atoms that the other side keeps: in ``water within 3.5 of ligand`` only
   water is tested.

One evaluator works out one expression and keeps what several of its parts
share (runs, group numbers, components, the masks of group words), so
``ligand`` costs little more after ``protein`` has been worked out.

Measured on copies of the haemoglobin structure 4HHB, each with its own
chains, on an Apple M-series laptop:

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * - Expression
     - 14,000 atoms
     - 100,000 atoms
     - 1,000,000 atoms
   * - ``within 5 of resname LIG`` (one ligand)
     - 1 ms
     - 4 ms
     - 36 ms
   * - ``byres within 5 of ligand`` (all hemes)
     - 5 ms
     - 34 ms
     - 390 ms
   * - ``water within 3.5 of ligand``
     - 2 ms
     - 14 ms
     - 145 ms
   * - ``protein and name CA``
     - 0.3 ms
     - 2 ms
     - 25 ms
   * - ``backbone``
     - 2 ms
     - 13 ms
     - 126 ms

Context
*******

Some words depend on more than the molecule. A
:class:`~miles.core.selection.Context` carries them: the selected atoms
(``selected``), a mask of the atoms not drawn (``visible``, ``hidden``),
named selections (each name and its expression) and extra numeric columns
by name, given as arrays or as callables that are only called when an
expression uses them. :meth:`miles.core.document.Document.selection_context`
fills it for a document: hidden atoms from focus mode and the **Hidden**
representation, and ``charge`` from :meth:`~miles.core.document.Document.charges`,
so charges are worked out only for expressions that compare them.

Named selections live on the document core as text
(``Document.selections``) and are saved in session
files under ``view["selections"]``. A name is parsed and worked out when an
expression uses it, with the evaluator's caches, and a name that refers to
itself, directly or through other names, is an error.
:meth:`~miles.core.document.Document.define` checks the name with
:func:`~miles.core.selection.check_name` and works the expression out once
before keeping it.

Using It from Plugins
*********************

Plugins usually go through the document they are given::

   def mark_pocket(app):
       doc = app.current
       pocket = doc.find("byres within 5 of ligand and not ligand")
       doc.set_representation("tube", pocket)
       app.message(f"{len(pocket)} pocket atoms")

``doc.find`` returns indices without touching the selection, ``doc.select``
selects, ``doc.count`` counts, and ``doc.define`` adds a name the user can
type too. Code that works on a bare :class:`~miles.core.molecule.Molecule`,
such as a batch script without the window, calls ``miles.find(mol,
expression)``; inside MILES, :func:`~miles.core.selection.evaluate` and
:func:`~miles.core.selection.find` take a context of their own.

A plugin that adds a numeric property column makes it comparable at once:
after ``mol.set_property("sasa", values)``, users can type ``sasa > 40``.
Values that are computed on demand belong in
the ``columns`` of a :class:`~miles.core.selection.Context` instead.

Extending the Language
**********************

To add a group word, add it to :data:`~miles.core.selection.FLAGS` (several
spellings may map to one group) and work it out in ``_Evaluator.flag`` as a
mask built from columns. Keep it vectorised: use ``where`` for string
columns and the bond arrays for anything about neighbours. A new field word
goes into :data:`~miles.core.selection.FIELDS` and ``_Evaluator.values``; a
new number into :data:`~miles.core.selection.COMPARABLE` and
``_Evaluator.numeric``. Every word in these tables is reserved
(:data:`~miles.core.selection.KEYWORDS`), so it can no longer name a
selection; prefer words that are unlikely to be residue or chain names.
Document the word in the user guide (:ref:`ug-selections`) and add cases to
``tests/test_selection.py``, which checks every word on a small hand-made
complex, the error messages and their positions, and the grid search
against brute-force distances.
