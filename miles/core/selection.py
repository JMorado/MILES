"""A language for choosing atoms: ``protein and within 5 of resname LIG``.

An expression combines words that pick atoms by what they are (``protein``,
``water``, ``ligand``, ``backbone``), by their columns (``name CA CB``,
``resid 10-20``, ``chain A``, ``element N O``), by numbers (``b_factor >
50``, ``charge < -0.3``), by distance (``within 5 of ligand``), by
structure (``byres``, ``bychain``, ``bonded to``, ``smarts "c1ccccc1"``)
and by the state of the view (``selected``, ``visible``), joined with
``and``, ``or``, ``not`` and parentheses. The full grammar is in the user
guide; :data:`KEYWORDS` lists the reserved words.

:func:`parse` turns the text into a tree and :func:`evaluate` works the
tree out on a :class:`~miles.core.molecule.Molecule` as one boolean mask.
Every step is array arithmetic on the molecule's columns: string columns
are compared run by run (atoms of one residue share their residue name, so
a protein of 100,000 atoms has about 12,000 runs to look at), and distances
use a grid of cells as wide as the cutoff. Nothing loops over atoms in
Python, so a selection of a large structure takes milliseconds.

Precedence, from loosest to tightest: ``or``, ``and``, then the prefixes
``not``, ``within D of``, ``around D of``, ``byres``, ``bychain`` and
``bonded to``, which apply to the term that follows them:
``byres within 5 of ligand and chain A`` is
``(byres (within 5 of ligand)) and chain A``.
"""

import difflib
import fnmatch
import functools
import math
import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from miles.core.elements import BY_Z, ELEMENTS, normalize_symbol
from miles.core.molecule import RESIDUE_KEY, Molecule, symbol_number
from miles.core.neighbors import near_mask
from miles.core.patterns import PatternError, match_pattern
from miles.core.protein import AMINO_ACIDS, NUCLEOTIDES, WATERS

# -- vocabulary ------------------------------------------------------------------------------

FLAGS: dict[str, str] = {
    "all": "all", "everything": "all", "none": "none",
    "protein": "protein", "nucleic": "nucleic", "water": "water", "waters": "water", "solvent": "water",
    "ion": "ion", "ions": "ion", "ligand": "ligand", "ligands": "ligand",
    "hydrogen": "hydrogen", "hydrogens": "hydrogen", "heavy": "heavy",
    "backbone": "backbone", "sidechain": "sidechain", "polar": "polar", "apolar": "apolar",
    "hetero": "hetero", "aromatic": "aromatic", "metal": "metal", "metals": "metal",
    "helix": "helix", "sheet": "sheet", "strand": "sheet", "loop": "loop", "coil": "loop",
    "selected": "selected", "sele": "selected", "visible": "visible", "hidden": "hidden",
}
"""Words that stand for a group of atoms, with the group each one means."""

FIELDS: dict[str, str] = {
    "element": "element", "elem": "element", "name": "name", "resname": "resname", "resn": "resname",
    "resid": "resid", "resi": "resid", "resseq": "resid", "chain": "chain", "index": "index",
}
"""Words followed by a list of values, with the column each one reads."""

COMPARABLE: dict[str, str] = {
    "x": "x", "y": "y", "z": "z", "resid": "resid", "resi": "resid", "resseq": "resid", "index": "index",
    "mass": "mass", "charge": "charge", "b_factor": "b_factor", "bfactor": "b_factor", "beta": "b_factor",
    "b": "b_factor", "occupancy": "occupancy", "occ": "occupancy", "q": "occupancy",
}
"""Built-in names that can be compared with a number (any numeric property column can too)."""

OPERATORS = frozenset({"and", "or", "not", "within", "around", "of", "byres", "bychain", "bonded", "to",
                       "smarts"})

KEYWORDS: frozenset[str] = frozenset(FLAGS) | frozenset(FIELDS) | OPERATORS
"""Reserved words; they cannot name a selection."""

BACKBONE_NAMES = frozenset({"N", "CA", "C", "O", "OXT"})
"""Atom names of the protein backbone."""

NUCLEIC_BACKBONE_NAMES = frozenset({"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O5'", "C5'", "C4'",
                                    "C3'", "O3'"})
"""Atom names of the sugar-phosphate backbone of nucleic acids."""

BACKBONE_HYDROGENS = frozenset({"H", "HN", "H1", "H2", "H3", "HA", "HA2", "HA3"})
"""Hydrogens on protein backbone atoms (not part of a side chain)."""

HALOGENS = frozenset({9, 17, 35, 53})
_METAL = np.array([z > 0 and BY_Z[z].is_metal for z in range(len(BY_Z))])
_MASS = np.array([e.mass or 0.0 for e in BY_Z])
_NAME_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_COMPARE = {"<": np.less, "<=": np.less_equal, ">": np.greater, ">=": np.greater_equal,
            "=": np.equal, "==": np.equal, "!=": np.not_equal}


class SelectionError(ValueError):
    """An expression that cannot be read or worked out.

    Attributes
    ----------
    message : str
        What is wrong, in words.
    text : str
        The expression.
    position : int or None
        Character (0-based) where the problem is, if it has a place.
    """

    def __init__(self, message: str, text: str = "", position: int | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.text = text
        self.position = position

    def __str__(self) -> str:
        if self.position is None:
            return self.message
        return f"{self.message} (at character {self.position + 1})"

    def pointer(self) -> str:
        """Return the expression with a caret under the problem, for monospaced output."""
        if self.position is None:
            return self.text
        return f"{self.text}\n{' ' * self.position}^"


# -- syntax tree -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Node:
    """A part of an expression; ``position`` is where it starts in the text."""

    position: int


@dataclass(frozen=True)
class Flag(Node):
    """A named group of atoms such as ``protein`` (``group`` is a value of :data:`FLAGS`)."""

    group: str


@dataclass(frozen=True)
class Values(Node):
    """``name CA CB``: a column and the values it may hold.

    ``values`` holds ``(kind, a, b)`` items: ``("text", value, quoted)`` for
    names (a quoted value is matched literally, without wildcards) and
    ``("range", low, high)`` for numbers.
    """

    field: str
    values: tuple[tuple[str, Any, Any], ...]


@dataclass(frozen=True)
class Compare(Node):
    """``b_factor > 50``: a numeric column compared with a number."""

    column: str
    operator: str
    value: float


@dataclass(frozen=True)
class Within(Node):
    """``within 5 of ...`` (``exclude`` False) or ``around 5 of ...`` (True)."""

    distance: float
    inner: Node
    exclude: bool


@dataclass(frozen=True)
class Expand(Node):
    """``byres``, ``bychain`` or ``bonded to`` (``to`` is ``"residue"``, ``"chain"`` or ``"bonded"``)."""

    to: str
    inner: Node


@dataclass(frozen=True)
class Smarts(Node):
    """``smarts "[OX2H]"``: atoms in matches of a SMARTS pattern."""

    pattern: str


@dataclass(frozen=True)
class Named(Node):
    """A word that is not a keyword: the name of a selection defined earlier."""

    name: str


@dataclass(frozen=True)
class Not(Node):
    """Atoms not in ``inner``."""

    inner: Node


@dataclass(frozen=True)
class Both(Node):
    """``left and right``."""

    left: Node
    right: Node


@dataclass(frozen=True)
class Either(Node):
    """``left or right``."""

    left: Node
    right: Node


# -- reading the text ------------------------------------------------------------------------

_TOKEN = re.compile(r"""\s*(?:
      (?P<paren>[()])
    | (?P<op><=|>=|==|!=|<|>|=|&|\||!)
    | (?P<string>"[^"]*"|'[^']*')
    | (?P<word>[^\s()<>=!&|"'][^\s()<>=!&|"]*)
    )""", re.VERBOSE)


@dataclass(frozen=True)
class _Token:
    kind: str   # "paren", "op", "string" or "word"
    text: str
    start: int

    @property
    def lower(self) -> str:
        return self.text.lower()


def _tokens(text: str) -> list[_Token]:
    found: list[_Token] = []
    pos = 0
    while pos < len(text):
        if text[pos:].strip() == "":
            break
        match = _TOKEN.match(text, pos)
        if match is None or match.end() == pos:
            start = len(text) - len(text[pos:].lstrip())
            quote = text[start]
            raise SelectionError(f"Missing closing {quote} for the quote that starts here", text, start)
        kind = match.lastgroup or "word"
        start = match.start(kind)
        found.append(_Token(kind, match.group(kind), start))
        pos = match.end()
    return found


class _Parser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.tokens = _tokens(text)
        self.i = 0

    def error(self, message: str, token: _Token | None = None) -> SelectionError:
        position = token.start if token is not None else len(self.text.rstrip())
        return SelectionError(message, self.text, position)

    def peek(self, ahead: int = 0) -> _Token | None:
        k = self.i + ahead
        return self.tokens[k] if k < len(self.tokens) else None

    def take(self) -> _Token | None:
        token = self.peek()
        if token is not None:
            self.i += 1
        return token

    def is_word(self, token: _Token | None, *words: str) -> bool:
        return token is not None and token.kind == "word" and token.lower in words

    def is_op(self, token: _Token | None, *ops: str) -> bool:
        return token is not None and token.kind == "op" and token.text in ops

    def parse(self) -> Node:
        if not self.tokens:
            raise SelectionError("Type a selection, such as protein and name CA", self.text, 0)
        node = self.either()
        token = self.peek()
        if token is not None:
            if token.kind == "paren":
                raise self.error("Unexpected ')' without a matching '('", token)
            if token.kind == "op":
                raise self.error(f"Unexpected '{token.text}'", token)
            raise self.misspelled() or self.error(f"Expected 'and' or 'or' before '{token.text}'", token)
        return node

    def misspelled(self) -> SelectionError | None:
        """Explain a word left over after an unknown one: a keyword typed wrong, as in ``chian A``."""
        previous = self.tokens[self.i - 1] if self.i else None
        if previous is None or previous.kind != "word" or previous.lower in KEYWORDS:
            return None
        hint = difflib.get_close_matches(previous.lower, [*FIELDS, *FLAGS, *OPERATORS], n=1)
        if not hint:
            return None
        return self.error(f"Unknown word '{previous.text}'. Did you mean '{hint[0]}'?", previous)

    def either(self) -> Node:
        left = self.both()
        while self.is_word(self.peek(), "or") or self.is_op(self.peek(), "|"):
            operator = self.take()
            right = self.both() if self.peek() is not None else self.missing(operator)
            left = Either(left.position, left, right)
        return left

    def both(self) -> Node:
        left = self.prefixed()
        while True:
            token = self.peek()
            if self.is_word(token, "and") or self.is_op(token, "&"):
                self.take()
                right = self.prefixed() if self.peek() is not None else self.missing(token)
            elif self.is_word(token, "within", "around"):
                right = self.prefixed()   # "water within 3.5 of ligand" reads as "water and within ..."
            else:
                return left
            left = Both(left.position, left, right)

    def missing(self, after: _Token | None) -> Node:
        word = after.text if after is not None else ""
        raise self.error(f"Expected a selection after '{word}'")

    def prefixed(self) -> Node:
        token = self.peek()
        if token is None:
            raise self.error("Expected a selection")
        if self.is_word(token, "not") or self.is_op(token, "!"):
            self.take()
            if self.peek() is None:
                raise self.error(f"Expected a selection after '{token.text}'")
            return Not(token.start, self.prefixed())
        if self.is_word(token, "within", "around"):
            self.take()
            distance = self.number(f"Expected a distance in angstrom after '{token.text}', "
                                   f"as in {token.lower} 5 of ligand")
            if distance < 0:
                raise self.error("The distance cannot be negative", self.tokens[self.i - 1])
            if not self.is_word(self.peek(), "of"):
                raise self.error(f"Expected 'of' after '{token.text} {self.tokens[self.i - 1].text}'",
                                 self.peek())
            of = self.take()
            if self.peek() is None:
                raise self.error(f"Expected a selection after '{of.text if of else 'of'}'")
            return Within(token.start, distance, self.prefixed(), token.lower == "around")
        if self.is_word(token, "byres", "bychain"):
            self.take()
            if self.peek() is None:
                raise self.error(f"Expected a selection after '{token.text}'")
            return Expand(token.start, "residue" if token.lower == "byres" else "chain", self.prefixed())
        if self.is_word(token, "bonded"):
            self.take()
            if not self.is_word(self.peek(), "to"):
                raise self.error("Expected 'to' after 'bonded', as in bonded to resname HEM", self.peek())
            self.take()
            if self.peek() is None:
                raise self.error("Expected a selection after 'bonded to'")
            return Expand(token.start, "bonded", self.prefixed())
        return self.primary()

    def number(self, message: str) -> float:
        token = self.peek()
        if token is None or token.kind != "word":
            raise self.error(message, token)
        try:
            value = float(token.text)
        except ValueError:
            raise self.error(message, token) from None
        if not math.isfinite(value):
            raise self.error(message, token)
        self.take()
        return value

    def primary(self) -> Node:
        token = self.take()
        assert token is not None
        if token.kind == "paren":
            if token.text == ")":
                raise self.error("Unexpected ')' without a matching '('", token)
            if self.peek() is None:
                raise self.error("Expected a selection after '('")
            node = self.either()
            close = self.take()
            if close is None or close.text != ")":
                raise self.error(f"Missing ')' for the '(' at character {token.start + 1}", close)
            return node
        if token.kind == "op":
            raise self.error(f"Unexpected '{token.text}'", token)
        if token.kind == "string":
            raise self.error(f"Expected a keyword such as name or resname before {token.text}", token)
        word = token.lower
        if self.is_op(self.peek(), *_COMPARE):
            return self.comparison(token)
        if word in FLAGS:
            return Flag(token.start, FLAGS[word])
        if word in FIELDS:
            return self.values(token, FIELDS[word])
        if word == "smarts":
            pattern = self.take()
            if pattern is None or pattern.kind not in ("string", "word"):
                raise self.error("Expected a SMARTS pattern after 'smarts', as in smarts \"c1ccccc1\"",
                                 pattern)
            text = pattern.text[1:-1] if pattern.kind == "string" else pattern.text
            return Smarts(pattern.start, text)
        if word in OPERATORS:
            raise self.error(f"Expected a selection before '{token.text}'", token)
        return Named(token.start, token.text)

    def comparison(self, token: _Token) -> Compare:
        operator = self.take()
        assert operator is not None
        value = self.number(f"Expected a number after '{token.text} {operator.text}'")
        return Compare(token.start, token.text, operator.text, value)

    def values(self, keyword: _Token, column: str) -> Values:
        numeric = column in ("resid", "index")
        items: list[tuple[str, Any, Any]] = []
        while True:
            token = self.peek()
            if token is None or token.kind in ("paren", "op"):
                break
            if token.kind == "word" and token.lower in KEYWORDS:
                break
            self.take()
            if numeric:
                low, high = self.number_range(token, column)
                if self.is_word(self.peek(), "to"):
                    self.take()
                    end = self.number_range(self.peek(), column)[1]
                    self.take()
                    low, high = min(low, end), max(high, end)
                items.append(("range", low, high))
            elif column == "element":
                typed = token.text[1:-1] if token.kind == "string" else token.text
                symbol = _element(typed)
                if symbol is None:
                    raise self.error(f"Unknown element '{typed}'", token)
                items.append(("text", symbol, True))
            elif token.kind == "string":
                items.append(("text", token.text[1:-1], True))
            else:
                items.append(("text", token.text, False))
        if not items:
            examples = {"element": "element C N", "name": "name CA CB", "resname": "resname HIS",
                        "resid": "resid 10-20 42", "chain": "chain A", "index": "index 1-50"}
            raise self.error(f"Expected a value after '{keyword.text}', as in {examples[column]}",
                             self.peek())
        return Values(keyword.start, column, tuple(items))

    def number_range(self, token: _Token | None, column: str) -> tuple[int, int]:
        what = "an atom number" if column == "index" else "a residue number"
        if token is None or token.kind != "word":
            raise self.error(f"Expected {what} or a range such as 10-20", token)
        match = re.fullmatch(r"(-?\d+)(?:(?:-|:)(-?\d+))?", token.text)
        if match is None:
            raise self.error(f"'{token.text}' is not {what} or a range such as 10-20", token)
        low = int(match.group(1))
        high = int(match.group(2)) if match.group(2) is not None else low
        return (low, high) if low <= high else (high, low)


def _element(text: str) -> str | None:
    """Element symbol for a typed symbol (any case) or atomic number; None if there is no such element."""
    if text.isdigit():
        return normalize_symbol(text)
    symbol = text.capitalize()
    return symbol if symbol in ELEMENTS or symbol == "X" else None


@functools.lru_cache(maxsize=256)
def parse(text: str) -> Node:
    """Read an expression into a syntax tree.

    Parameters
    ----------
    text : str
        The expression, e.g. ``"byres within 5 of ligand"``.

    Returns
    -------
    Node
        The root of the tree (trees are immutable and may be shared).

    Raises
    ------
    SelectionError
        If the text is not a valid expression; the error says where.
    """
    return _Parser(text).parse()


def check_name(name: str) -> None:
    """Make sure a word can name a selection.

    Raises
    ------
    ValueError
        If it is empty, not a plain word (letters, digits and underscores,
        not starting with a digit) or one of the :data:`KEYWORDS`.
    """
    if not _NAME_RE.match(name):
        raise ValueError(f"'{name}' cannot name a selection: use letters, digits and underscores")
    if name.lower() in KEYWORDS:
        raise ValueError(f"'{name}' is a keyword of the selection language")


# -- working it out --------------------------------------------------------------------------

@dataclass
class Context:
    """What an expression may refer to besides the molecule.

    Attributes
    ----------
    selected : iterable of int or None
        Atoms selected in the window, for ``selected``.
    hidden : numpy.ndarray or None
        Boolean mask of the atoms that are not drawn, for ``visible`` and ``hidden``.
    names : mapping of str to str
        Named selections: each name and the expression it stands for.
    columns : mapping of str to array or callable
        Numeric per-atom values to compare, in addition to the molecule's
        columns, such as ``"charge"``. A callable is called only when the
        expression uses the name; it may return None when there are no values.
    """

    selected: Iterable[int] | None = None
    hidden: np.ndarray | None = None
    names: Mapping[str, str] = field(default_factory=dict)
    columns: Mapping[str, Any] = field(default_factory=dict)


def evaluate(mol: Molecule, expression: str | Node, context: Context | None = None) -> np.ndarray:
    """Work out which atoms an expression picks.

    Parameters
    ----------
    mol : Molecule
        Structure.
    expression : str or Node
        The expression, or a tree from :func:`parse`.
    context : Context, optional
        Selection, hidden atoms, named selections and extra columns.

    Returns
    -------
    numpy.ndarray
        Boolean mask, one value per atom.

    Raises
    ------
    SelectionError
        If the expression is not valid, names an unknown selection or
        property, or holds an invalid SMARTS pattern.
    """
    text = expression if isinstance(expression, str) else ""
    node = parse(expression) if isinstance(expression, str) else expression
    return _Evaluator(mol, context or Context(), text).mask(node)


def find(mol: Molecule, expression: str | Node, context: Context | None = None) -> np.ndarray:
    """Like :func:`evaluate`, but return the indices (0-based, ascending) of the atoms picked."""
    return np.flatnonzero(evaluate(mol, expression, context))


def residue_ids(mol: Molecule) -> np.ndarray:
    """Give each atom the number of its residue; atoms of one residue share the number.

    Residues are told apart by chain, residue number, insertion code and
    residue name (:data:`~miles.core.molecule.RESIDUE_KEY`); atoms without
    residue data are grouped by molecule (connected component).

    Returns
    -------
    numpy.ndarray
        One integer per atom; the numbers start at 0 but are not in file order.
    """
    return _Evaluator(mol, Context(), "").residue_ids()


class _Evaluator:
    """Works out one expression on one molecule, keeping what several parts of it share."""

    def __init__(self, mol: Molecule, context: Context, text: str) -> None:
        self.mol = mol
        self.n = len(mol)
        self.context = context
        self.text = text
        self.cache: dict[Any, Any] = {}
        self.names_in_use: list[str] = []

    def error(self, message: str, node: Node | None) -> SelectionError:
        return SelectionError(message, self.text, node.position if node is not None and self.text else None)

    def none(self) -> np.ndarray:
        return np.zeros(self.n, dtype=bool)

    def mask(self, node: Node) -> np.ndarray:
        match node:
            case Either(left=left, right=right):
                return self.mask(left) | self.mask(right)
            case Both(left=left, right=right):
                # a distance search next to "and" only needs to test the atoms the other side keeps
                if isinstance(right, Within):
                    first = self.mask(left)
                    return first & self.within(self.mask(right.inner), right.distance, right.exclude, first)
                if isinstance(left, Within):
                    second = self.mask(right)
                    return self.within(self.mask(left.inner), left.distance, left.exclude, second) & second
                return self.mask(left) & self.mask(right)
            case Not(inner=inner):
                return ~self.mask(inner)
            case Flag(group=group):
                return self.flag(group, node)
            case Values():
                return self.values(node)
            case Compare():
                return self.compare(node)
            case Within(distance=distance, inner=inner, exclude=exclude):
                return self.within(self.mask(inner), distance, exclude)
            case Expand(to=to, inner=inner):
                return self.expand(self.mask(inner), to)
            case Smarts(pattern=pattern):
                return self.smarts(pattern, node)
            case Named(name=name):
                return self.named(name, node)
        raise TypeError(f"not a selection node: {node!r}")

    # -- columns, run by run -----------------------------------------------------------------

    def column(self, name: str) -> np.ndarray:
        key = ("column", name)
        if key not in self.cache:
            self.cache[key] = self.mol.column(name)
        return self.cache[key]

    def runs(self, name: str) -> tuple[np.ndarray, np.ndarray, list[Any]]:
        """Return the runs of equal values of a column: run of each atom, first atom of each run, values."""
        key = ("runs", name)
        if key not in self.cache:
            values = self.column(name)
            change = np.ones(self.n, dtype=bool)
            if self.n > 1:
                change[1:] = values[1:] != values[:-1]
            starts = np.flatnonzero(change)
            self.cache[key] = (np.cumsum(change) - 1, starts, values[starts].tolist())
        return self.cache[key]

    def where(self, name: str, test: Callable[[Any], bool]) -> np.ndarray:
        """Atoms whose value of a string column passes a test, deciding each distinct value once."""
        if not self.n:
            return self.none()
        run, _, values = self.runs(name)
        verdicts: dict[Any, bool] = {}
        per_run = np.fromiter((verdicts[v] if v in verdicts else verdicts.setdefault(v, bool(test(v)))
                               for v in values), dtype=bool, count=len(values))
        return per_run[run]

    def resname_in(self, names: frozenset[str] | set[str]) -> np.ndarray:
        return self.where("resname", names.__contains__)

    def has_residue(self) -> np.ndarray:
        key = "has_residue"
        if key not in self.cache:
            self.cache[key] = self.where("resname", bool)
        return self.cache[key]

    def group_ids(self, names: tuple[str, ...]) -> np.ndarray:
        """Give atoms that share the values of some columns one group number (by runs, then by value)."""
        if not self.n:
            return np.zeros(0, dtype=np.int64)
        change = np.zeros(self.n, dtype=bool)
        change[0] = True
        for name in names:
            values = self.column(name)
            change[1:] |= values[1:] != values[:-1]
        starts = np.flatnonzero(change)
        run = np.cumsum(change) - 1
        keys = zip(*(self.column(name)[starts].tolist() for name in names), strict=True)
        seen: dict[tuple[Any, ...], int] = {}
        per_run = np.fromiter((seen.setdefault(k, len(seen)) for k in keys), dtype=np.int64,
                              count=len(starts))
        return per_run[run]

    def components(self) -> np.ndarray:
        key = "components"
        if key not in self.cache:
            self.cache[key] = self.mol.component_labels()
        return self.cache[key]

    def grouped(self, names: tuple[str, ...]) -> np.ndarray:
        """Group numbers over residue atoms by ``names``; atoms without residue data by molecule."""
        key = ("grouped", names)
        if key not in self.cache:
            ids = self.group_ids(names)
            loose = ~self.has_residue()
            if loose.any():
                ids = ids.copy()
                ids[loose] = ids.max(initial=-1) + 1 + self.components()[loose]
            self.cache[key] = ids
        return self.cache[key]

    def residue_ids(self) -> np.ndarray:
        return self.grouped(RESIDUE_KEY)

    def degrees(self) -> np.ndarray:
        key = "degrees"
        if key not in self.cache:
            self.cache[key] = self.mol.degrees()
        return self.cache[key]

    def bonded_to(self, inner: np.ndarray) -> np.ndarray:
        pairs = self.mol.bond_pairs
        out = self.none()
        if len(pairs):
            a, b = pairs[:, 0], pairs[:, 1]
            out[b[inner[a]]] = True
            out[a[inner[b]]] = True
        return out

    # -- groups ------------------------------------------------------------------------------

    def flag(self, group: str, node: Node) -> np.ndarray:
        key = ("flag", group)
        if key in self.cache:
            return self.cache[key]
        numbers = self.mol.numbers
        if group == "all":
            out = np.ones(self.n, dtype=bool)
        elif group == "none":
            out = self.none()
        elif group == "protein":
            out = self.resname_in(AMINO_ACIDS)
        elif group == "nucleic":
            out = self.resname_in(NUCLEOTIDES)
        elif group == "water":
            out = self.resname_in(WATERS)
        elif group == "hydrogen":
            out = numbers == 1
        elif group == "heavy":
            out = numbers > 1
        elif group == "metal":
            out = _METAL[numbers]
        elif group == "hetero":
            out = self.column("hetero").copy()
        elif group == "ion":
            out = self.ions()
        elif group == "ligand":
            polymer = self.flag("protein", node) | self.flag("nucleic", node)
            out = ~polymer & ~self.flag("water", node) & ~self.flag("ion", node) & (numbers > 0)
        elif group == "backbone":
            protein, nucleic = self.flag("protein", node), self.flag("nucleic", node)
            out = ((protein & self.where("name", BACKBONE_NAMES.__contains__))
                   | (nucleic & self.where("name", NUCLEIC_BACKBONE_NAMES.__contains__)))
        elif group == "sidechain":
            out = (self.flag("protein", node) & ~self.flag("backbone", node)
                   & ~self.where("name", BACKBONE_HYDROGENS.__contains__))
        elif group == "polar":
            nitrogen_oxygen = (numbers == 7) | (numbers == 8)
            out = nitrogen_oxygen | ((numbers == 1) & self.bonded_to(nitrogen_oxygen))
        elif group == "apolar":
            out = ~self.flag("polar", node)
        elif group == "aromatic":
            out = self.none()
            aromatic = self.mol.bond_pairs[self.mol.bond_orders == 1.5]
            out[aromatic.reshape(-1)] = True
        elif group in ("helix", "sheet", "loop"):
            code = {"helix": "H", "sheet": "E", "loop": ""}[group]
            out = self.where("ss", lambda value: value == code)
            if group == "loop":
                out &= self.flag("protein", node)
        elif group == "selected":
            out = self.none()
            if self.context.selected is not None:
                chosen = np.fromiter((int(i) for i in self.context.selected), dtype=np.intp)
                out[chosen[(chosen >= 0) & (chosen < self.n)]] = True
        elif group in ("visible", "hidden"):
            hidden = self.context.hidden
            if hidden is None or len(hidden) != self.n:
                hidden = self.none()
            hidden = np.asarray(hidden, dtype=bool)
            out = ~hidden if group == "visible" else hidden.copy()
        else:
            raise self.error(f"Unknown group '{group}'", node)
        self.cache[key] = out
        return out

    def ions(self) -> np.ndarray:
        """Metal and halide atoms that stand alone: no bonds, or alone in their residue."""
        numbers = self.mol.numbers
        candidates = _METAL[numbers] | np.isin(numbers, list(HALOGENS))
        if not candidates.any():
            return candidates
        alone = self.degrees() == 0
        residues = self.residue_ids()
        heavy = numbers > 1
        sizes = np.bincount(residues[heavy], minlength=int(residues.max(initial=0)) + 1)
        alone |= self.has_residue() & (sizes[residues] == 1)
        return candidates & alone

    # -- values and numbers ------------------------------------------------------------------

    def values(self, node: Values) -> np.ndarray:
        if node.field == "element":
            wanted = np.zeros(len(BY_Z), dtype=bool)
            for _, symbol, _ in node.values:
                wanted[symbol_number(symbol)] = True
            return wanted[self.mol.numbers]
        if node.field in ("resid", "index"):
            numbers = (self.column("resseq") if node.field == "resid"
                       else np.arange(1, self.n + 1))
            out = self.none()
            for _, low, high in node.values:
                out |= (numbers >= low) & (numbers <= high)
            return out
        forms = self.exact_forms(node)
        if node.field == "name" and forms is not None and len(forms) <= 4:
            # atom names change from one atom to the next: compare the whole column directly
            column = self.column("name")
            out = self.none()
            for form in forms:
                out |= column == form
            return out
        return self.where(node.field, self.matcher(node))

    def exact_forms(self, node: Values) -> set[str] | None:
        """Return the values a name column may hold, or None when a value has wildcards."""
        forms: set[str] = set()
        for _, text, quoted in node.values:
            if not quoted and ("*" in text or "?" in text):
                return None
            forms |= {text, text.upper()} if node.field != "chain" else {text}
        return forms

    def matcher(self, node: Values) -> Callable[[Any], bool]:
        """Test for one value of a name column: exact, or with * and ? as wildcards.

        Names and residue names also match the capitals of what was typed
        (``resname his`` finds HIS); chain identifiers are matched exactly,
        since large structures use both ``A`` and ``a``.
        """
        exact: set[str] = set()
        patterns: list[re.Pattern[str]] = []
        loose = node.field != "chain"
        for _, text, quoted in node.values:
            forms = {text, text.upper()} if loose else {text}
            if not quoted and ("*" in text or "?" in text):
                patterns.extend(re.compile(fnmatch.translate(form)) for form in forms)
            else:
                exact |= forms
        return lambda value: value in exact or any(p.match(value) for p in patterns)

    def numeric(self, name: str, node: Node) -> np.ndarray:
        column = COMPARABLE.get(name.lower(), name)
        if column in ("x", "y", "z"):
            return self.mol.coords[:, "xyz".index(column)]
        if column == "resid":
            return self.column("resseq")
        if column == "index":
            return np.arange(1, self.n + 1)
        if column == "mass":
            return _MASS[self.mol.numbers]
        if column in self.context.columns:
            values = self.context.columns[column]
            values = values() if callable(values) else values
            if values is None:
                raise self.error(f"There are no {column} values for this structure", node)
            return np.asarray(values, dtype=float)
        if column in ("charge", "b_factor", "occupancy") and not self.mol.has_property(column):
            raise self.error(f"There are no {column} values for this structure", node)
        if self.mol.has_property(column):
            values = self.mol.column(column)
            if not np.issubdtype(values.dtype, np.number) and values.dtype != np.bool_:
                raise self.error(f"'{name}' holds text, not numbers", node)
            return values
        numeric = (p for p in self.mol.properties
                   if np.issubdtype(self.mol.property_spec(p).dtype, np.number))
        known = sorted({*COMPARABLE.values(), *self.context.columns, *numeric})
        hint = difflib.get_close_matches(name.lower(), known, n=1)
        suggestion = (f" Did you mean '{hint[0]}'?" if hint
                      else f" Numbers can be compared for: {', '.join(known)}.")
        raise self.error(f"No number called '{name}'.{suggestion}", node)

    def compare(self, node: Compare) -> np.ndarray:
        values = self.numeric(node.column, node)
        with np.errstate(invalid="ignore"):
            return np.asarray(_COMPARE[node.operator](values, node.value), dtype=bool)

    # -- distances and structure -------------------------------------------------------------

    def within(self, inner: np.ndarray, distance: float, exclude: bool,
               domain: np.ndarray | None = None) -> np.ndarray:
        """Atoms within ``distance`` of ``inner`` (only those in ``domain`` are tested when it is given)."""
        out = inner.copy()
        if inner.any() and distance > 0:
            coords = self.mol.coords
            centers = coords[inner]
            low, high = centers.min(axis=0) - distance, centers.max(axis=0) + distance
            near = ~inner if domain is None else domain & ~inner
            near &= np.all((coords >= low) & (coords <= high), axis=1)
            candidates = np.flatnonzero(near)
            if len(candidates):
                out[candidates[near_mask(coords[candidates], centers, distance)]] = True
        if exclude:
            out &= ~inner
        return out

    def expand(self, inner: np.ndarray, to: str) -> np.ndarray:
        if to == "bonded":
            return self.bonded_to(inner) & ~inner
        if not inner.any():
            return inner
        ids = self.residue_ids() if to == "residue" else self.grouped(("chain",))
        hit = np.zeros(int(ids.max(initial=-1)) + 1, dtype=bool)
        hit[ids[inner]] = True
        return hit[ids]

    def smarts(self, pattern: str, node: Node) -> np.ndarray:
        try:
            found = match_pattern(self.mol, pattern)
        except PatternError as exc:
            raise self.error(str(exc), node) from None
        out = self.none()
        out[np.asarray(found.atoms, dtype=np.intp)] = True
        return out

    def named(self, name: str, node: Node) -> np.ndarray:
        names = self.context.names
        key = name if name in names else next((k for k in names if k.lower() == name.lower()), None)
        if key is None:
            choices = sorted({*FLAGS, *FIELDS, *names})
            hint = difflib.get_close_matches(name.lower(), [c.lower() for c in choices], n=1)
            suggestion = f" Did you mean '{hint[0]}'?" if hint else ""
            raise self.error(f"Unknown word '{name}'.{suggestion}", node)
        if key in self.names_in_use:
            raise self.error(f"The selection '{key}' refers to itself", node)
        cache_key = ("named", key)
        if cache_key not in self.cache:
            try:
                tree = parse(names[key])
            except SelectionError as exc:
                raise self.error(f"The selection '{key}' is not valid: {exc.message}", node) from None
            self.names_in_use.append(key)
            try:
                inner = _Evaluator(self.mol, self.context, names[key])
                inner.cache = self.cache
                inner.names_in_use = self.names_in_use
                try:
                    self.cache[cache_key] = inner.mask(tree)
                except SelectionError as exc:
                    cycle = exc.message.endswith("refers to itself")
                    message = exc.message if cycle else f"In the selection '{key}': {exc.message}"
                    raise self.error(message, node) from None
            finally:
                self.names_in_use.pop()
        return self.cache[cache_key]
