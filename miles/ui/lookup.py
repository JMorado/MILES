"""Looking up structures by SMILES, or through structure sources, without blocking the interface.

SMILES strings are built here. Everything else, such as names and database
codes, goes to a :class:`StructureSource` that a plugin added: the PDB and
PubChem sources are the ``structure_lookup`` built-in plugin.
"""

from collections.abc import Callable
from dataclasses import dataclass

from PySide6.QtCore import QObject, Signal

from miles.core.chem import ChemError, from_smiles, is_smiles
from miles.core.molecule import Molecule
from miles.files.common import FileFormatError, JobInfo
from miles.ui.tasks import executor


@dataclass(frozen=True)
class StructureSource:
    """Somewhere structures can be fetched from by a code or a name.

    Attributes
    ----------
    name : str
        Name of the source, e.g. ``"PDB"``.
    fetch : callable
        ``fetch(text) -> (Molecule, JobInfo)``; runs on a worker thread, so
        it must not touch the interface. Raise
        :class:`~miles.files.common.FileFormatError` with a message for
        the user when nothing is found.
    matches : callable or None
        ``matches(text) -> bool``: the texts this source is for, tried
        before SMILES (e.g. four-character PDB codes). None makes it a
        fallback for text that is neither matched by another source nor
        SMILES, such as compound names.
    what : str
        What users type, for the hint in the search box, e.g. ``"PDB code"``.
    example : str
        An example of it, e.g. ``"1crn"``.
    whole : bool
        The result opens in a tab of its own even from the fragment search;
        otherwise it becomes the fragment to place there (and opens in a new
        tab from the command line and from scripts).
    message : callable or None
        ``message(text) -> str`` shown while fetching, e.g. ``"Fetching
        1CRN from the PDB ..."``; the default is ``"Building <text> ..."``.
    """

    name: str
    fetch: Callable[[str], tuple[Molecule, JobInfo]]
    matches: Callable[[str], bool] | None = None
    what: str = ""
    example: str = ""
    whole: bool = True
    message: Callable[[str], str] | None = None


def choose_source(text: str, sources: list[StructureSource]) -> StructureSource | None:
    """Return the source that handles a text, or None for SMILES and for text no source takes.

    Parameters
    ----------
    text : str
        What was typed.
    sources : list of StructureSource
        The sources, in the order they were added.

    Returns
    -------
    StructureSource or None
        The first source whose ``matches`` accepts the text; else, unless
        the text is SMILES, the first source without ``matches``.
    """
    text = text.strip()
    for source in sources:
        if source.matches is not None and source.matches(text):
            return source
    if is_smiles(text):
        return None
    return next((source for source in sources if source.matches is None), None)


def search_hint(sources: list[StructureSource]) -> str:
    """Return the placeholder of a search box: what can be typed, with examples."""
    kinds = [("SMILES", "CC(=O)O")] + [(s.what, s.example) for s in sources if s.what]
    names = [what for what, _ in kinds]
    joined = names[0] if len(names) == 1 else ", ".join(names[:-1]) + " or " + names[-1]
    return f"{joined}, e.g. {', '.join(example for _, example in kinds if example)}"


def look_up(text: str, sources: list[StructureSource]) -> tuple[Molecule, JobInfo | None, int | None]:
    """Build a SMILES string or fetch a text from a source, in the calling thread.

    Parameters
    ----------
    text : str
        SMILES, or what a source takes.
    sources : list of StructureSource
        The sources to choose from.

    Returns
    -------
    mol : Molecule
        The structure.
    info : JobInfo or None
        For a source that gives whole structures, its information.
    attach : int or None
        Hydrogen standing in for ``*`` in SMILES, if any.

    Raises
    ------
    FileFormatError
        If the text cannot be built or found.
    """
    text = text.strip()
    if not text:
        raise FileFormatError("Nothing to build")
    source = choose_source(text, sources)
    if source is not None:
        mol, info = source.fetch(text)
        return mol, info if source.whole else None, None
    if not is_smiles(text):
        raise FileFormatError(f'"{text}" is not SMILES, and no structure source can look it up')
    try:
        mol, attach = from_smiles(text)
    except ChemError as exc:
        raise FileFormatError(str(exc)) from exc
    return mol, None, attach


class StructureLookup(QObject):
    """Builds or fetches a structure on a worker thread.

    Parameters
    ----------
    sources : list of StructureSource
        The window's sources; the list is shared, so sources added later
        are used too.

    Signals
    -------
    fragment(mol, attach, smiles, text)
        A small molecule from SMILES or a name, with its attachment hydrogen.
    structure(mol, info)
        A whole structure, e.g. an entry of the PDB.
    failed(message)
        Nothing could be found.
    """

    fragment = Signal(object, object, str, str)
    structure = Signal(object, object)
    failed = Signal(str)

    def __init__(self, sources: list[StructureSource] | None = None) -> None:
        super().__init__()
        self.sources = sources if sources is not None else []

    def start(self, text: str) -> None:
        """Look up ``text`` on a worker thread; one of the signals follows on the interface thread."""
        executor().submit(look_up, text, list(self.sources), owner=self, failed=self.failed.emit,
                          done=lambda found: self._found(text, *found))

    def _found(self, text: str, mol: Molecule, info: JobInfo | None, attach: int | None) -> None:
        if info is not None:
            self.structure.emit(mol, info)
        else:
            self.fragment.emit(mol, attach, text.strip(), text)
