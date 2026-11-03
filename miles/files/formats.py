"""Reading and writing structures through a registry of formats.

Sessions (``.miles``) have their own reader and writer, Gaussian input
(``.gjf``, ``.com``) its own reader and PDB files (``.pdb``, ``.ent``) their
own writer; every other extension falls through to Open Babel.
"""

import re
from collections.abc import Callable

from miles.core.chem import ChemError, from_smiles
from miles.core.molecule import Molecule
from miles.files import gaussian, openbabel, pdb, session
from miles.files.common import FileFormatError, JobInfo, extension

Reader = Callable[[str], tuple[Molecule, JobInfo]]
Writer = Callable[[str, Molecule, JobInfo], None]


class FormatRegistry:
    """Readers and writers keyed by file extension, with Open Babel as fallback."""

    def __init__(self) -> None:
        self._readers: dict[str, Reader] = {}
        self._writers: dict[str, Writer] = {}

    def register(self, ext: str, reader: Reader | None = None, writer: Writer | None = None) -> None:
        """Add a reader and/or writer for an extension.

        Parameters
        ----------
        ext : str
            Extension without the dot.
        reader, writer : callable, optional
            Functions following the ``Reader`` and ``Writer`` signatures.
        """
        if reader is not None:
            self._readers[ext] = reader
        if writer is not None:
            self._writers[ext] = writer

    def unregister(self, ext: str, reader: Reader | None = None, writer: Writer | None = None) -> None:
        """Remove a reader and/or writer of an extension, if they are the ones registered.

        Parameters
        ----------
        ext : str
            Extension without the dot.
        reader, writer : callable, optional
            The functions given to :meth:`register`; a later registration
            for the extension is left alone.
        """
        if reader is not None and self._readers.get(ext) is reader:
            del self._readers[ext]
        if writer is not None and self._writers.get(ext) is writer:
            del self._writers[ext]

    def reader(self, path: str) -> Reader:
        """Reader for a path (Open Babel unless another is registered)."""
        return self._readers.get(extension(path), openbabel.read)

    def writer(self, path: str) -> Writer:
        """Writer for a path (Open Babel unless another is registered)."""
        return self._writers.get(extension(path), openbabel.write)

    def read_extensions(self) -> list[str]:
        """Every readable extension."""
        return sorted(set(self._readers) | set(openbabel.read_extensions()))

    def write_extensions(self) -> list[str]:
        """Every writable extension."""
        return sorted(set(self._writers) | set(openbabel.write_extensions()))

    def own_writers(self) -> list[str]:
        """Extensions written by a registered writer rather than Open Babel (e.g. by plugins)."""
        return sorted(self._writers)


REGISTRY = FormatRegistry()
REGISTRY.register(session.EXTENSION, session.read, session.write)
for _ext in ("gjf", "com"):
    REGISTRY.register(_ext, reader=gaussian.read)   # Open Babel only writes Gaussian input
for _ext in ("pdb", "ent"):
    REGISTRY.register(_ext, writer=pdb.write)       # keeps residues, B-factors and alternate locations

COMMON_WRITE = ["xyz", "sdf", "mol", "mol2", "pdb", "cif", "cml"]


def read_filter() -> str:
    """File dialog filter listing every readable format."""
    exts = " ".join(f"*.{e}" for e in REGISTRY.read_extensions())
    return f"Structures ({exts} *.log *.out);;Sessions (*.miles);;All files (*)"


def write_filter() -> str:
    """File dialog filter with the common writable formats first, then those added by plugins."""
    common = ";;".join(f"{e.upper()} (*.{e})" for e in COMMON_WRITE)
    added = "".join(f";;{e.upper()} (*.{e})" for e in REGISTRY.own_writers()
                    if e != session.EXTENSION and e not in COMMON_WRITE)
    return f"{common}{added};;Session with styles (*.miles);;Any other format (*)"


def read_file(path: str) -> tuple[Molecule, JobInfo]:
    """Read a structure file.

    Parameters
    ----------
    path : str
        File to read; the extension selects the format.

    Returns
    -------
    mol : Molecule
        The structure.
    info : JobInfo
        Title, charge, frames and (for sessions) the saved view.
    """
    return REGISTRY.reader(path)(path)


def write_file(path: str, mol: Molecule, info: JobInfo | None = None) -> None:
    """Write a structure file in the format given by its extension.

    Parameters
    ----------
    path : str
        Destination file.
    mol : Molecule
        Structure to write.
    info : JobInfo, optional
        Title, frames and view (frames and view are only kept by sessions).
    """
    REGISTRY.writer(path)(path, mol, info or JobInfo())


def sniff_format(text: str) -> str | None:
    """Guess the Open Babel format of pasted text.

    Parameters
    ----------
    text : str
        Text from the clipboard.

    Returns
    -------
    str or None
        Format code, ``"lines"`` for bare ``symbol x y z`` rows, ``"smiles"``
        for a single token, or ``None`` when nothing fits.
    """
    t = text.strip()
    first = t.splitlines()[0].split() if t else []
    if len(first) == 1 and first[0].isdigit():
        return "xyz"
    if "@<TRIPOS>" in t:
        return "mol2"
    if "V2000" in t or "V3000" in t:
        return "sdf"
    if re.search(r"^(ATOM  |HETATM)", t, re.M):
        return "pdb"
    if t.startswith("data_"):
        return "mmcif"
    if t.startswith("<") and "molecule" in t[:400]:
        return "cml"
    if t.startswith(("%", "#")):
        return "gjf"
    if re.search(r"^\s*[A-Za-z]{1,3}\d*[\s,]+-?\d", t, re.M):
        return "lines"
    if len(t.splitlines()) == 1 and " " not in t:
        return "smiles"
    return None


def parse_text(text: str) -> Molecule:
    """Structure from pasted text.

    Parameters
    ----------
    text : str
        XYZ, MOL/SDF, MOL2, PDB, mmCIF, CML, Gaussian input, plain
        ``symbol x y z`` lines (extra columns ignored) or a SMILES string.

    Returns
    -------
    Molecule
        The structure, with bonds perceived where the text had none.

    Raises
    ------
    FileFormatError
        If the text is not a structure.
    """
    fmt = sniff_format(text)
    if fmt is None:
        raise FileFormatError("No structure found in the pasted text")
    if fmt == "smiles":
        try:
            return from_smiles(text.strip())[0]
        except ChemError as exc:
            raise FileFormatError(str(exc)) from exc
    if fmt == "gjf":
        return gaussian.read_text(text)[0]
    if fmt == "lines":
        rows = []
        for line in text.strip().splitlines():
            parts = line.replace(",", " ").split()
            numbers = [p for p in parts[1:] if _is_number(p)]
            if parts and len(numbers) >= 3:
                rows.append(f"{parts[0].rstrip('0123456789')} {' '.join(numbers[-3:])}")
        text = f"{len(rows)}\n\n" + "\n".join(rows) + "\n"
        fmt = "xyz"
    return openbabel.read_text(text, fmt)


def _is_number(token: str) -> bool:
    try:
        float(token.replace("D", "E"))
    except ValueError:
        return False
    return True
