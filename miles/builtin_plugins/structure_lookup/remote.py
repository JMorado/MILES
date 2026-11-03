"""Structures from the internet: RCSB PDB entries and PubChem compounds."""

import json
import re
import urllib.error
import urllib.parse
import urllib.request

from miles import ChemError, FileFormatError, JobInfo, Molecule, from_smiles, is_smiles, read_text

RCSB = "https://files.rcsb.org/download/{code}.{ext}"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/{what}"
PDB_CODE = re.compile(r"[0-9][A-Za-z0-9]{3}")


def is_pdb_code(text: str) -> bool:
    """Whether text looks like a 4-character PDB identifier such as ``1crn``.

    Parameters
    ----------
    text : str
        Candidate code.

    Returns
    -------
    bool
        True for a digit followed by three letters or digits.
    """
    return PDB_CODE.fullmatch(text.strip()) is not None


def _download(url: str, timeout: float) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def fetch_pdb(code: str, timeout: float = 20.0) -> tuple[Molecule, JobInfo]:
    """Download an entry from the RCSB Protein Data Bank.

    Parameters
    ----------
    code : str
        PDB identifier, e.g. ``"4hhb"``.
    timeout : float, optional
        Network timeout in seconds.

    Returns
    -------
    mol : Molecule
        The first model, with residue data.
    info : JobInfo
        Title set to the code.

    Raises
    ------
    FileFormatError
        If the code is malformed, unknown or the download fails.
    """
    code = code.strip().upper()
    if not is_pdb_code(code):
        raise FileFormatError(f"{code} is not a PDB code")
    for ext, fmt in (("pdb", "pdb"), ("cif", "mmcif")):
        try:
            text = _download(RCSB.format(code=code, ext=ext), timeout)
        except urllib.error.HTTPError:
            continue
        except OSError as exc:
            raise FileFormatError(f"Could not download {code}: {exc}") from exc
        return read_text(text, fmt), JobInfo(title=code)
    raise FileFormatError(f"{code} was not found in the PDB")


def fetch_compound(name: str, timeout: float = 8.0) -> Molecule:
    """3D structure of a compound by name from PubChem.

    PubChem's computed 3D conformer is used when there is one; otherwise the
    structure is built from PubChem's SMILES.

    Parameters
    ----------
    name : str
        Common or systematic name, e.g. ``"caffeine"``.
    timeout : float, optional
        Network timeout in seconds.

    Returns
    -------
    Molecule
        The compound with hydrogens.

    Raises
    ------
    FileFormatError
        If PubChem does not know the name or cannot be reached.
    """
    quoted = urllib.parse.quote(" ".join(name.split()))
    try:
        sdf = _download(PUBCHEM.format(name=quoted, what="SDF?record_type=3d"), timeout)
        return read_text(sdf, "sdf")
    except (urllib.error.HTTPError, FileFormatError):
        pass
    except OSError as exc:
        raise FileFormatError(f"Could not reach PubChem: {exc}") from exc
    try:
        url = PUBCHEM.format(name=quoted, what="property/IsomericSMILES/JSON")
        props = json.loads(_download(url, timeout))
        record = props["PropertyTable"]["Properties"][0]
        smiles = record.get("IsomericSMILES") or record.get("SMILES") or record["CanonicalSMILES"]
    except (OSError, KeyError, IndexError, ValueError) as exc:
        raise FileFormatError(f'Could not find "{name}" on PubChem') from exc
    try:
        return from_smiles(smiles)[0]
    except ChemError as exc:
        raise FileFormatError(str(exc)) from exc


def resolve(text: str) -> tuple[Molecule, int | None, str]:
    """Structure for a SMILES string or a compound name.

    Parameters
    ----------
    text : str
        SMILES (``*`` marks an attachment point) or a name.

    Returns
    -------
    mol : Molecule
        The structure.
    attach : int or None
        Hydrogen standing in for ``*``, if any.
    label : str
        The SMILES string or the name, for display.

    Raises
    ------
    FileFormatError
        If the text is neither valid SMILES nor a known name.
    """
    text = text.strip()
    if not text:
        raise FileFormatError("Nothing to build")
    if is_smiles(text):
        try:
            mol, attach = from_smiles(text)
        except ChemError as exc:
            raise FileFormatError(str(exc)) from exc
        return mol, attach, text
    return fetch_compound(text), None, text
