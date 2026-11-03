"""Writing PDB files from the molecule's columns.

Open Babel's PDB writer renumbers residues and leaves out alternate
locations, occupancies and B-factors, so PDB files are written here from
the columns directly. Every atom becomes an ATOM record (or HETATM for
hetero residues and atoms without residue data) with its residue, chain,
insertion code, alternate location, occupancy, B-factor, element and
formal charge; TER closes each polymer chain, and CONECT records give the
bonds that residue templates cannot: those of hetero residues and those
between residues other than the links of a polymer chain (disulfides, for
example). A double or triple bond is written by repeating the partner, as
readers such as Open Babel expect.

Atoms without residue data are written as residue ``UNL`` (one residue per
molecule) with names made from the element and a count, as in ``C1``.
Reading PDB files stays with Open Babel (see :mod:`miles.files.openbabel`).
"""

import math
from collections import Counter

import numpy as np

from miles.core.charges import formal_charges
from miles.core.molecule import Molecule
from miles.core.protein import AMINO_ACIDS, NUCLEOTIDES
from miles.files.common import FileFormatError, JobInfo

MAX_ATOMS = 99_999
"""Atoms a PDB file can number; larger structures need mmCIF."""

_POLYMER = AMINO_ACIDS | NUCLEOTIDES
_LINKS = {("C", "N"), ("O3'", "P"), ("O3*", "P")}
"""Atom names of the bonds that join consecutive residues of a chain (residue templates cover them)."""


def _atom_name(name: str, symbol: str) -> str:
    """Four-character atom name field: one-letter elements start in column 14, as the format asks."""
    if len(name) >= 4:
        return name[:4]
    if len(symbol) == 1 and name.upper().startswith(symbol.upper()):
        return f" {name:<3s}"
    return f"{name:<4s}"


def _charge_field(charge: int) -> str:
    return "  " if charge == 0 else f"{abs(charge)}{'+' if charge > 0 else '-'}"


def pdb_text(mol: Molecule, title: str = "") -> str:
    """Write a structure as the text of a PDB file.

    Parameters
    ----------
    mol : Molecule
        Structure to write; its residue data, alternate locations, the
        ``occupancy`` and ``b_factor`` properties (from a PDB file) and the
        formal charges are kept.
    title : str, optional
        Written as a TITLE record.

    Returns
    -------
    str
        The file contents.

    Raises
    ------
    FileFormatError
        If the structure has more atoms, or residue numbers, than the
        format can hold.
    """
    n = len(mol)
    if n > MAX_ATOMS:
        raise FileFormatError(f"A PDB file holds at most {MAX_ATOMS:,} atoms; this structure has {n:,}. "
                              "Save it as mmCIF (.cif) instead.")
    symbols = mol.symbols
    names = mol.column("name").tolist()
    resnames = mol.column("resname").tolist()
    resseqs = mol.column("resseq").tolist()
    if resseqs and (max(resseqs) > 9999 or min(resseqs) < -999):
        raise FileFormatError("A PDB file holds residue numbers from -999 to 9999. "
                              "Save it as mmCIF (.cif) instead.")
    chains = mol.column("chain").tolist()
    icodes = mol.column("icode").tolist()
    altlocs = mol.column("altloc").tolist()
    hetero = mol.column("hetero").tolist()
    occupancy = mol.column("occupancy") if mol.has_property("occupancy") else np.full(n, math.nan)
    b_factor = mol.column("b_factor") if mol.has_property("b_factor") else np.full(n, math.nan)
    charges = formal_charges(mol).tolist() if n else []
    molecules = mol.component_labels().tolist() if n else []

    # atoms without residue data: residue UNL per molecule, names from element and count
    counts: Counter[tuple[int, str]] = Counter()
    for i in range(n):
        if not resnames[i]:
            counts[(molecules[i], symbols[i])] += 1
            names[i] = f"{symbols[i]}{counts[(molecules[i], symbols[i])]}"
            resnames[i], resseqs[i], hetero[i] = "UNL", molecules[i] + 1, True

    lines = [f"TITLE     {title[:70]}"] if title else []
    serial: dict[int, int] = {}
    number = 0
    for i, (x, y, z) in enumerate(mol.coords.tolist()):
        number += 1
        serial[i] = number
        record = "HETATM" if hetero[i] or resnames[i] not in _POLYMER else "ATOM  "
        occ = 1.0 if math.isnan(occupancy[i]) else float(occupancy[i])
        bf = 0.0 if math.isnan(b_factor[i]) else float(b_factor[i])
        chain, icode, altloc = (chains[i] or " ")[:1], (icodes[i] or " ")[:1], (altlocs[i] or " ")[:1]
        lines.append(f"{record}{number:5d} {_atom_name(names[i], symbols[i])}{altloc}"
                     f"{resnames[i][:3]:>3s} {chain}{resseqs[i]:4d}{icode}   "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}          "
                     f"{symbols[i].upper():>2s}{_charge_field(charges[i])}")
        last = i + 1 == n
        end_of_chain = record == "ATOM  " and (
            last or chains[i + 1] != chains[i] or hetero[i + 1] or resnames[i + 1] not in _POLYMER)
        if end_of_chain:
            number += 1
            lines.append(f"TER   {number:5d}      {resnames[i][:3]:>3s} {(chains[i] or ' ')[:1]}"
                         f"{resseqs[i]:4d}{(icodes[i] or ' ')[:1]}")

    partners: dict[int, list[int]] = {}
    residue = list(zip(chains, resseqs, icodes, resnames, strict=True))
    for (i, j), order in zip(mol.bond_pairs.tolist(), mol.bond_orders.tolist(), strict=True):
        templated = (not hetero[i] and not hetero[j] and resnames[i] in _POLYMER and resnames[j] in _POLYMER
                     and (residue[i] == residue[j] or (names[i], names[j]) in _LINKS
                          or (names[j], names[i]) in _LINKS))
        if templated:
            continue
        repeat = 1 if order == 1.5 else max(1, round(order))
        partners.setdefault(i, []).extend([j] * repeat)
        partners.setdefault(j, []).extend([i] * repeat)
    for i in sorted(partners):
        others = [serial[j] for j in partners[i]]
        for k in range(0, len(others), 4):
            lines.append(f"CONECT{serial[i]:5d}" + "".join(f"{s:5d}" for s in others[k:k + 4]))
    lines.append("END")
    return "\n".join(lines) + "\n"


def write(path: str, mol: Molecule, info: JobInfo) -> None:
    """Write a PDB file (the format registry's writer for ``.pdb`` and ``.ent``).

    Parameters
    ----------
    path : str
        Destination.
    mol : Molecule
        Structure to write.
    info : JobInfo
        Supplies the title.

    Raises
    ------
    FileFormatError
        If the structure does not fit the format, or the file cannot be written.
    """
    text = pdb_text(mol, info.title)
    try:
        with open(path, "w", encoding="ascii", errors="replace") as fh:
            fh.write(text)
    except OSError as exc:
        raise FileFormatError(f"Could not write {path}: {exc}") from exc
