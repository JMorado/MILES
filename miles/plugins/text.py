"""Text about atoms that plugins write in panels, tables and labels."""

from miles.core.molecule import Molecule


def atom_name(mol: Molecule, i: int) -> str:
    """Name an atom: ``A/ASP25:OD1`` for an atom of a residue, else its element and number, ``C12``.

    The chain is left out when the residue has none, and an insertion
    code follows the residue number (``A/GLY52A:CA``).
    """
    atom = mol.atoms[i]
    if not atom.resname:
        return f"{atom.symbol}{i + 1}"
    chain = f"{atom.chain}/" if atom.chain else ""
    icode = str(getattr(atom, "icode", "") or "").strip()
    name = atom.name.strip() or f"{atom.symbol}{i + 1}"
    return f"{chain}{atom.resname}{atom.resseq}{icode}:{name}"
