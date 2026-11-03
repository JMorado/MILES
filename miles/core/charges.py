"""Partial atomic charges and the dipole moment.

Charges come from a calculation file (Mulliken charges of a quantum
chemistry output, or the charges stored in a MOL2 file) or are computed
with one of Open Babel's charge models. Both computed models depend only on
the atoms and bonds, not on the coordinates.

Formal charges come from the molecule's ``formal_charge`` column where a
file or SMILES stated them; elsewhere they are worked out from the bonding:
an ammonium nitrogen with four bonds is +1, a carboxylate oxygen with one
bond is -1, and so on. That needs the hydrogens, so it is only done in
molecules that have them (or were built in MILES) and for isolated ions;
add hydrogens before computing partial charges.
"""

from enum import StrEnum

import numpy as np
from openbabel import openbabel as ob

from miles.core.chem import ChemError, silence_openbabel, to_openbabel
from miles.core.molecule import Molecule
from miles.core.units import DEBYE_PER_E_ANGSTROM


class ChargeModel(StrEnum):
    """Where partial charges come from."""

    FILE = "file"            # read with the structure (Mulliken, MOL2 charges)
    GASTEIGER = "gasteiger"  # Gasteiger-Marsili electronegativity equalisation
    MMFF94 = "mmff94"        # charges of the MMFF94 force field

    @property
    def label(self) -> str:
        """Name of the model in menus and messages, e.g. "Gasteiger"."""
        return {"file": "From File", "gasteiger": "Gasteiger", "mmff94": "MMFF94"}[self.value]


def formal_charges(mol: Molecule) -> np.ndarray:
    """Formal charge of each atom: as a file or SMILES stated it, otherwise worked out from its bonds.

    Worked out, nitrogen and phosphorus with four bonds are +1, nitrogen
    with two bonds -1; oxygen, sulfur and selenium with one bond -1 and with
    three bonds +1; boron with four bonds -1; isolated alkali, alkaline
    earth and zinc atoms are cations and isolated halogens anions. Aromatic
    rings are kekulized first, so ring atoms are judged by single and double
    bonds. Carbon is always neutral. Molecules read from a file without
    their hydrogens are left neutral apart from isolated ions, since their
    bonds cannot tell a hydroxyl from an anion.

    These are the charges every writer and Open Babel tool is given (see
    :func:`miles.core.chem.to_openbabel`).

    Parameters
    ----------
    mol : Molecule
        Structure, best with all its hydrogens.

    Returns
    -------
    numpy.ndarray
        Integer charges, shape ``(n_atoms,)``.
    """
    obmol = to_openbabel(mol)
    return np.array([a.GetFormalCharge() for a in ob.OBMolAtomIter(obmol)], dtype=int).reshape(-1)


def partial_charges(mol: Molecule, model: ChargeModel) -> np.ndarray:
    """Compute partial atomic charges.

    Parameters
    ----------
    mol : Molecule
        Structure with all its hydrogens (not modified).
    model : ChargeModel
        ``GASTEIGER`` or ``MMFF94``. The Gasteiger model of Open Babel
        leaves the formal charges out, so they are added to the atoms that
        carry them; both models then sum to the net charge.

    Returns
    -------
    numpy.ndarray
        Charges in elementary charges, shape ``(n_atoms,)``.

    Raises
    ------
    ValueError
        For ``ChargeModel.FILE``, whose charges are read, not computed.
    ChemError
        If the model has no parameters for some element of the structure.
    """
    if model == ChargeModel.FILE:
        raise ValueError("Charges from a file are read with the structure, not computed")
    if not len(mol):
        return np.zeros(0)
    obmol = to_openbabel(mol)   # with the formal charges, see formal_charges()
    formal = np.array([a.GetFormalCharge() for a in ob.OBMolAtomIter(obmol)], dtype=float)
    charge_model = ob.OBChargeModel.FindType(str(model))
    if charge_model is None:
        raise ChemError(f"Open Babel has no {model.label} charge model")
    with silence_openbabel():
        ok = charge_model.ComputeCharges(obmol)
    charges = np.array([a.GetPartialCharge() for a in ob.OBMolAtomIter(obmol)], dtype=float)
    if not ok or not np.all(np.isfinite(charges)):
        elements = ", ".join(sorted(set(mol.symbols) - {"H", "C", "N", "O"})) or "this structure"
        raise ChemError(f"{model.label} charges are not defined for {elements}")
    if model == ChargeModel.GASTEIGER:
        charges += formal
    return charges


def dipole_moment(mol: Molecule, charges: np.ndarray) -> np.ndarray:
    """Dipole moment of point charges on the atoms.

    The moment is taken about the centre of mass. For a neutral molecule the
    choice of origin does not matter; for an ion it does, and the result
    describes the charge distribution about the centre of mass.

    Parameters
    ----------
    mol : Molecule
        Structure supplying the positions (angstrom) and masses.
    charges : numpy.ndarray
        Charge of each atom in elementary charges.

    Returns
    -------
    numpy.ndarray
        Dipole vector in debye, pointing from the negative towards the
        positive charge, shape ``(3,)``.
    """
    if not len(mol):
        return np.zeros(3)
    q = np.asarray(charges, dtype=float)
    return (q[:, None] * (mol.coords - mol.center_of_mass())).sum(axis=0) * DEBYE_PER_E_ANGSTROM
