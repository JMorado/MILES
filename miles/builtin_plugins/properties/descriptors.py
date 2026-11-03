"""Molecular properties: composition, identifiers and drug-likeness descriptors.

The descriptors are Open Babel's: logP and molar refractivity from the atom
contributions of Wildman and Crippen, the topological polar surface area
(TPSA) from Ertl's fragment contributions, and rotatable bonds from Open
Babel's rotor rules. Hydrogen-bond donors and acceptors are counted as in
Lipinski's rule of five: donors are nitrogen and oxygen atoms that carry
hydrogen, acceptors are all nitrogen and oxygen atoms.

Identifiers and descriptors are computed only for small molecules with
all their hydrogens. Proteins, nucleic acids and other large structures get
their composition only; so do structures with missing hydrogens, whose
SMILES and charges would describe radicals.
"""

from dataclasses import dataclass

from openbabel import openbabel as ob

from miles import Molecule, element, formal_charges, silence_openbabel, to_openbabel
from miles.builtin_plugins.shared.molecules import hydrogens_missing, is_small_molecule

LIPINSKI_LIMITS = (("molecular weight", 500.0), ("logP", 5.0), ("H-bond donors", 5), ("H-bond acceptors", 10))
"""Lipinski's rule of five: an orally active drug usually exceeds at most one of these."""


@dataclass
class Properties:
    """Properties of a structure.

    Attributes
    ----------
    formula : str
        Molecular formula in Hill order.
    weight : float
        Molecular weight in g/mol, from standard atomic weights.
    exact_mass : float
        Monoisotopic mass in g/mol (most abundant isotope of each element).
    atoms : int
        Number of atoms.
    heavy_atoms : int
        Number of atoms other than hydrogen.
    charge : int or None
        Net formal charge worked out from the bonding; None when hydrogens
        are missing.
    smiles, inchi, inchikey : str or None
        Canonical SMILES (with stereochemistry from the geometry), standard
        InChI and InChIKey; None when not computed.
    logp : float or None
        Octanol-water partition coefficient (Wildman-Crippen).
    tpsa : float or None
        Topological polar surface area in square angstrom (Ertl).
    molar_refractivity : float or None
        Molar refractivity (Wildman-Crippen).
    donors, acceptors : int or None
        Hydrogen-bond donors and acceptors as counted for the rule of five.
    rotatable_bonds : int or None
        Number of rotatable bonds.
    small_molecule : bool
        False for polymers and large structures, which get no identifiers
        or descriptors.
    hydrogens_missing : bool
        Some carbon atoms have fewer bonds than they need, so hydrogens are
        probably missing; identifiers, charge and descriptors are then left
        out.
    """

    formula: str
    weight: float
    exact_mass: float
    atoms: int
    heavy_atoms: int
    charge: int | None = None
    smiles: str | None = None
    inchi: str | None = None
    inchikey: str | None = None
    logp: float | None = None
    tpsa: float | None = None
    molar_refractivity: float | None = None
    donors: int | None = None
    acceptors: int | None = None
    rotatable_bonds: int | None = None
    small_molecule: bool = True
    hydrogens_missing: bool = False

    @property
    def lipinski_violations(self) -> list[str] | None:
        """Names of the rule-of-five limits exceeded, or None without descriptors."""
        values = (self.weight, self.logp, self.donors, self.acceptors)
        if any(v is None for v in values):
            return None
        return [name for (name, limit), value in zip(LIPINSKI_LIMITS, values, strict=True)
                if value is not None and value > limit]


def _exact_mass(mol: Molecule) -> float:
    """Sum of the masses of each element's most abundant isotope."""
    masses: dict[str, float] = {}
    for symbol in mol.symbols:
        if symbol not in masses:
            masses[symbol] = ob.GetExactMass(element(symbol).number) if element(symbol).number else 0.0
    return float(sum(masses[s] for s in mol.symbols))


def _descriptor(obmol: "ob.OBMol", name: str) -> float:
    with silence_openbabel():
        return float(ob.OBDescriptor.FindType(name).Predict(obmol))


def _write(obmol: "ob.OBMol", fmt: str) -> str:
    conv = ob.OBConversion()
    conv.SetOutFormat(fmt)
    with silence_openbabel():
        text = conv.WriteString(obmol)
    return text.split()[0] if text.strip() else ""


def molecular_properties(mol: Molecule) -> Properties:
    """Compute the properties of a structure.

    Parameters
    ----------
    mol : Molecule
        Structure, with all its hydrogens for identifiers and descriptors.

    Returns
    -------
    Properties
        Composition always; charge, identifiers and descriptors for small
        molecules with complete hydrogens.
    """
    props = Properties(formula=mol.formula(), weight=float(mol.masses().sum()) if len(mol) else 0.0,
                       exact_mass=_exact_mass(mol), atoms=len(mol),
                       heavy_atoms=int((mol.numbers != 1).sum()),
                       small_molecule=is_small_molecule(mol), hydrogens_missing=hydrogens_missing(mol))
    if not len(mol) or props.hydrogens_missing:
        return props
    charges = formal_charges(mol)
    props.charge = int(charges.sum())
    if not props.small_molecule:
        return props
    obmol = to_openbabel(mol)
    for atom, charge in zip(ob.OBMolAtomIter(obmol), charges, strict=True):
        atom.SetFormalCharge(int(charge))
    ob.StereoFrom3D(obmol)
    props.smiles = _write(obmol, "can") or None
    props.inchi = _write(obmol, "inchi") or None
    props.inchikey = _write(obmol, "inchikey") or None
    props.logp = _descriptor(obmol, "logP")
    props.tpsa = _descriptor(obmol, "TPSA")
    props.molar_refractivity = _descriptor(obmol, "MR")
    props.donors = round(_descriptor(obmol, "HBD"))
    props.acceptors = round(_descriptor(obmol, "HBA2"))
    props.rotatable_bonds = round(_descriptor(obmol, "rotors"))
    return props
