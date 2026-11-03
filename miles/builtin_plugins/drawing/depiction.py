"""2D structure drawings (skeletal formulae) as SVG, made with Open Babel.

Open Babel lays the molecule out in two dimensions and draws it. Before
that, the hydrogens on carbon are made implicit, except one on each
stereocentre, which Open Babel needs to draw the wedge or hash bond that
shows the configuration taken from the 3D geometry. Formal charges are
worked out from the bonding, so charged groups are drawn as such.

Open Babel writes a fixed-size picture with the drawing nested inside;
:func:`depiction_svg` returns the inner drawing as a standalone SVG at its
natural size, which also suits renderers that do not allow nested SVG
elements, such as Qt's.
"""

import xml.etree.ElementTree as ET

from openbabel import openbabel as ob

from miles import ChemError, Molecule, formal_charges, silence_openbabel, to_openbabel
from miles.builtin_plugins.shared.molecules import MAX_HEAVY_ATOMS, hydrogens_missing

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
_STYLE = ("font-family", "stroke", "stroke-width", "stroke-linecap")


def _layout(mol: Molecule) -> "ob.OBMol":
    """Open Babel molecule with 2D coordinates and implicit hydrogens, ready to draw."""
    obmol = to_openbabel(mol)
    if not hydrogens_missing(mol):
        for atom, charge in zip(ob.OBMolAtomIter(obmol), formal_charges(mol), strict=True):
            atom.SetFormalCharge(int(charge))
    ob.StereoFrom3D(obmol)
    facade = ob.OBStereoFacade(obmol)
    centres = {a.GetId() for a in ob.OBMolAtomIter(obmol) if facade.HasTetrahedralStereo(a.GetId())}
    kept: set[int] = set()
    implicit = []
    for atom in ob.OBMolAtomIter(obmol):
        if atom.GetAtomicNum() != 1:
            continue
        neighbours = list(ob.OBAtomAtomIter(atom))
        if len(neighbours) == 1 and neighbours[0].GetId() in centres and neighbours[0].GetId() not in kept:
            kept.add(neighbours[0].GetId())   # one explicit H per stereocentre carries the wedge
            continue
        if len(neighbours) == 1 and neighbours[0].GetAtomicNum() != 1:
            implicit.append(atom)
    for atom in implicit:
        heavy = next(ob.OBAtomAtomIter(atom))
        heavy.SetImplicitHCount(heavy.GetImplicitHCount() + 1)
        obmol.DeleteHydrogen(atom)
    with silence_openbabel():
        ob.OBOp.FindType("gen2D").Do(obmol)
    return obmol


def depiction_svg(mol: Molecule, title: str = "") -> str:
    """Draw a structure as a 2D skeletal formula.

    Parameters
    ----------
    mol : Molecule
        Structure to draw; it is not changed. Its 3D geometry sets the
        stereochemistry shown.
    title : str, optional
        Text of the SVG ``<title>`` element (not drawn).

    Returns
    -------
    str
        A standalone SVG document with a transparent background, sized to
        the drawing (one unit per pixel).

    Raises
    ------
    ChemError
        If the structure is empty or has more than
        :data:`~miles.builtin_plugins.shared.molecules.MAX_HEAVY_ATOMS` heavy atoms.
    """
    heavy = sum(a.symbol != "H" for a in mol.atoms)
    if not len(mol):
        raise ChemError("No atoms to draw")
    if heavy > MAX_HEAVY_ATOMS:
        raise ChemError(f"2D drawings are for molecules of up to {MAX_HEAVY_ATOMS} heavy atoms; "
                        "select a ligand to draw it")
    obmol = _layout(mol)
    conv = ob.OBConversion()
    conv.SetOutFormat("svg")
    conv.AddOption("d", ob.OBConversion.OUTOPTIONS)            # no molecule name
    conv.AddOption("C", ob.OBConversion.OUTOPTIONS)            # skeletal ends: no CH3 labels
    conv.AddOption("b", ob.OBConversion.OUTOPTIONS, "none")    # transparent background
    with silence_openbabel():
        text = conv.WriteString(obmol)
    return _standalone(text, title)


def _standalone(text: str, title: str) -> str:
    """Lift the drawing out of Open Babel's fixed-size outer SVG element."""
    ET.register_namespace("", SVG_NAMESPACE)
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ChemError(f"Open Babel wrote an unreadable drawing: {exc}") from exc
    inner = root.find(f".//{{{SVG_NAMESPACE}}}svg")
    if inner is None or not inner.get("viewBox"):
        raise ChemError("Open Babel wrote no drawing")
    view = inner.get("viewBox", "")
    _, _, width, height = view.split()
    out = ET.Element(f"{{{SVG_NAMESPACE}}}svg", {"version": "1.1", "width": width, "height": height,
                                                 "viewBox": view})
    for key in _STYLE:
        if inner.get(key):
            out.set(key, inner.get(key, ""))
    ET.SubElement(out, f"{{{SVG_NAMESPACE}}}title").text = title
    out.extend(child for child in inner if child.tag != f"{{{SVG_NAMESPACE}}}title")
    return '<?xml version="1.0"?>\n' + ET.tostring(out, encoding="unicode") + "\n"
