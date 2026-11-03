"""MILES: a molecule builder and 3D viewer for proteins, ligands, binding pockets and small molecules.

``import miles`` gives the chemistry library that scripts and plugins
compute with; :mod:`miles.plugin` gives what a plugin needs to work with the
window. These two are the public, stable surface, versioned together by
``REQUIRES_API``; every other module is internal and may change.

Names are loaded when first used, so ``import miles`` is quick and imports
neither Qt nor Open Babel; :func:`launch` opens the window only when called.

Examples
--------
>>> import miles
>>> mol, info = miles.read_file("4hhb.pdb")
>>> pocket = miles.find(mol, "byres within 5 of ligand")
>>> app = miles.launch("ligand.sdf")      # the window, from IPython or Jupyter

What it provides, by purpose:

Structures and files
    ``Molecule``, ``Atom``, ``AtomLike``, ``Representation``, ``Geometry``,
    ``JobInfo``, ``FileFormatError``, ``read_file``, ``write_file``,
    ``parse_text``, ``read_text``, ``write_cube``
Selections
    ``find``, ``SelectionError``, ``residue_ids``
Geometry, neighbours and trajectories
    ``distance``, ``angle``, ``dihedral``, ``kabsch``, ``rotation_matrix``,
    ``unit_vector``, ``place_atom``, ``set_bond_length``, ``set_angle``,
    ``set_dihedral``, ``pairs_within``, ``cross_pairs``, ``near_mask``,
    ``align_frames``, ``rmsd_atoms``, ``rmsd_series``
Elements
    ``element``, ``normalize_symbol``
Open Babel, SMILES and force fields
    ``to_openbabel``, ``from_openbabel``, ``openbabel_coordinates``,
    ``silence_openbabel``, ``from_smiles``, ``to_smiles``, ``is_smiles``,
    ``ChemError``, ``FORCE_FIELDS``, ``minimize``, ``conformers``,
    ``clean_geometry``, ``clean_selection``
Building and perception
    ``add_hydrogens``, ``clash_score``, ``amino_acid_codes``,
    ``amino_acid_template``, ``rebond``, ``infer_atom_types``,
    ``infer_geometry``
Proteins and nucleic acids
    ``AMINO_ACIDS``, ``NUCLEOTIDES``, ``ONE_LETTER``, ``residues``, ``is_polymer``,
    ``has_polymer``, ``is_water``, ``water_mask``
Quantum chemistry data and charges
    ``VibrationalMode``, ``OrbitalSet``, ``Spin``, ``Volume``,
    ``mapped_isosurface``, ``ChargeModel``, ``dipole_moment``,
    ``formal_charges``
Units
    ``BOHR``, ``HARTREE``, ``HARTREE_EV``, ``ELECTRONVOLT``, ``DEBYE_PER_AU``
"""

import importlib
import os
import sys
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from PySide6.QtCore import QSettings

    from miles.api import App
    from miles.core.building import (
        add_hydrogens as add_hydrogens,
    )
    from miles.core.building import (
        clash_score as clash_score,
    )
    from miles.core.charges import (
        ChargeModel as ChargeModel,
    )
    from miles.core.charges import (
        dipole_moment as dipole_moment,
    )
    from miles.core.charges import (
        formal_charges as formal_charges,
    )
    from miles.core.chem import (
        FORCE_FIELDS as FORCE_FIELDS,
    )
    from miles.core.chem import (
        ChemError as ChemError,
    )
    from miles.core.chem import (
        conformers as conformers,
    )
    from miles.core.chem import (
        from_openbabel as from_openbabel,
    )
    from miles.core.chem import (
        from_smiles as from_smiles,
    )
    from miles.core.chem import (
        is_smiles as is_smiles,
    )
    from miles.core.chem import (
        minimize as minimize,
    )
    from miles.core.chem import (
        openbabel_coordinates as openbabel_coordinates,
    )
    from miles.core.chem import (
        silence_openbabel as silence_openbabel,
    )
    from miles.core.chem import (
        to_openbabel as to_openbabel,
    )
    from miles.core.chem import (
        to_smiles as to_smiles,
    )
    from miles.core.elements import (
        element as element,
    )
    from miles.core.elements import (
        normalize_symbol as normalize_symbol,
    )
    from miles.core.forcefield import (
        clean_geometry as clean_geometry,
    )
    from miles.core.forcefield import (
        clean_selection as clean_selection,
    )
    from miles.core.fragments import (
        amino_acid_codes as amino_acid_codes,
    )
    from miles.core.fragments import (
        amino_acid_template as amino_acid_template,
    )
    from miles.core.geometry import (
        Geometry as Geometry,
    )
    from miles.core.geometry import (
        angle as angle,
    )
    from miles.core.geometry import (
        dihedral as dihedral,
    )
    from miles.core.geometry import (
        distance as distance,
    )
    from miles.core.geometry import (
        kabsch as kabsch,
    )
    from miles.core.geometry import (
        place_atom as place_atom,
    )
    from miles.core.geometry import (
        rotation_matrix as rotation_matrix,
    )
    from miles.core.geometry import (
        unit_vector as unit_vector,
    )
    from miles.core.internal import (
        set_angle as set_angle,
    )
    from miles.core.internal import (
        set_bond_length as set_bond_length,
    )
    from miles.core.internal import (
        set_dihedral as set_dihedral,
    )
    from miles.core.jobinfo import (
        JobInfo as JobInfo,
    )
    from miles.core.molecule import (
        Atom as Atom,
    )
    from miles.core.molecule import (
        AtomLike as AtomLike,
    )
    from miles.core.molecule import (
        Molecule as Molecule,
    )
    from miles.core.molecule import (
        Representation as Representation,
    )
    from miles.core.neighbors import (
        cross_pairs as cross_pairs,
    )
    from miles.core.neighbors import (
        near_mask as near_mask,
    )
    from miles.core.neighbors import (
        pairs_within as pairs_within,
    )
    from miles.core.orbitals import (
        OrbitalSet as OrbitalSet,
    )
    from miles.core.orbitals import (
        Spin as Spin,
    )
    from miles.core.perception import (
        infer_atom_types as infer_atom_types,
    )
    from miles.core.perception import (
        infer_geometry as infer_geometry,
    )
    from miles.core.perception import (
        rebond as rebond,
    )
    from miles.core.protein import (
        AMINO_ACIDS as AMINO_ACIDS,
    )
    from miles.core.protein import (
        NUCLEOTIDES as NUCLEOTIDES,
    )
    from miles.core.protein import (
        ONE_LETTER as ONE_LETTER,
    )
    from miles.core.protein import (
        has_polymer as has_polymer,
    )
    from miles.core.protein import (
        is_polymer as is_polymer,
    )
    from miles.core.protein import (
        is_water as is_water,
    )
    from miles.core.protein import (
        residues as residues,
    )
    from miles.core.protein import (
        water_mask as water_mask,
    )
    from miles.core.selection import (
        SelectionError as SelectionError,
    )
    from miles.core.selection import (
        find as find,
    )
    from miles.core.selection import (
        residue_ids as residue_ids,
    )
    from miles.core.trajectory import (
        align_frames as align_frames,
    )
    from miles.core.trajectory import (
        rmsd_atoms as rmsd_atoms,
    )
    from miles.core.trajectory import (
        rmsd_series as rmsd_series,
    )
    from miles.core.units import (
        BOHR as BOHR,
    )
    from miles.core.units import (
        DEBYE_PER_AU as DEBYE_PER_AU,
    )
    from miles.core.units import (
        ELECTRONVOLT as ELECTRONVOLT,
    )
    from miles.core.units import (
        HARTREE as HARTREE,
    )
    from miles.core.units import (
        HARTREE_EV as HARTREE_EV,
    )
    from miles.core.vibrations import (
        VibrationalMode as VibrationalMode,
    )
    from miles.core.volume import (
        Volume as Volume,
    )
    from miles.core.volume import (
        mapped_isosurface as mapped_isosurface,
    )
    from miles.files.common import (
        FileFormatError as FileFormatError,
    )
    from miles.files.cube import (
        write_cube as write_cube,
    )
    from miles.files.formats import (
        parse_text as parse_text,
    )
    from miles.files.formats import (
        read_file as read_file,
    )
    from miles.files.formats import (
        write_file as write_file,
    )
    from miles.files.openbabel import (
        read_text as read_text,
    )

__version__ = "0.1.0"


def _find_openbabel_data() -> None:
    """Point Open Babel at its data files when the environment does not.

    Open Babel reads force field parameters and descriptor tables from the
    folder in ``BABEL_DATADIR``, or else from the path built into it. On
    Windows that path is not updated when conda installs the package, and
    conda sets ``BABEL_DATADIR`` only when an environment is activated, so
    Python started directly (by a Start menu shortcut, or in a shell that
    does not run the activation scripts) would find no data. The data is in
    ``share/openbabel`` of the environment, or in a folder named after the
    version inside it.
    """
    if os.environ.get("BABEL_DATADIR"):
        return
    root = os.path.join(sys.prefix, "share", "openbabel")
    if not os.path.isdir(root):
        return
    for folder in [root, *sorted(os.path.join(root, name) for name in os.listdir(root))]:
        if os.path.isfile(os.path.join(folder, "mmff94.ff")):
            os.environ["BABEL_DATADIR"] = folder
            return


_find_openbabel_data()

_LAZY = {
    # structures and files
    "Molecule": "miles.core.molecule",
    "Atom": "miles.core.molecule",
    "AtomLike": "miles.core.molecule",
    "Representation": "miles.core.molecule",
    "Geometry": "miles.core.geometry",
    "JobInfo": "miles.core.jobinfo",
    "FileFormatError": "miles.files.common",
    "read_file": "miles.files.formats",
    "write_file": "miles.files.formats",
    "parse_text": "miles.files.formats",
    "read_text": "miles.files.openbabel",
    "write_cube": "miles.files.cube",
    # selections
    "find": "miles.core.selection",
    "SelectionError": "miles.core.selection",
    "residue_ids": "miles.core.selection",
    # geometry, neighbours and trajectories
    "distance": "miles.core.geometry",
    "angle": "miles.core.geometry",
    "dihedral": "miles.core.geometry",
    "kabsch": "miles.core.geometry",
    "rotation_matrix": "miles.core.geometry",
    "unit_vector": "miles.core.geometry",
    "place_atom": "miles.core.geometry",
    "set_bond_length": "miles.core.internal",
    "set_angle": "miles.core.internal",
    "set_dihedral": "miles.core.internal",
    "pairs_within": "miles.core.neighbors",
    "cross_pairs": "miles.core.neighbors",
    "near_mask": "miles.core.neighbors",
    "align_frames": "miles.core.trajectory",
    "rmsd_atoms": "miles.core.trajectory",
    "rmsd_series": "miles.core.trajectory",
    # elements
    "element": "miles.core.elements",
    "normalize_symbol": "miles.core.elements",
    # open babel, smiles and force fields
    "to_openbabel": "miles.core.chem",
    "from_openbabel": "miles.core.chem",
    "openbabel_coordinates": "miles.core.chem",
    "silence_openbabel": "miles.core.chem",
    "from_smiles": "miles.core.chem",
    "to_smiles": "miles.core.chem",
    "is_smiles": "miles.core.chem",
    "ChemError": "miles.core.chem",
    "FORCE_FIELDS": "miles.core.chem",
    "minimize": "miles.core.chem",
    "conformers": "miles.core.chem",
    "clean_geometry": "miles.core.forcefield",
    "clean_selection": "miles.core.forcefield",
    # building and perception
    "add_hydrogens": "miles.core.building",
    "clash_score": "miles.core.building",
    "amino_acid_codes": "miles.core.fragments",
    "amino_acid_template": "miles.core.fragments",
    "rebond": "miles.core.perception",
    "infer_atom_types": "miles.core.perception",
    "infer_geometry": "miles.core.perception",
    # proteins and nucleic acids
    "AMINO_ACIDS": "miles.core.protein",
    "NUCLEOTIDES": "miles.core.protein",
    "ONE_LETTER": "miles.core.protein",
    "residues": "miles.core.protein",
    "is_polymer": "miles.core.protein",
    "has_polymer": "miles.core.protein",
    "is_water": "miles.core.protein",
    "water_mask": "miles.core.protein",
    # quantum chemistry data and charges
    "VibrationalMode": "miles.core.vibrations",
    "OrbitalSet": "miles.core.orbitals",
    "Spin": "miles.core.orbitals",
    "Volume": "miles.core.volume",
    "mapped_isosurface": "miles.core.volume",
    "ChargeModel": "miles.core.charges",
    "dipole_moment": "miles.core.charges",
    "formal_charges": "miles.core.charges",
    # units
    "BOHR": "miles.core.units",
    "HARTREE": "miles.core.units",
    "HARTREE_EV": "miles.core.units",
    "ELECTRONVOLT": "miles.core.units",
    "DEBYE_PER_AU": "miles.core.units",
}
"""Where each public name lives."""

__all__ = [*_LAZY, "__version__", "launch"]


def __getattr__(name: str) -> Any:
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module 'miles' has no attribute {name!r}")
    value = getattr(importlib.import_module(module), name)
    globals()[name] = value   # later lookups skip this function
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))


def launch(*paths: str, plugins: bool = True, appearance: str | None = None,
           settings: "QSettings | None" = None) -> "App":
    """Open the MILES window from IPython, Jupyter or a script and return its ``app``.

    In IPython and Jupyter the window keeps working between cells; a script
    run with plain ``python`` ends with ``app.wait_until_closed()``. See
    :func:`miles.interactive.launch` for the parameters.

    Examples
    --------
    >>> import miles
    >>> app = miles.launch("ligand.sdf")
    >>> app.current.style(mode="tube")
    >>> app.current.image()          # a picture below the cell
    """
    from miles.interactive import launch as start

    return start(*paths, plugins=plugins, appearance=appearance, settings=settings)
