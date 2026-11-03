"""Selection Info: mass, centre and size of the selected atoms.

An example MILES plugin. Copy this file into the plugin folder
(Help > Plugins... > Open Plugin Folder) and restart MILES. It adds:

- a Selection Info panel (View > Selection Info) with the formula, mass,
  centre of mass and radius of gyration of the selected atoms, or of the
  whole structure when nothing is selected, kept up to date as the
  selection and the structure change;
- Tools > Selection Info > Select Heavy Atoms and Tools > Selection Info >
  Move Selection to Origin, two commands that work through the scripting
  interface and can be undone;
- reading and writing ``.atoms`` files, plain text with one
  ``symbol x y z`` line per atom (lines starting with ``#`` are comments).
"""

import numpy as np
from PySide6.QtWidgets import QLabel, QWidget

from miles import Atom, FileFormatError, JobInfo, Molecule, normalize_symbol, rebond
from miles.plugin import App, DocumentView, PluginContext

NAME = "Selection Info"
VERSION = "1.0"
DESCRIPTION = "Mass, centre and radius of gyration of the selection, and .atoms files"


def summary(doc: DocumentView) -> str:
    """Describe the selected atoms, or every atom when nothing is selected."""
    mol = doc.molecule
    if not len(mol):
        return "No atoms"
    atoms = doc.selection or list(range(len(mol)))
    part = mol.subset(atoms)
    masses = part.masses()
    centre = part.center_of_mass()
    radius = float(np.sqrt((masses * ((part.coords - centre) ** 2).sum(axis=1)).sum() / masses.sum()))
    which = f"{len(atoms)} selected atoms" if doc.selection else f"All {len(atoms)} atoms"
    return (f"{which}\n"
            f"Formula   {part.formula()}\n"
            f"Mass      {masses.sum():.3f} g/mol\n"
            f"Centre    {centre[0]:.3f}  {centre[1]:.3f}  {centre[2]:.3f}\n"
            f"Rg        {radius:.3f} A")


def select_heavy_atoms(app: App) -> None:
    """Select every atom that is not a hydrogen."""
    doc = app.current
    if doc is not None:
        atoms = doc.select([i for i, s in enumerate(doc.molecule.symbols) if s != "H"])
        app.message(f"Selected {len(atoms)} heavy atoms")


def move_selection_to_origin(app: App) -> None:
    """Translate the structure so the centre of mass of the selection is at the origin."""
    doc = app.current
    if doc is None or not len(doc.molecule):
        return
    atoms = doc.selection or list(range(len(doc.molecule)))
    with doc.editing() as mol:   # one undoable step
        mol.coords = mol.coords - mol.subset(atoms).center_of_mass()
    app.message("Moved the centre of mass of the selection to the origin")


def read_atoms(path: str) -> tuple[Molecule, JobInfo]:
    """Read a ``.atoms`` file: one ``symbol x y z`` line per atom."""
    atoms, coords = [], []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            parts = line.split("#", 1)[0].split()
            if not parts:
                continue
            symbol = normalize_symbol(parts[0])
            if symbol is None or len(parts) < 4:
                raise FileFormatError(f"Line {number} is not 'symbol x y z': {line.strip()}")
            atoms.append(Atom(symbol))
            coords.append([float(value) for value in parts[1:4]])
    if not atoms:
        raise FileFormatError("The file has no atoms")
    mol = Molecule.from_atoms(atoms, coords)
    rebond(mol)   # the format has no bonds, so they are worked out from the distances
    return mol, JobInfo()


def write_atoms(path: str, mol: Molecule, info: JobInfo) -> None:
    """Write a ``.atoms`` file."""
    with open(path, "w", encoding="utf-8") as handle:
        if info.title:
            handle.write(f"# {info.title}\n")
        for symbol, (x, y, z) in zip(mol.symbols, mol.coords, strict=True):
            handle.write(f"{symbol:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def register(context: PluginContext) -> None:
    """Add the panel, the commands and the file format."""
    text = QLabel()

    def build(app: App) -> QWidget:
        text.setStyleSheet("font-family: Menlo, Consolas, monospace;")
        return text

    def refresh(doc: DocumentView) -> None:
        if doc == context.app.current:
            text.setText(summary(doc))

    panel = context.add_panel("Selection Info", build, shortcut="Ctrl+Alt+I")
    context.on_selection_changed(refresh)
    context.on_structure_changed(refresh)
    context.on_current_changed(refresh)
    if panel is not None and context.app.current is not None:
        refresh(context.app.current)

    def has_atoms(app: App) -> bool:
        return app.current is not None and len(app.current.molecule) > 0

    context.add_command("Tools > Selection Info > Select Heavy Atoms", select_heavy_atoms, enabled=has_atoms)
    context.add_command("Tools > Selection Info > Move Selection to Origin", move_selection_to_origin,
                        enabled=has_atoms)
    context.add_file_reader("atoms", read_atoms)
    context.add_file_writer("atoms", write_atoms)
