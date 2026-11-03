"""Orbitals and Grids: isosurfaces of grids, molecular orbitals and densities, in the grid panel.

The panel opens by itself, above the timeline, for a structure with grids
(cube files) or molecular orbitals (Molden files, GAMESS output). It shows
a grid or an orbital as an isosurface, works out orbitals and densities
when they are first chosen, colours a density by a second grid such as the
electrostatic potential, combines two grids, adds the grids of another
file and saves a grid as a cube file. The isosurface of each tab is saved
in session files.
"""

import os
import re

from miles import FileFormatError, Volume, read_file, write_cube
from miles.builtin_plugins.grids.panel import VolumePanel, request_to_saved
from miles.plugin import App, DocumentView, DrawContext, MeshItem, PluginContext, SceneItems, read_filter
from miles.plugin.render import Mesh, volume_mesh

NAME = "Orbitals and Grids"
VERSION = "1.0"
DESCRIPTION = "Isosurfaces of grids, molecular orbitals and densities, with the grid panel"


def cube_file_name(volume: Volume) -> str:
    """Return a file name for a grid, made from its name, e.g. ``"MO_12_alpha_HOMO.cube"``."""
    stem = re.sub(r"[^A-Za-z0-9+-]+", "_", volume.name).strip("_")
    return f"{stem or 'grid'}.cube"


def add_grid(app: App, panel: VolumePanel, path: str) -> bool:
    """Add the grids of a file to the current document, e.g. a potential for a density.

    The structures must match: same atoms, same coordinates within 0.01
    angstrom. A density on screen is then coloured by the first new grid.

    Returns
    -------
    bool
        Whether any grid was added.
    """
    doc = app.current
    if doc is None:
        return False
    name = os.path.basename(path)
    try:
        mol, info = read_file(path)
    except (OSError, FileFormatError, ValueError, IndexError) as exc:
        app.warn("Could not add grid", f"{name}\n\n{exc}")
        return False
    same = len(mol) == len(doc.molecule)
    if same and len(mol):
        same = float(abs(mol.coords - doc.molecule.coords).max()) < 0.01
    if not info.volumes or not same:
        reason = "has no grids" if not info.volumes else "is for a different structure"
        app.warn("Could not add grid", f"{name} {reason}.")
        return False
    names = {v.name for v in doc.info.volumes}
    first = None
    for volume in info.volumes:
        if volume.name in names:
            volume.name = f"{volume.name} ({name})"
        index = doc.add_volume(volume)   # undoable, like any edit
        first = index if first is None else first
    if first is not None:
        panel.map_grid(first)
    return True


def save_grid(app: App, panel: VolumePanel, path: str) -> bool:
    """Write the grid chosen in the panel, with the current atoms, as a cube file.

    Returns
    -------
    bool
        Whether the file was written.
    """
    current = app.current
    volume = panel.selected_volume()
    if current is None or volume is None:
        return False
    try:
        write_cube(path, current.molecule, volume, current.name)
    except (OSError, ValueError) as exc:
        app.warn("Could not save grid", f"{os.path.basename(path)}\n\n{exc}")
        return False
    app.message(f"Saved {os.path.basename(path)}")
    return True


def register(context: PluginContext) -> None:
    """Add the grid panel and draw the isosurface it chooses for each tab."""
    app = context.app
    panel = VolumePanel()
    context.add_document_panel("grid", panel, place="bottom")

    def isosurface(doc: DocumentView, drawing: DrawContext) -> SceneItems | None:
        shown = panel.volume_shown(doc)
        if shown is None:
            return None
        volume, request = shown
        mapped = panel.mapped_volume(doc)
        wanted = (id(volume), request.level, id(mapped) if mapped is not None else None, request.mapped_range)
        cache = context.document_cache(doc)
        cached: tuple[object, Mesh | None] | None = cache.get("mesh")
        if cached is None or cached[0] != wanted:
            built = volume_mesh(volume, request.level, mapped, request.mapped_range)
            cached = cache["mesh"] = (wanted, built)
        mesh = cached[1]
        return SceneItems(meshes=[MeshItem(mesh, request.alpha)]) if mesh is not None else None

    layer = context.add_scene_items(isosurface)

    def store(doc: DocumentView) -> None:
        """Keep the document's isosurface in its session state and draw it."""
        if doc in app.documents:
            context.document_state(doc)["volume"] = request_to_saved(panel.request(doc))
            layer.update(doc)

    def add_grid_dialog() -> None:
        if app.current is not None:
            path = app.ask_open_path("Add Grid", read_filter())
            if path:
                add_grid(app, panel, path)

    def save_grid_dialog() -> None:
        volume = panel.selected_volume()
        if app.current is not None and volume is not None:
            path = app.ask_save_path("Save Grid", cube_file_name(volume), "Cube files (*.cube *.cub)")
            if path:
                save_grid(app, panel, path)

    panel.volumeChanged.connect(store)
    panel.addGridRequested.connect(add_grid_dialog)
    panel.saveGridRequested.connect(save_grid_dialog)
    for doc in app.documents:
        store(doc)
    context.on_document_opened(store)
