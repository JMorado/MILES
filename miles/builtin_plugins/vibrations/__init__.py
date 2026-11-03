"""Vibrations: the IR spectrum and normal modes of a frequency calculation, animated or as arrows.

The panel opens by itself for a structure with vibrational modes and from
Tools > Vibrations. Choosing a mode plays it in the timeline, Arrows draws
how each atom moves, and the right-click menu of the spectrum copies or
exports it.
"""

import numpy as np

from miles.builtin_plugins.vibrations.panel import VibrationsPanel
from miles.plugin import DocumentView, DrawContext, PluginContext, SceneItems
from miles.plugin.render import arrow_rows

NAME = "Vibrations"
VERSION = "1.0"
DESCRIPTION = "IR spectrum and normal modes of frequency calculations, animated or drawn as arrows"

ARROW_RADIUS = 0.055   # angstrom, shafts of the arrows of a mode


def register(context: PluginContext) -> None:
    """Add the panel, its arrows and Tools > Vibrations, which opens and closes it."""
    app = context.app
    panel = VibrationsPanel()
    context.add_document_panel("vibrations", panel, "Tools > Vibrations", shortcut="Ctrl+Shift+V",
                               before="protein tools")
    panel.exportRequested.connect(lambda kind: app.export_plot(panel, kind))
    panel.playRequested.connect(lambda: app.current.play() if app.current is not None else None)

    def arrows(doc: DocumentView, drawing: DrawContext) -> SceneItems | None:
        """Arrows from the atoms as shown along the displacements of the chosen mode."""
        vectors = panel.arrows_for(doc)
        mol = doc.molecule
        if vectors is None or len(vectors) != len(mol) or not len(mol):
            return None
        moving = doc.shown_atoms() & (np.linalg.norm(vectors, axis=1) > 0)
        if not moving.any():
            return None
        start = mol.coords[moving]
        tips = start + vectors[moving]
        return SceneItems(arrows=arrow_rows(start, tips, ARROW_RADIUS, drawing.theme.accent))

    layer = context.add_scene_items(arrows)
    panel.arrowsChanged.connect(lambda: layer.update())
