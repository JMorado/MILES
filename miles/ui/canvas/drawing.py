"""What a canvas draws: the arrays of a frame, what plugins draw with the molecule, and the export overlay."""

from dataclasses import fields
from typing import TYPE_CHECKING, Any

import numpy as np
from PySide6.QtGui import QImage, QPainter

from miles.core.neighbors import cross_pairs
from miles.render.arrows import COLUMNS
from miles.render.layers import Annotations, DrawContext, Label, Line, MeshItem, SceneItems
from miles.render.picking import AtomHit
from miles.render.scene import Highlight
from miles.render.theme import Theme
from miles.ui.canvas.extensions import Layer
from miles.ui.canvas.overlay import paint_annotations
from miles.ui.canvas.scene_cache import Frame

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

RIBBON_ATOMS = frozenset({"N", "CA", "C", "O", "OXT", "P", "OP1", "OP2", "O1P", "O2P", "O5'", "C5'", "C4'",
                          "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"})
"""Backbone atoms that a cartoon replaces; annotation lines to them can end on the ribbon."""

RIBBON_REACH = 4.0          # angstrom; how far from such an atom the ribbon is looked for

RIBBON_MIN_LENGTH = 0.3     # angstrom; shorter lines to the ribbon are left out

RESOLVED_LAYERS = 32        # annotation layers whose resolved lines are kept between paints

ResolvedLine = tuple[Line, np.ndarray, np.ndarray]
"""A plugin's line with its two ends in angstrom."""


class Drawing:
    """What one canvas draws: the frame, and what plugins draw with the molecule.

    Plugin layers are called with the canvas and their results cached
    until the document, the display options, the theme or the layer
    changes; exported pictures always call them again.

    Parameters
    ----------
    canvas : Canvas
        The canvas whose drawing this is.
    """

    def __init__(self, canvas: "Canvas") -> None:
        self.canvas = canvas
        self._layer_cache: dict[Layer, tuple[Any, Any]] = {}
        # resolved lines of each annotation layer: id -> (annotations, cartoon mesh, lines)
        self._resolved_lines: dict[int, tuple[Annotations, Any, list[ResolvedLine]]] = {}

    def forget_layers(self) -> None:
        """Call every plugin layer again before the next frame (e.g. the charges changed)."""
        self._layer_cache.clear()

    def hover_highlights(self) -> dict[int, Highlight]:
        """Highlight of the atom or bond under the cursor, as the current tool shows it."""
        canvas = self.canvas
        hit = canvas.hover
        if hit is None:
            return {}
        highlight = canvas.tool.highlight(canvas, hit)
        return dict.fromkeys([hit.index] if isinstance(hit, AtomHit) else [hit.i, hit.j], highlight)

    def highlights(self) -> dict[int, Highlight]:
        """Highlight of each selected or hovered atom (a selected atom stays selected under the cursor)."""
        result = dict.fromkeys(self.canvas.selection, Highlight.SELECTED)
        for atom, highlight in self.hover_highlights().items():
            result.setdefault(atom, highlight)
        return result

    def frame(self) -> Frame:
        """Arrays to draw now."""
        canvas = self.canvas
        doc = canvas.doc
        items = self.scene_items(dragging=canvas.dragging())
        return canvas.scenes.frame(doc.mol, doc.revision, canvas.display, canvas.theme,
                                   self.hover_highlights(), canvas.selection, atom_colors=items.atom_colors,
                                   arrows=items.arrows, meshes=items.meshes, display=doc.display,
                                   layout_revision=doc.layout_revision)

    def layer_results(self, layers: list[Layer], context: DrawContext) -> list[Any]:
        """Results of plugin layers for this canvas, leaving out the empty ones.

        On screen each result is kept until the document, the display
        options, the theme or the layer changes; exported pictures always
        call the layers again.
        """
        canvas = self.canvas
        results = []
        style = tuple(getattr(canvas.display, f.name) for f in fields(canvas.display))
        for layer in layers:
            if context.export:
                result = layer.function(canvas, context)
            else:
                key = (id(canvas.mol), canvas.doc.revision, style, context.theme.name, layer.generation,
                       context.dragging)
                cached = self._layer_cache.get(layer)
                if cached is None or cached[0] != key:
                    cached = (key, layer.function(canvas, context))
                    self._layer_cache[layer] = cached
                result = cached[1]
            if result is not None:
                results.append(result)
        return results

    def scene_items(self, theme: Theme | None = None, export: bool = False,
                    dragging: bool = False) -> SceneItems:
        """Arrows, meshes and atom colours that plugins draw in 3D, merged.

        Parameters
        ----------
        theme : Theme, optional
            Colours (default: the canvas theme).
        export : bool, optional
            For an exported picture or movie frame.
        dragging : bool, optional
            Atoms are being dragged (see :class:`~miles.render.layers.DrawContext`).

        Returns
        -------
        SceneItems
            Every layer's arrows and meshes, and the atom colours of the last
            layer that gives colours for all atoms.
        """
        context = DrawContext(theme or self.canvas.theme, export, dragging)
        items = self.layer_results(self.canvas.extensions.scene_items, context)
        arrows = [np.asarray(item.arrows, np.float32) for item in items if len(item.arrows)]
        colors = next((item.atom_colors for item in reversed(items)
                       if item.atom_colors is not None and len(item.atom_colors) == len(self.canvas.mol)),
                      None)
        meshes: list[MeshItem] = [mesh for item in items for mesh in item.meshes]
        return SceneItems(np.concatenate(arrows) if arrows else np.zeros((0, COLUMNS), np.float32), meshes,
                          colors)

    def annotations(self, theme: Theme | None = None, export: bool = False) -> list[Annotations]:
        """Lines and labels that plugins draw over the molecule, one entry per layer."""
        context = DrawContext(theme or self.canvas.theme, export)
        return self.layer_results(self.canvas.extensions.annotations, context)

    def resolve_lines(self, annotations: Annotations) -> list[tuple[Line, np.ndarray, np.ndarray]]:
        """Return the lines of a layer with their two ends in angstrom, leaving out those not drawn.

        Ends anchored to atoms take the atoms' coordinates. A line with
        ``drawn_only`` is left out when an atom is not shown; one with
        ``ribbon`` moves an end at a backbone atom drawn only as part of a
        cartoon (:data:`RIBBON_ATOMS`) to the nearest point of the ribbon,
        and is left out when an end is any other atom that is not drawn or
        when it would be shorter than :data:`RIBBON_MIN_LENGTH`.
        """
        canvas = self.canvas
        mesh = canvas.scenes.base.mesh if canvas.scenes.base is not None else None
        cache = self._resolved_lines
        cached = cache.get(id(annotations))
        if cached is not None and cached[0] is annotations and cached[1] is mesh:
            return cached[2]
        mol = canvas.mol
        n = len(mol)
        drawn = canvas.drawn_atoms() if any(line.ribbon for line in annotations.lines) else None
        shown = canvas.shown_atoms() if any(line.drawn_only for line in annotations.lines) else None
        ribbon = mesh is not None and not mesh.empty
        resolved: list[list[Any]] = []
        queries: list[tuple[int, int]] = []
        for line in annotations.lines:
            ends: list[np.ndarray] = []
            for side, anchor in enumerate((line.start, line.end)):
                if isinstance(anchor, int | np.integer):
                    i = int(anchor)
                    if not 0 <= i < n:
                        break
                    if line.ribbon and drawn is not None and not drawn[i]:
                        if not (ribbon and mol.atoms[i].name in RIBBON_ATOMS):
                            break
                        queries.append((len(resolved), side))
                    elif line.drawn_only and shown is not None and not shown[i]:
                        break
                    ends.append(mol.coords[i].astype(float))
                else:
                    ends.append(np.asarray(anchor, dtype=float).reshape(3))
            if len(ends) == 2:
                resolved.append([line, ends[0], ends[1]])
        if queries and mesh is not None:
            starts = np.array([resolved[k][side + 1] for k, side in queries])
            nearest = _nearest_vertices(starts, mesh.vertices)
            for (k, side), point in zip(queries, nearest, strict=True):
                if np.isnan(point).any():
                    resolved[k][0] = None
                else:
                    resolved[k][side + 1] = point
        result = [(line, a, b) for line, a, b in resolved if line is not None
                  and (not line.ribbon or float(np.linalg.norm(a - b)) > RIBBON_MIN_LENGTH)]
        if len(cache) >= RESOLVED_LAYERS:
            cache.clear()
        cache[id(annotations)] = (annotations, mesh, result)
        return result

    def resolve_label(self, label: Label) -> np.ndarray | None:
        """Point in angstrom of a label's anchor, or None if the label is not drawn."""
        return self.anchor_point(label.at, label.drawn_only)

    def anchor_point(self, anchor: Any, drawn_only: bool = False) -> np.ndarray | None:
        """Point in angstrom of an atom index or a point, or None for an atom that is missing or not shown."""
        canvas = self.canvas
        if isinstance(anchor, int | np.integer):
            i = int(anchor)
            if not 0 <= i < len(canvas.mol) or (drawn_only and not canvas.shown_atoms()[i]):
                return None
            return canvas.mol.coords[i].astype(float)
        return np.asarray(anchor, dtype=float).reshape(3)

    def paint_export_overlay(self, image: QImage, theme: Theme) -> None:
        """Draw the annotations (measurements, stereo labels, what plugins add) onto an exported picture."""
        canvas = self.canvas
        annotations = self.annotations(theme, export=True)
        if not (canvas.mol.measurements or canvas.display.stereo or annotations):
            return
        p = QPainter(image)
        p.scale(image.width() / max(canvas.width(), 1), image.height() / max(canvas.height(), 1))
        paint_annotations(p, canvas, theme, annotations)
        p.end()


def _nearest_vertices(queries: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """Nearest vertex to each query point within :data:`RIBBON_REACH`; NaN where there is none."""
    out = np.full((len(queries), 3), np.nan)
    i, j = cross_pairs(queries, vertices, RIBBON_REACH)
    if len(i):
        distance = np.linalg.norm(queries[i] - vertices[j], axis=1)
        order = np.lexsort((distance, i))
        first = np.unique(i[order], return_index=True)[1]
        out[i[order][first]] = vertices[j[order][first]]
    return out
