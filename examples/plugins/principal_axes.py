"""Principal Axes: the axes of inertia of a structure, drawn as arrows through its centre of mass.

An example of a plugin that draws. View > Principal Axes switches the
arrows on and off for each tab (and is saved in sessions). Each arrow is
labelled with its moment of inertia, and the status line says what kind of
rotor the structure is. The arrows are 3D scene items and the labels
annotations, so both follow the structure while it is edited or played,
and appear in exported pictures and movies.

Copy this file into the plugin folder (Help > Plugins... > Open Plugin
Folder) and restart MILES.
"""

import numpy as np

from miles import Molecule
from miles.plugin import Annotations, DocumentView, DrawContext, Label, PluginContext, SceneItems
from miles.plugin.render import arrow_rows

NAME = "Principal Axes"
VERSION = "1.0"

COLORS = ["#D9544D", "#4DA35F", "#4A7FD6"]   # smallest, middle and largest moment
RADIUS = 0.05     # angstrom, shaft of the arrows
MARGIN = 0.8      # angstrom the arrows reach beyond the outermost atom
SAME = 0.01       # moments within 1 % of each other count as equal


def principal_axes(mol: Molecule) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the centre of mass, the moments of inertia (amu angstrom^2, ascending) and the axes as rows."""
    centre = mol.center_of_mass()
    r = mol.coords - centre
    masses = mol.masses()
    squared = float(np.sum(masses * np.einsum("ij,ij->i", r, r)))
    tensor = np.eye(3) * squared - np.einsum("i,ij,ik->jk", masses, r, r)
    moments, axes = np.linalg.eigh(tensor)
    return centre, moments, axes.T


def rotor(moments: np.ndarray) -> str:
    """Return the kind of rotor for three moments of inertia in ascending order."""
    a, b, c = (float(m) for m in moments)
    scale = max(c, 1e-9)
    if a / scale < SAME:
        return "linear rotor"
    if (c - a) / scale < SAME:
        return "spherical top"
    if (b - a) / scale < SAME or (c - b) / scale < SAME:
        return "symmetric top"
    return "asymmetric top"


def register(context: PluginContext) -> None:
    """Add the toggle and the arrows, labels and status text it shows."""

    def axes(doc: DocumentView) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
        """Return the centre, the moments and the ends of each arrow, or None when the axes are off."""
        mol = doc.molecule
        if not toggle.value(doc) or len(mol) < 2:
            return None
        centre, moments, directions = principal_axes(mol)
        reach = np.abs((mol.coords - centre) @ directions.T).max(axis=0) + MARGIN
        tails = centre - directions * reach[:, None]
        tips = centre + directions * reach[:, None]
        return centre, moments, tails, tips

    def arrows(doc: DocumentView, drawing: DrawContext) -> SceneItems | None:
        found = axes(doc)
        if found is None:
            return None
        _, _, tails, tips = found
        rows = [arrow_rows([tail], [tip], RADIUS, tuple(int(color[k:k + 2], 16) / 255 for k in (1, 3, 5)))
                for tail, tip, color in zip(tails, tips, COLORS, strict=True)]
        return SceneItems(arrows=np.concatenate(rows))

    def labels(doc: DocumentView, drawing: DrawContext) -> Annotations | None:
        found = axes(doc)
        if found is None:
            return None
        centre, moments, _, tips = found
        return Annotations(labels=[Label(f"{moment:.1f}", tip, color, away_from=centre)
                                   for moment, tip, color in zip(moments, tips, COLORS, strict=True)])

    def status(doc: DocumentView) -> str | None:
        found = axes(doc)
        return rotor(found[1]) if found is not None else None

    layers = [context.add_scene_items(arrows), context.add_annotations(labels), context.add_status(status)]

    def switched(doc: DocumentView, on: bool) -> None:
        for layer in layers:
            layer.update(doc)

    toggle = context.add_toggle("View > Principal Axes", "principal_axes", switched, after="Hydrogen Bonds")
