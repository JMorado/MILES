"""What plugins add to the picture: annotations drawn over it and items drawn in 3D.

A plugin describes these with plain data, anchored to atoms or to points in
angstrom; the canvas projects and draws them on screen, in exported
pictures and in movie frames.

* :class:`Annotations` are lines, arcs and text labels drawn over the
  molecule, like the pinned measurements, with optional tooltips.
* :class:`SceneItems` are drawn with the molecule: arrows
  (:mod:`miles.render.arrows` rows), extra meshes such as isosurfaces,
  and colours that replace the element colours of the atoms.
"""

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from miles.render.arrows import COLUMNS
from miles.render.cartoon import Mesh
from miles.render.theme import Theme

Anchor = int | Sequence[float] | np.ndarray
"""Where an annotation is: an atom index, or a point in angstrom."""


@dataclass(frozen=True)
class Line:
    """A line between two anchors, drawn over the molecule.

    Attributes
    ----------
    start, end : Anchor
        Atom indices or points in angstrom.
    color : str
        Colour such as ``"#3A9BDC"`` (``"#AARRGGBB"`` with transparency);
        empty for the theme's text colour.
    width : float
        Pen width in pixels.
    dashes : tuple of float
        Dash and gap lengths in pen widths, e.g. ``(4, 2)``; empty for a
        solid line.
    round_caps : bool
        Round the ends of the line (and of each dash).
    drawn_only : bool
        Leave the line out when an atom it is anchored to is not shown.
    tip : str
        Shown when the cursor rests on the line.
    ribbon : bool
        An end at a backbone atom drawn only as part of a cartoon moves to
        the nearest point of the ribbon; lines to other atoms that are not
        drawn are left out, as are lines that end up shorter than 0.3
        angstrom. This is how hydrogen bonds to a cartoon are drawn.
    """

    start: Anchor
    end: Anchor
    color: str = ""
    width: float = 1.5
    dashes: tuple[float, ...] = ()
    round_caps: bool = True
    drawn_only: bool = False
    ribbon: bool = False
    tip: str = ""


@dataclass(frozen=True)
class Label:
    """Text centred on an anchor, with a halo in the background colour so it stays readable.

    Attributes
    ----------
    text : str
        What to write.
    at : Anchor
        Atom index or point in angstrom.
    color : str
        Colour, or empty for the theme's text colour.
    size : float
        Font size in points.
    bold, italic : bool
        Font style.
    offset : tuple of float
        Shift in pixels, ``(right, down)``.
    away_from : Anchor or None
        When given, the label moves along the line on screen from this
        anchor to ``at`` until its near edge is ``gap`` pixels past ``at``,
        e.g. just beyond the tip of an arrow whatever the view, without
        covering the tip.
    gap : float
        Distance in pixels for ``away_from``, from ``at`` to the edge of the text.
    halo : float
        Opacity of the halo, 0 for none.
    drawn_only : bool
        Leave the label out when its atom is not shown.
    tip : str
        Shown when the cursor rests on the label.
    """

    text: str
    at: Anchor
    color: str = ""
    size: float = 10.5
    bold: bool = True
    italic: bool = False
    offset: tuple[float, float] = (0.0, 0.0)
    away_from: Anchor | None = None
    gap: float = 8.0
    halo: float = 0.85
    drawn_only: bool = False
    tip: str = ""


@dataclass(frozen=True)
class Arc:
    """An arc around an anchor, from the direction of one anchor to the direction of another.

    Drawn on screen with a fixed radius, the shorter way round, it marks an
    angle such as a bond angle.

    Attributes
    ----------
    center, start, end : Anchor
        The vertex and the two anchors whose directions the arc joins.
    color : str
        Colour, or empty for the theme's text colour.
    radius : float
        Radius in pixels.
    width : float
        Pen width in pixels.
    dashes : tuple of float
        Dash pattern in pen widths; empty for a solid arc.
    drawn_only : bool
        Leave the arc out when an atom it is anchored to is not shown.
    tip : str
        Shown when the cursor rests on the arc.
    """

    center: Anchor
    start: Anchor
    end: Anchor
    color: str = ""
    radius: float = 18.0
    width: float = 1.5
    dashes: tuple[float, ...] = ()
    drawn_only: bool = False
    tip: str = ""


@dataclass
class Annotations:
    """Lines, arcs and labels a plugin draws over the molecule of one tab."""

    lines: list[Line] = field(default_factory=list)
    labels: list[Label] = field(default_factory=list)
    arcs: list[Arc] = field(default_factory=list)


@dataclass(frozen=True)
class MeshItem:
    """A triangle mesh drawn with the molecule, like the molecular surface.

    Attributes
    ----------
    mesh : Mesh
        Vertices in angstrom, normals, colours and triangles
        (:class:`~miles.render.cartoon.Mesh`).
    alpha : float
        Opacity; below 1 the mesh is see-through and drawn after
        everything opaque.
    """

    mesh: Mesh
    alpha: float = 1.0


@dataclass
class SceneItems:
    """What a plugin draws in 3D for one tab.

    Attributes
    ----------
    arrows : numpy.ndarray
        Arrow rows from :func:`miles.render.arrows.arrow_rows`.
    meshes : list of MeshItem
        Extra meshes. Return the same :class:`MeshItem` while nothing
        changed, so it is not uploaded to the graphics card again.
    atom_colors : numpy.ndarray or None
        RGB colour of every atom, shape ``(n_atoms, 3)``, in place of the
        element colours.
    """

    arrows: np.ndarray = field(default_factory=lambda: np.zeros((0, COLUMNS), np.float32))
    meshes: list[MeshItem] = field(default_factory=list)
    atom_colors: np.ndarray | None = None


@dataclass(frozen=True)
class DrawContext:
    """What a plugin's drawing function is drawing for.

    Attributes
    ----------
    theme : Theme
        Colours: the window's theme on screen, the export theme for
        pictures and movies (e.g. on a white background).
    export : bool
        True for exported pictures and movie frames, where every frame
        must be exact; on screen a function may reuse slow results for a
        moment while frames play.
    dragging : bool
        Atoms are being dragged: a function may keep showing a slow
        result (a surface, say) and work it out again when the drag ends,
        when it is called once more.
    """

    theme: Theme
    export: bool = False
    dragging: bool = False
