"""Meshes, arrows and surfaces that plugins draw with the molecule (provisional).

These build what :class:`~miles.plugin.SceneItems` holds. They are public,
but provisional: they follow the renderer of MILES and may change between
releases. The stable interface is :mod:`miles.plugin`.

``Mesh`` (vertices, normals, colours and triangles); ``arrow_rows``
(arrows for ``SceneItems.arrows``); ``charge_colors`` (atom colours from
partial charges); ``surface_mesh`` with ``SurfaceRequest`` and
``SurfaceColoring`` (a molecular surface) and ``ionisable_charges`` (the
charges of ionisable groups that colour a protein's potential);
``volume_mesh`` with
``VolumeRequest``, ``MAPPED_COLORS`` and ``symmetric_range`` (an
isosurface of a grid, and the colour range of a grid mapped onto it);
``LabelMode`` and ``RIBBON_ATOMS`` (the backbone atoms a cartoon replaces).
"""

import numpy as np
from numpy.typing import ArrayLike

from miles.core.surface import SurfaceColoring, ionisable_charges
from miles.core.volume import symmetric_range
from miles.plugins.colors import ColorLike, color_rgb
from miles.render.arrows import arrow_rows as _arrow_rows
from miles.render.cartoon import Mesh
from miles.render.scene import LabelMode, charge_colors
from miles.render.surface import MAPPED_COLORS, SurfaceRequest, VolumeRequest, surface_mesh, volume_mesh
from miles.ui.canvas.drawing import RIBBON_ATOMS


def arrow_rows(tails: ArrayLike, tips: ArrayLike, radius: float, color: ColorLike) -> np.ndarray:
    """Rows for ``SceneItems.arrows``: arrows of one radius and colour.

    Parameters
    ----------
    tails, tips : array_like
        Start and end of each arrow in angstrom, shape ``(n, 3)``.
    radius : float
        Shaft radius in angstrom.
    color : str, QColor or tuple
        Colour in any form :func:`miles.plugin.color` reads.

    Returns
    -------
    numpy.ndarray
        One row per arrow.
    """
    return _arrow_rows(tails, tips, radius, color_rgb(color))


__all__ = [
    "MAPPED_COLORS",
    "RIBBON_ATOMS",
    "LabelMode",
    "Mesh",
    "SurfaceColoring",
    "SurfaceRequest",
    "VolumeRequest",
    "arrow_rows",
    "charge_colors",
    "ionisable_charges",
    "surface_mesh",
    "symmetric_range",
    "volume_mesh",
]
