"""The molecule canvas.

:mod:`~miles.ui.canvas.canvas` holds the canvas widget, which owns its
parts: :mod:`~miles.ui.canvas.view` (the camera),
:mod:`~miles.ui.canvas.scene_cache` (the cached drawing arrays),
:mod:`~miles.ui.canvas.drawing` (the frame and what plugins draw),
:mod:`~miles.ui.canvas.input` (mouse and keyboard) and
:mod:`~miles.ui.canvas.commands` (what the menus do).
:mod:`~miles.ui.canvas.surfaces` draws the molecule, ``tools`` holds what
each tool does with clicks and drags, ``overlay`` the text and marks over
the molecule, ``context_menu`` the right-click menu and ``view_state``
what sessions keep of a canvas.
"""

from miles.ui.canvas.canvas import Canvas
from miles.ui.canvas.surfaces import (
    SOFTWARE_ENV,
    GLView,
    PainterSurface,
    Surface,
    configure_gl,
    gl_available,
    make_surface,
)

__all__ = ["SOFTWARE_ENV", "Canvas", "GLView", "PainterSurface", "Surface", "configure_gl", "gl_available",
           "make_surface"]
