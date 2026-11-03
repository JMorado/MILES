"""What plugins add to the canvases of a window: drawing layers, status text and context menu items.

The main window owns one :class:`Extensions` and hands it to every canvas
it creates; :class:`~miles.plugins.PluginManager` fills it. The canvas
calls the functions here and caches their results, so it needs to know
nothing about the features that provide them.
"""

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from miles.render.layers import DrawContext

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

SECTIONS = ("atom", "bond", "selection")
"""Places in the canvas's right-click menu where plugins can add items."""


class Layer:
    """A function that the canvases call to draw or describe their document.

    Parameters
    ----------
    function : callable
        ``function(canvas, context)`` returning the layer's data for that
        canvas (:class:`~miles.render.layers.Annotations`,
        :class:`~miles.render.layers.SceneItems` or status text), or
        None. ``context`` is a :class:`~miles.render.layers.DrawContext`
        (None for status text).

    Attributes
    ----------
    generation : int
        Raised by :meth:`invalidate`; part of the key under which canvases
        cache the results.
    """

    def __init__(self, function: Callable[["Canvas", DrawContext | None], Any]) -> None:
        self.function = function
        self.generation = 0

    def invalidate(self) -> None:
        """Make every canvas call the function again before it next draws."""
        self.generation += 1


@dataclass(frozen=True)
class MenuItem:
    """An item a plugin adds to the canvas's right-click menu.

    Attributes
    ----------
    title : str
        Text of the item.
    action : callable or None
        Called without arguments when the item is chosen.
    items : sequence of MenuItem
        When not empty, the item is a submenu holding these.
    enabled : bool
        Whether the item can be chosen.
    """

    title: str
    action: Callable[[], Any] | None = None
    items: Sequence["MenuItem"] = ()
    enabled: bool = True


@dataclass(frozen=True)
class MenuTarget:
    """What was right-clicked.

    Attributes
    ----------
    atom : int or None
        The atom under the cursor, or the first atom of the bond under it.
    bond : tuple of int or None
        The bond under the cursor.
    """

    atom: int | None = None
    bond: tuple[int, int] | None = None


MenuFunction = Callable[["Canvas", MenuTarget], Sequence[MenuItem]]


class Extensions:
    """Layers, status text and menu items added by plugins, in the order they were added."""

    def __init__(self) -> None:
        self.annotations: list[Layer] = []
        self.scene_items: list[Layer] = []
        self.status: list[Layer] = []
        self.menu_items: dict[str, list[MenuFunction]] = {section: [] for section in SECTIONS}
