"""Editing state shared by every window: current tool, fragment and mouse mode."""

import time
from enum import StrEnum

import numpy as np
from PySide6.QtCore import QObject, Signal

from miles.core.elements import ELEMENTS
from miles.core.fragments import (
    Category,
    Fragment,
    default_element_fragment,
    element_fragment,
    element_types,
)

TYPING_PAUSE = 0.9
"""Seconds within which a second letter extends an element symbol (``c`` then ``l`` is Cl)."""


class Tool(StrEnum):
    """What a click on the canvas does."""

    SELECT = "select"
    BUILD = "build"
    DELETE = "delete"


class MouseMode(StrEnum):
    """What a plain left drag on empty space does."""

    ROTATE = "rotate"
    TRANSLATE = "translate"
    SCALE = "scale"


MOUSE_MODE_KEYS = {"r": MouseMode.ROTATE, "t": MouseMode.TRANSLATE, "s": MouseMode.SCALE}
"""Letters that switch the mouse mode, as in VMD."""


class Builder(QObject):
    """Current tool, fragment and mouse mode.

    Attributes
    ----------
    fragment : Fragment
        What the build tool places or attaches.
    tool : Tool
        Active tool.
    preview_rotation : numpy.ndarray or None
        Orientation of the fragment preview, also used when placing it.
    last_element : str
        Most recent element fragment, for the bonding-type shortcuts.
    mouse_mode : MouseMode
        What dragging on empty space does.
    """

    fragmentChanged = Signal()
    toolChanged = Signal()
    mouseModeChanged = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.fragment: Fragment = default_element_fragment("C")
        self.tool = Tool.BUILD
        self.preview_rotation: np.ndarray | None = None
        self.last_element = "C"
        self.mouse_mode = MouseMode.ROTATE
        self._typed = ""
        self._typed_at = 0.0
        self._mode_before: MouseMode | None = None   # restored if the letter starts a symbol

    def set_tool(self, tool: Tool) -> None:
        """Switch tools."""
        if tool != self.tool:
            self.tool = tool
            self.toolChanged.emit()

    def set_mouse_mode(self, mode: MouseMode) -> None:
        """Choose what dragging on empty space does."""
        if mode != self.mouse_mode:
            self.mouse_mode = mode
            self.mouseModeChanged.emit()

    def set_fragment(self, fragment: Fragment) -> None:
        """Choose a fragment and switch to the build tool."""
        self.fragment = fragment
        if fragment.category == Category.ELEMENT:
            self.last_element = fragment.symbol
        self.preview_rotation = None
        self.fragmentChanged.emit()
        self.set_tool(Tool.BUILD)

    def set_attach(self, index: int) -> bool:
        """Choose the terminal atom of the fragment that the new bond replaces.

        Parameters
        ----------
        index : int
            A terminal atom, or a heavy atom whose first hydrogen (or terminal
            neighbour) is taken.

        Returns
        -------
        bool
            Whether the attachment point changed.
        """
        mol = self.fragment.mol
        if index >= len(mol):
            return False
        if mol.degree(index) != 1:
            ends = mol.hydrogens_on(index) or mol.terminal_neighbors(index)
            if not ends:
                return False
            index = ends[0]
        self.fragment.attach = index
        self.fragmentChanged.emit()
        return True

    # -- keyboard shortcuts --------------------------------------------------------------------

    def type_letter(self, letter: str) -> Fragment | None:
        """Handle a typed letter: element symbols choose an element, r, t and s the mouse mode.

        A second letter typed quickly completes a two-letter symbol (``c``
        then ``l`` is Cl); if the first letter had switched the mouse mode,
        the previous mode comes back. With the build tool, ``s`` is sulfur.

        Returns
        -------
        Fragment or None
            The element fragment chosen, if any.
        """
        now = time.monotonic()
        if self._typed and now - self._typed_at < TYPING_PAUSE:
            two = (self._typed + letter).capitalize()
            if two in ELEMENTS:
                self._typed = ""
                if self._mode_before is not None:
                    self.set_mouse_mode(self._mode_before)
                return self._use(default_element_fragment(two))
        self._typed, self._typed_at = letter, now
        self._mode_before = None
        one = letter.upper()
        mode = MOUSE_MODE_KEYS.get(letter.lower())
        if one in ELEMENTS and (mode is None or self.tool == Tool.BUILD):
            return self._use(default_element_fragment(one))
        if mode is not None:
            self._mode_before = self.mouse_mode
            self.set_mouse_mode(mode)
        return None

    def end_typing(self) -> None:
        """Forget the last letter, so the next one starts a new symbol (e.g. after a click)."""
        self._typed = ""
        self._mode_before = None

    def choose_bonding_type(self, n: int) -> Fragment | None:
        """Use the ``n``-th bonding type (1-based) of the current element."""
        symbol = self.fragment.symbol if self.fragment.category == Category.ELEMENT else self.last_element
        types = element_types(symbol)
        return self._use(element_fragment(symbol, types[n - 1])) if 1 <= n <= len(types) else None

    def _use(self, fragment: Fragment) -> Fragment:
        self.set_fragment(fragment)
        return fragment
