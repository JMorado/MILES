"""Light and dark colour themes, including the element palette."""

from dataclasses import dataclass, field, replace
from enum import StrEnum

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QGuiApplication

from miles.core.elements import SYMBOLS, element

# Softer versions of the usual element colours; elements not listed get a
# desaturated CPK colour.
_ELEMENT_COLORS = {
    "H": "#F4F5F7", "C": "#6E747E", "N": "#3D6FE3", "O": "#E5484D", "F": "#5DC28C",
    "Cl": "#33A86A", "Br": "#B4553D", "I": "#7A52D6", "S": "#EDB43D", "P": "#F08A2A",
    "B": "#F2A07E", "Si": "#C9A57D", "Se": "#E79A3A", "Li": "#A57BE0", "Na": "#9C6BDB",
    "K": "#8A5BD1", "Mg": "#6CC07A", "Ca": "#5DB06C", "Al": "#B7A3A6", "Fe": "#D9733F",
    "Co": "#E58AA0", "Ni": "#5FBF8F", "Cu": "#C98A4B", "Zn": "#8C93B8", "Ru": "#3E8C8F",
    "Rh": "#2F7F8A", "Pd": "#2F6E85", "Ag": "#B9BEC6", "Pt": "#A9AEC0", "Au": "#E8B931",
    "Hg": "#A5A9C2", "Ti": "#AEB4BC", "Mn": "#9C7AC7", "Cr": "#8A99C7", "Sn": "#7D8F8F",
    "Xe": "#4FA3B6", "He": "#BFE8EE", "Ne": "#A9D9EC", "Ar": "#7FC6DA",
}


class Appearance(StrEnum):
    """Which theme the user asked for."""

    SYSTEM = "system"
    LIGHT = "light"
    DARK = "dark"


class Background(StrEnum):
    """Background of exported pictures."""

    THEME = "theme"
    WHITE = "white"
    TRANSPARENT = "transparent"


def _desaturate(hex_color: str, amount: float = 0.22) -> QColor:
    c = QColor(hex_color)
    c.setHsvF(max(c.hsvHueF(), 0.0), c.hsvSaturationF() * (1 - amount), c.valueF() * 0.94, c.alphaF())
    return c


def rgb(color: QColor) -> tuple[float, float, float]:
    """Red, green and blue of a colour as floats in [0, 1]."""
    return color.redF(), color.greenF(), color.blueF()


@dataclass
class Theme:
    """Colours of the interface and the molecule.

    Attributes
    ----------
    name : str
        Identifier, part of cache keys.
    background, outline, text, muted : QColor
        Canvas background, atom outlines, text and secondary text.
    panel, panel_border : QColor
        Floating panels.
    accent, hover, attach, danger : QColor
        Selection, hover, attachment point and delete highlights.
    carbon, hydrogen : QColor
        Element colours that differ between light and dark themes.
    """

    name: str
    background: QColor
    outline: QColor
    text: QColor
    muted: QColor
    panel: QColor
    panel_border: QColor
    accent: QColor
    hover: QColor
    attach: QColor
    danger: QColor
    carbon: QColor
    hydrogen: QColor
    _atom_colors: dict[str, QColor] = field(default_factory=dict, repr=False)
    _color_table: np.ndarray | None = field(default=None, init=False, repr=False, compare=False)

    def element_colors(self) -> np.ndarray:
        """Colour of every element as floats in [0, 1], indexed by atomic number, shape ``(n, 3)``."""
        if self._color_table is None:
            self._color_table = np.array([self.rgb(symbol) for symbol in SYMBOLS], dtype=np.float32)
        return self._color_table

    @property
    def dark(self) -> bool:
        """Whether the background is dark."""
        return self.background.lightness() < 128

    def atom_color(self, symbol: str) -> QColor:
        """Colour of an element."""
        color = self._atom_colors.get(symbol)
        if color is None:
            if symbol == "C":
                color = QColor(self.carbon)
            elif symbol == "H":
                color = QColor(self.hydrogen)
            elif symbol in _ELEMENT_COLORS:
                color = QColor(_ELEMENT_COLORS[symbol])
            else:
                color = _desaturate(element(symbol).color)
            self._atom_colors[symbol] = color
        return color

    def rgb(self, symbol: str) -> tuple[float, float, float]:
        """Colour of an element as floats in [0, 1]."""
        return rgb(self.atom_color(symbol))

    def highlight_colors(self) -> list[QColor]:
        """Colour of each highlight, indexed by ``scene.Highlight``."""
        return [self.outline, self.accent, self.hover, self.attach, self.danger]


LIGHT = Theme(
    name="light",
    background=QColor("#F6F7F9"),
    outline=QColor("#1D2026"),
    text=QColor("#1D2026"),
    muted=QColor("#8B909A"),
    panel=QColor(255, 255, 255),
    panel_border=QColor(0, 0, 0, 28),
    accent=QColor("#3A7BFF"),
    hover=QColor("#7FA9FF"),
    attach=QColor("#22B573"),
    danger=QColor("#EF4C4C"),
    carbon=QColor("#6E747E"),
    hydrogen=QColor("#F4F5F7"),
)

DARK = Theme(
    name="dark",
    background=QColor("#141619"),
    outline=QColor("#050607"),
    text=QColor("#E7E9EC"),
    muted=QColor("#7D838D"),
    panel=QColor(32, 35, 40),
    panel_border=QColor(255, 255, 255, 24),
    accent=QColor("#5B93FF"),
    hover=QColor("#8DB4FF"),
    attach=QColor("#35C985"),
    danger=QColor("#FF6363"),
    carbon=QColor("#9AA0AA"),
    hydrogen=QColor("#E9ECF0"),
)

_WHITE = replace(LIGHT, name="white", background=QColor(255, 255, 255), _atom_colors={})


def system_prefers_dark() -> bool:
    """Whether the operating system uses a dark appearance."""
    if QGuiApplication.instance() is None:
        return False
    return QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark


def theme_for(appearance: str) -> Theme:
    """Theme for an :class:`Appearance` choice."""
    if appearance == Appearance.DARK or (appearance == Appearance.SYSTEM and system_prefers_dark()):
        return DARK
    return LIGHT


def export_theme(theme: Theme, background: str) -> Theme:
    """Theme for an exported picture with a :class:`Background` choice."""
    return _WHITE if background == Background.WHITE else theme
