"""Colours as plugins give them: a name, a QColor or a tuple, turned into what drawing needs.

Every place of the plugin interface that takes a colour accepts any of
these forms (see :func:`color`), so a plugin can pass a theme colour, a
hex string or numbers straight through.
"""

from collections.abc import Sequence
from dataclasses import fields, is_dataclass, replace
from typing import Any

from PySide6.QtGui import QColor

ColorLike = str | QColor | Sequence[float] | Sequence[int]
"""A colour: a name such as ``"#3a7bd5"`` or ``"tomato"``, a :class:`QColor`,
three or four floats from 0 to 1, or three or four integers from 0 to 255."""


def color(value: ColorLike) -> QColor:
    """Return a colour given in any of the accepted forms (see :data:`ColorLike`) as a QColor.

    Raises
    ------
    ValueError
        If the value is not a colour.
    """
    if isinstance(value, QColor):
        return QColor(value)
    if isinstance(value, str):
        found = QColor(value)
        if not found.isValid():
            raise ValueError(f"{value!r} is not a colour")
        return found
    parts = list(value)
    if len(parts) not in (3, 4):
        raise ValueError(f"A colour needs three or four components, not {len(parts)}")
    if all(isinstance(p, int) and not isinstance(p, bool) for p in parts):
        return QColor(*[int(p) for p in parts])
    return QColor.fromRgbF(*[min(1.0, max(0.0, float(p))) for p in parts])


def color_name(value: ColorLike | None) -> str:
    """Return a colour as ``#rrggbb`` (``#aarrggbb`` when it is see-through); empty for none."""
    if value is None or (isinstance(value, str) and not value):
        return ""
    found = color(value)
    return found.name(QColor.NameFormat.HexArgb if found.alpha() < 255 else QColor.NameFormat.HexRgb)


def color_rgb(value: ColorLike) -> tuple[float, float, float]:
    """Return red, green and blue of a colour as floats from 0 to 1."""
    found = color(value)
    return found.redF(), found.greenF(), found.blueF()


def named_colors(item: Any) -> Any:
    """Return a drawing item with every ``color`` field given as a name (items without one are kept)."""
    if not is_dataclass(item) or isinstance(item, type):
        return item
    names = {f.name for f in fields(item)}
    value = getattr(item, "color", None)
    if "color" not in names or isinstance(value, str):
        return item
    return replace(item, color=color_name(value))
