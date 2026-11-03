"""What sessions keep of a canvas: its camera, display options, charge model and plugin settings."""

import contextlib
import copy
from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import numpy as np

from miles.core.charges import ChargeModel
from miles.core.molecule import Representation
from miles.render.camera import Camera
from miles.render.cartoon import CartoonColor
from miles.render.scene import LabelMode

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

STYLE_CONVERTERS: dict[str, Callable[[Any], Any]] = {
    "mode": Representation, "labels": LabelMode, "cartoon_color": CartoonColor,
    "show_hydrogens": bool, "depth_fade": bool, "stereo": bool,
}
"""Display options sessions keep, with how to read each back."""


def view_state(canvas: "Canvas") -> dict[str, Any]:
    """Camera, display options, charge model and plugin state of a canvas, as JSON-ready data."""
    cam, d = canvas.camera, canvas.display
    return {
        "camera": {"rotation": cam.rotation.tolist(), "center": cam.center.tolist(), "scale": cam.scale,
                   "pan": cam.pan.tolist()},
        "style": {"mode": str(d.mode), "show_hydrogens": d.show_hydrogens, "labels": str(d.labels),
                  "depth_fade": d.depth_fade, "cartoon_color": str(d.cartoon_color), "stereo": d.stereo},
        "charge_model": str(canvas.doc.charge_model),
        "plugins": copy.deepcopy(canvas.plugin_state),
    }


def apply_view_state(canvas: "Canvas", state: dict[str, Any]) -> None:
    """Restore what :func:`view_state` saved; unknown or missing entries are skipped.

    Sessions saved before some features became plugins keep their
    settings in ``style`` (e.g. ``hydrogen_bonds``); the plugins read
    them from ``doc.info.view`` themselves.
    """
    cam = state.get("camera")
    if cam:
        camera = Camera()
        camera.rotation = np.array(cam["rotation"], dtype=float)
        camera.center = np.array(cam["center"], dtype=float)
        camera.scale = float(cam["scale"])
        camera.pan = np.array(cam["pan"], dtype=float)
        canvas.view.stop()
        canvas.view.camera = camera
        canvas.view.needs_fit = False
    style = state.get("style") or {}
    for name, convert in STYLE_CONVERTERS.items():
        if name in style:
            with contextlib.suppress(ValueError):
                setattr(canvas.display, name, convert(style[name]))
    if isinstance(state.get("plugins"), dict):
        canvas.plugin_state.update(copy.deepcopy(state["plugins"]))
    doc = canvas.doc
    with contextlib.suppress(ValueError):
        model = ChargeModel(state.get("charge_model", doc.charge_model))
        if model != ChargeModel.FILE or doc.info.charges is not None:
            doc.charge_model = model   # as the session left it: opening a file is not an edit
    canvas.update_status()
