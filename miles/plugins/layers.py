"""Drawings and status text of a plugin, worked out again on demand (:class:`PluginLayer`)."""

from typing import TYPE_CHECKING

from miles.api import DocumentView
from miles.ui.canvas.extensions import Layer

if TYPE_CHECKING:
    from miles.plugins.discovery import PluginInfo
    from miles.plugins.manager import PluginManager


class PluginLayer:
    """Annotations, scene items or status text a plugin added, worked out again by :meth:`update`.

    On screen a drawing is kept until the tab's structure, selection, frame,
    display options or theme change, or :meth:`update` is called. When the
    function raises, the error is reported and the function is not called
    again for that tab until the document changes or :meth:`update` is
    called, so a bug does not repeat on every redraw.
    """

    def __init__(self, manager: "PluginManager", kind: str, info: "PluginInfo | None" = None) -> None:
        self._manager = manager
        self.kind = kind
        self.info = info
        self.failed_at: dict[int, tuple[int, int, int]] = {}   # canvas id: (molecule, revision, generation)
        self.layer = Layer(lambda canvas, context: None)

    def update(self, document: DocumentView | None = None) -> None:
        """Call the function again for every tab and redraw (or update only this tab's status line now).

        Parameters
        ----------
        document : DocumentView, optional
            The tab whose drawing or status changed; the others are worked
            out again when they are next drawn.
        """
        self.failed_at.clear()
        self.layer.invalidate()
        canvases = [document.canvas] if document is not None else list(self._manager.window.canvases)
        for canvas in canvases:
            if self.kind == "status":
                canvas.update_status()
            else:
                canvas.redraw()
