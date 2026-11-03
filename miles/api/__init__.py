"""Scripting interface to a running MILES window.

This module is what scripts, plugins and the Python console use to drive
the application. It wraps the main window in two small classes:

:class:`App`
    The window: open and create documents, run any menu command by name,
    show a message.
:class:`DocumentView`
    One tab: its :class:`~miles.core.molecule.Molecule`, selection,
    frames, display options and camera, with undoable edits and picture
    export.

The chemistry library :mod:`miles` (building, force fields, geometry,
selections...) works on the molecules these objects hand out; edit them
inside :meth:`DocumentView.edit` so the change can be undone and the view
follows. Calls must come from the interface thread (the console and
plugin commands run there).

Examples
--------
>>> doc = app.current
>>> doc.select(smarts="[OX2H]")          # hydroxyl oxygens
>>> doc.style(mode="ball", labels="symbol")
>>> app.run("Clean Geometry")
>>> doc.export_image("figure.png", scale=3)

The classes live in :mod:`miles.api.app`, :mod:`miles.api.document` and
:mod:`miles.api.picture`; plugins and scripts import them from
:mod:`miles.plugin`.
"""

from miles.api.app import App, Severity, current_app
from miles.api.document import DocumentView, document_view
from miles.api.picture import Picture

__all__ = ["App", "DocumentView", "Picture", "Severity", "current_app", "document_view"]
