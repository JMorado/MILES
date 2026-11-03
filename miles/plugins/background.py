"""Slow drawings: results worked out on a worker thread for each tab (:class:`BackgroundResults`)."""

from collections.abc import Callable, Hashable
from typing import TYPE_CHECKING, Any

from miles.api import DocumentView, Severity

if TYPE_CHECKING:
    from miles.plugins.context import PluginContext


class BackgroundResults:
    """Results of slow work for each tab, worked out off the interface thread and shown when they arrive.

    A layer whose drawing takes long, such as a molecular surface, asks
    for its result with :meth:`get` instead of computing it: on screen the
    work goes to a worker thread (with the busy bar and its Cancel link),
    the layer draws the previous result or nothing meanwhile, and it is
    drawn again when the result arrives. Exported pictures and movie frames
    pass ``now=True`` and get the exact result at once; while atoms are
    dragged (``drawing.dragging``) a layer shows :meth:`latest` and asks
    for nothing. Only the latest
    request of a tab is worked out; older ones are dropped. A panel uses it
    the same way, with its own function to call when a result arrives.
    Create one with :meth:`~miles.plugins.PluginContext.background_results`.
    """

    def __init__(self, context: "PluginContext", ready: Callable[[DocumentView], Any], name: str,
                 message: str) -> None:
        self._context = context
        self._ready = ready
        self._name = f"background:{name}"
        self._message = message

    def _store(self, document: DocumentView) -> dict[str, Any]:
        return self._context.document_cache(document).setdefault(self._name, {})

    def get(self, document: DocumentView, key: Hashable, job: Callable[..., Any], *args: Any,
            now: bool = False) -> Any:
        """Return the result of ``job(*args)`` for a tab, or the previous result while it is worked out.

        Parameters
        ----------
        document : DocumentView
            The tab.
        key : hashable
            Says what the result depends on, e.g. ``(doc.layout_revision,
            doc.revision, setting)``; a new key starts new work.
        job : callable
            Called as ``job(*args)`` on a worker thread: give it copies of
            what it needs (``doc.molecule.copy()``), never the document.
        *args
            Its arguments.
        now : bool, optional
            Work it out at once on this thread (for exported pictures and
            movie frames, ``drawing.export``).

        Returns
        -------
        object
            The result for ``key``; while it is worked out, the latest result
            of an earlier key, or None.
        """
        store = self._store(document)
        if store.get("key") == key:
            return store.get("value")
        if now:
            value = job(*args)
            store.update(key=key, value=value)
            return value
        task = store.get("task")
        if task is not None and task.cancelled:   # the user cancelled: do not start it again by itself
            store.pop("task")
            store["cancelled"] = store.pop("running", None)
        if store.get("cancelled") == key:
            return store.get("value")
        store["wanted"] = (key, job, args)
        if "running" not in store:
            self._start(document, store)
        return store.get("value")

    def has(self, document: DocumentView, key: Hashable) -> bool:
        """Whether the result for ``key`` has arrived for a tab."""
        return bool(self._store(document).get("key") == key)

    def latest(self, document: DocumentView) -> Any:
        """Return the latest result that arrived for a tab (None if none has)."""
        return self._store(document).get("value")

    def forget(self, document: DocumentView) -> None:
        """Drop what is kept for a tab, so that the next :meth:`get` works it out again."""
        self._context.document_cache(document).pop(self._name, None)

    def _start(self, document: DocumentView, store: dict[str, Any]) -> None:
        key, job, args = store.pop("wanted")
        store["running"] = key
        app = self._context.app

        def finished(value: Any, error: str = "") -> None:
            store.pop("running", None)
            store.pop("task", None)
            if document not in app.documents:
                return
            store.update(key=key, value=value)
            if error:
                app.message(f"{self._message} failed: {error}", severity=Severity.WARNING)
            if "wanted" in store:
                if store["wanted"][0] == key:
                    store.pop("wanted")
                else:
                    self._start(document, store)
            self._ready(document)

        store["task"] = app.run_in_background(job, *args, message=self._message, done=finished,
                                              failed=lambda text: finished(None, text))
