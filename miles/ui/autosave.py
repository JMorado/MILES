"""Autosave: copies of modified documents written every few minutes, offered again after a crash.

Each running MILES window writes its copies into a folder of its own
under :func:`autosave_folder`, which holds a lock file for as long as the
window is open. A folder whose lock can be taken belongs to a MILES
that is no longer running, so its documents are offered for recovery at
the next start; folders of windows still open are left alone. The copies
are session files (see :mod:`miles.files.session`), so frames, the
view and hidden atoms come back as well. A copy is deleted when its
document is saved or closed, and the whole folder when the window closes
normally.
"""

import gzip
import json
import os
import shutil
import sys
import time
import uuid
from collections.abc import Callable, Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QLockFile, QObject, QSettings, QStandardPaths, QTimer, Signal

from miles.core.molecule import Molecule
from miles.files import session
from miles.files.common import JobInfo
from miles.ui.canvas.view_state import view_state
from miles.ui.document import Change, Document

if TYPE_CHECKING:
    from miles.ui.canvas import Canvas

FOLDER_ENV = "MILES_AUTOSAVE_DIR"
"""Environment variable that moves the autosave folder (used by the tests)."""
INTERVAL_KEY = "autosave/interval"
"""Settings key: minutes between autosaves, 0 for never."""
DEFAULT_INTERVAL = 1
INTERVALS = (1, 2, 5, 10, 0)
"""Choices offered in Settings, in minutes (0 switches autosave off)."""
LOCK_NAME = "lock"
EMPTY_AGE = 600   # seconds before an empty folder of a stopped window is deleted
COMPRESSION = 1   # fast gzip level: the copies are short-lived, speed matters more than size


def autosave_folder() -> Path:
    """Folder that holds the autosave folders of all running windows (it may not exist yet).

    Returns
    -------
    pathlib.Path
        ``$MILES_AUTOSAVE_DIR`` if set, otherwise ``autosave`` in the
        application data folder, next to the plugin folder.
    """
    override = os.environ.get(FOLDER_ENV)
    if override:
        return Path(override)
    base = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
    return Path(base or Path.home() / ".miles") / "autosave"


def interval_text(minutes: int) -> str:
    """How an autosave interval reads in Settings, e.g. "Every 5 minutes" or "Off"."""
    if minutes <= 0:
        return "Off"
    return "Every minute" if minutes == 1 else f"Every {minutes} minutes"


@dataclass
class Leftover:
    """A document autosaved by a MILES that stopped before it was saved or closed.

    Attributes
    ----------
    name : str
        Name of its tab.
    path : str or None
        File it was opened from or last saved to, None for a new document.
    time : float
        When it was autosaved, in seconds since the epoch.
    file : pathlib.Path
        The autosaved session file.
    """

    name: str
    path: str | None
    time: float
    file: Path


@dataclass
class _Entry:
    """What the autosave knows about one open document."""

    serial: int
    dirty: bool = True     # changed since the last copy was written
    written: bool = False  # a copy exists on disk


@dataclass
class _Claim:
    """Autosave folder of a MILES that is gone, locked while its documents are offered."""

    folder: Path
    lock: QLockFile
    leftovers: list[Leftover] = field(default_factory=list)


def snapshot(canvas: "Canvas") -> dict[str, Any]:
    """Session data of a tab's document with its frames and view, safe to write on another thread."""
    doc = canvas.doc
    return session.to_dict(doc.mol, doc.session_info(view_state(canvas)))


def _write_atomically(path: Path, payload: bytes) -> None:
    """Replace ``path`` with ``payload`` so that a reader never sees half a file."""
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "wb") as fh:
        fh.write(payload)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temporary, path)


def _write_copy(folder: Path, serial: int, data: dict[str, Any], about: dict[str, Any]) -> None:
    """Write a document's session file and its description (runs on the writer thread)."""
    folder.mkdir(parents=True, exist_ok=True)
    payload = gzip.compress(json.dumps(data, separators=(",", ":")).encode("utf-8"), COMPRESSION)
    _write_atomically(folder / f"{serial}.{session.EXTENSION}", payload)
    _write_atomically(folder / f"{serial}.json", json.dumps(about).encode("utf-8"))


def _remove_copy(folder: Path, serial: int) -> None:
    for suffix in (session.EXTENSION, "json"):
        (folder / f"{serial}.{suffix}").unlink(missing_ok=True)


def _leftovers(folder: Path) -> list[Leftover]:
    """Documents autosaved in a folder, newest first; unreadable descriptions are skipped."""
    found = []
    for about_file in folder.glob("*.json"):
        copy = about_file.with_suffix(f".{session.EXTENSION}")
        if not copy.exists():
            continue
        try:
            about = json.loads(about_file.read_text(encoding="utf-8"))
            found.append(Leftover(str(about.get("name") or "Untitled"), about.get("path") or None,
                                  float(about.get("time") or copy.stat().st_mtime), copy))
        except (OSError, ValueError, AttributeError):
            continue
    return sorted(found, key=lambda leftover: leftover.time, reverse=True)


def _age(folder: Path) -> float:
    """Seconds since anything in a folder changed."""
    try:
        times = [folder.stat().st_mtime, *(p.stat().st_mtime for p in folder.iterdir())]
    except OSError:
        return 0.0
    return time.time() - max(times)


def read_leftover(leftover: Leftover) -> tuple[Molecule, JobInfo]:
    """Read an autosaved document back.

    Raises
    ------
    miles.files.common.FileFormatError
        If the copy cannot be read.
    """
    return session.read(str(leftover.file))


class Autosave(QObject):
    """Writes copies of a window's modified documents and finds those left by a crash.

    Copies are written every :attr:`interval` minutes, only for documents
    changed since their last copy. The session data is put together on the
    interface thread, which takes milliseconds; compressing and writing it,
    the slow part for large structures, happens on a writer thread.

    Parameters
    ----------
    canvases : callable
        Returns the window's open tabs.
    settings : QSettings
        Where the interval is kept.
    """

    failed = Signal(str)   # a copy could not be written; the message says why

    def __init__(self, canvases: Callable[[], list["Canvas"]], settings: QSettings) -> None:
        super().__init__()
        self._canvases = canvases
        self.settings = settings
        self.folder = self._new_folder()
        self._lock: QLockFile | None = None   # taken when the folder is first written to
        self._entries: dict[Document, _Entry] = {}
        self._serial = 0
        self._writer = ThreadPoolExecutor(max_workers=1, thread_name_prefix="autosave")
        self._pending: list[Future[None]] = []
        self._claims: list[_Claim] = []
        self._reported: set[str] = set()
        self._closed = False
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.save_now)
        self._apply_interval()

    # -- interval ---------------------------------------------------------------------------

    @property
    def interval(self) -> int:
        """Minutes between autosaves; 0 means autosave is off."""
        try:
            return max(0, int(str(self.settings.value(INTERVAL_KEY, DEFAULT_INTERVAL))))
        except ValueError:
            return DEFAULT_INTERVAL

    def set_interval(self, minutes: int) -> None:
        """Autosave every ``minutes`` minutes from now on, or never (0), and remember it.

        Switching autosave off also deletes the copies already written.
        """
        self.settings.setValue(INTERVAL_KEY, int(minutes))
        self._apply_interval()
        if minutes <= 0:
            for doc in list(self._entries):
                self._drop_copy(doc)

    def _apply_interval(self) -> None:
        minutes = self.interval
        if minutes > 0 and not self._closed:
            self.timer.start(minutes * 60_000)
        else:
            self.timer.stop()

    # -- documents --------------------------------------------------------------------------

    def track(self, doc: Document) -> None:
        """Start watching a new tab's document."""
        self._serial += 1
        self._entries[doc] = _Entry(self._serial)
        doc.changed.connect(lambda kind, d=doc: self._changed(d, kind))
        doc.renamed.connect(lambda d=doc: self._renamed(d))

    def mark_changed(self, doc: Document) -> None:
        """Write a copy of ``doc`` at the next autosave even though no edit was announced."""
        entry = self._entries.get(doc)
        if entry is not None:
            entry.dirty = True

    def _changed(self, doc: Document, kind: str) -> None:
        if kind != Change.FRAME:   # showing another frame is not an edit
            self.mark_changed(doc)

    def _renamed(self, doc: Document) -> None:
        if not doc.modified:   # just saved: the file holds the changes now
            self._drop_copy(doc)

    def forget(self, doc: Document) -> None:
        """Stop watching a document whose tab was closed and delete its copy."""
        self._drop_copy(doc)
        self._entries.pop(doc, None)

    def _drop_copy(self, doc: Document) -> None:
        entry = self._entries.get(doc)
        if entry is None:
            return
        entry.dirty = True
        if entry.written:
            entry.written = False
            self._submit(_remove_copy, self.folder, entry.serial)

    def save_now(self) -> int:
        """Write a copy of every modified document changed since its last copy.

        Returns
        -------
        int
            Number of copies started (they finish on the writer thread).
        """
        if self._closed or self.interval <= 0:
            return 0
        started = 0
        for canvas in self._canvases():
            doc = canvas.doc
            entry = self._entries.get(doc)
            if entry is None:
                continue
            if not doc.modified or not len(doc.mol):
                if entry.written:
                    self._drop_copy(doc)
                continue
            if entry.written and not entry.dirty:
                continue
            if not self._take_lock():
                return started
            about = {"name": doc.name, "path": doc.path, "time": time.time()}
            self._submit(_write_copy, self.folder, entry.serial, snapshot(canvas), about)
            entry.dirty = False
            entry.written = True
            started += 1
        return started

    def wait(self) -> None:
        """Block until every copy started so far is on disk (for tests and shutdown)."""
        for future in list(self._pending):
            future.result()
        self._pending = [f for f in self._pending if not f.done()]

    @staticmethod
    def _new_folder() -> Path:
        return autosave_folder() / f"{os.getpid()}-{uuid.uuid4().hex[:8]}"

    def _take_lock(self) -> bool:
        if self._lock is not None:
            return True
        for _ in range(3):   # another window may lock a folder in the instant between mkdir and lock
            try:
                self.folder.mkdir(parents=True, exist_ok=True)
            except OSError as exc:
                self._report(f"Could not autosave: {exc}")
                return False
            lock = QLockFile(str(self.folder / LOCK_NAME))
            lock.setStaleLockTime(0)
            if lock.tryLock(0):
                self._lock = lock
                return True
            self.folder = self._new_folder()
        self._report("Could not autosave: the autosave folder is in use")
        return False

    def _submit(self, func: Callable[..., None], *args: Any) -> None:
        if self._closed:
            return
        future = self._writer.submit(func, *args)
        future.add_done_callback(self._done)
        self._pending = [f for f in self._pending if not f.done()]
        self._pending.append(future)

    def _done(self, future: "Future[None]") -> None:
        exc = future.exception()
        if exc is not None:
            self.failed.emit(f"Could not autosave: {exc}")   # queued to the interface thread

    def _report(self, message: str) -> None:
        if message not in self._reported:
            self._reported.add(message)
            print(message, file=sys.stderr)
            self.failed.emit(message)

    def close(self) -> None:
        """Wait for the writer and delete this window's copies, as the window closed normally."""
        if self._closed:
            return
        self._closed = True
        self.timer.stop()
        self._writer.shutdown(wait=True)
        self._pending.clear()
        if self._lock is not None:
            self._lock.unlock()
            self._lock = None
        shutil.rmtree(self.folder, ignore_errors=True)
        self.release_leftovers()

    # -- recovery ---------------------------------------------------------------------------

    def find_leftovers(self) -> list[Leftover]:
        """Documents left by MILES windows that are no longer running, newest first.

        The folders they are in stay locked until :meth:`discard_leftovers`
        or :meth:`release_leftovers`, so another window starting at the same
        time does not offer them too. Empty folders of stopped windows are
        deleted on the way.
        """
        root = autosave_folder()
        if not root.is_dir():
            return []
        claimed = {claim.folder for claim in self._claims}
        own = f"{os.getpid()}-"   # this process is running, so its windows did not crash
        for folder in sorted(root.iterdir()):
            if not folder.is_dir() or folder.name.startswith(own) or folder in claimed:
                continue
            age = _age(folder)   # before locking it, which touches it
            lock = QLockFile(str(folder / LOCK_NAME))
            lock.setStaleLockTime(0)   # only a process that is gone gives its lock up
            if not lock.tryLock(0):
                continue   # that MILES is still running
            leftovers = _leftovers(folder)
            if not leftovers:
                lock.unlock()
                if age > EMPTY_AGE:   # a younger one may belong to a window just starting
                    shutil.rmtree(folder, ignore_errors=True)
                continue
            self._claims.append(_Claim(folder, lock, leftovers))
        found = [leftover for claim in self._claims for leftover in claim.leftovers]
        return sorted(found, key=lambda leftover: leftover.time, reverse=True)

    def discard_leftovers(self, keep: Sequence[Leftover] = ()) -> None:
        """Delete the documents found by :meth:`find_leftovers` (after restoring or discarding them).

        Parameters
        ----------
        keep : sequence of Leftover, optional
            Documents to offer again at the next start, such as copies that
            could not be read.
        """
        kept = {leftover.file for leftover in keep}
        for claim in self._claims:
            claim.lock.unlock()
            if not any(leftover.file in kept for leftover in claim.leftovers):
                shutil.rmtree(claim.folder, ignore_errors=True)
                continue
            for leftover in claim.leftovers:
                if leftover.file not in kept:
                    leftover.file.unlink(missing_ok=True)
                    leftover.file.with_suffix(".json").unlink(missing_ok=True)
        self._claims.clear()

    def release_leftovers(self) -> None:
        """Keep the documents found for the next start, unlocking their folders."""
        for claim in self._claims:
            claim.lock.unlock()
        self._claims.clear()
