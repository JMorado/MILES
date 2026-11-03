"""Slow work on worker threads, for the whole application, with results delivered on the interface thread.

:func:`executor` returns the one :class:`Executor` of the application: a
few worker threads and a queue. A job's result (or error) comes back on the
interface thread, where the ``done`` (or ``failed``) function given with it
is called, unless:

* its :class:`CancelToken` was cancelled (Cancel in the busy bar, or a
  newer request with the same ``key``): the result is dropped, and a job
  that takes the token can stop early by checking it;
* its owner, a Qt object such as the panel that asked for it, was deleted.

No Qt object is ever handed to a worker thread: jobs get plain arguments,
owners are kept only as weak references, and every Qt call happens on the
interface thread. The executor keeps every job until its result has been
delivered, so the last reference to a job (and to what its functions hold)
is always dropped on the interface thread, never on a worker.

:class:`TaskRunner` shows the commands' jobs in the busy bar; panels that
recompute what they show use :meth:`Executor.submit` with a ``key``
directly.
"""

import contextlib
import itertools
import os
import queue
import threading
import weakref
from collections.abc import Callable, Hashable
from dataclasses import dataclass, field
from typing import Any

import shiboken6
from PySide6.QtCore import QCoreApplication, QObject, Signal

WORKERS = max(2, min(4, (os.cpu_count() or 2) - 1))
"""Worker threads of the application."""


class CancelToken:
    """Asks a job to stop. Jobs that take it check :attr:`cancelled` between steps and return early."""

    def __init__(self) -> None:
        self._event = threading.Event()

    def cancel(self) -> None:
        """Ask the job to stop; its result, if it still comes, is dropped."""
        self._event.set()

    @property
    def cancelled(self) -> bool:
        """Whether :meth:`cancel` was called."""
        return self._event.is_set()


@dataclass(eq=False)
class Job:
    """A function waiting for or running on a worker thread, and what to do with its result.

    Attributes
    ----------
    token : CancelToken
        Cancels the job; checked before the result is delivered.
    key : hashable or None
        A newer job with the same key cancels this one.
    """

    func: Callable[..., Any]
    args: tuple[Any, ...]
    token: CancelToken
    done: Callable[[Any], Any] | None = None
    failed: Callable[[str], Any] | None = None
    owner: "weakref.ref[QObject] | None" = None
    key: Hashable | None = None
    pass_token: bool = False
    ended: list[Callable[["Job"], Any]] = field(default_factory=list)   # called once it is over

    def cancel(self) -> None:
        """Drop the result and ask the function to stop (if it takes the token)."""
        self.token.cancel()

    def wanted(self) -> bool:
        """Whether the result is still to be delivered."""
        if self.token.cancelled:
            return False
        if self.owner is not None:
            owner = self.owner()
            return owner is not None and shiboken6.isValid(owner)
        return True


class Executor(QObject):
    """Worker threads shared by the application; see the module documentation.

    It lives on the interface thread for as long as the application; the
    workers only call functions and send their results back.
    """

    _result = Signal(int, object, object)   # job number, result, error message (None when it succeeded)

    def __init__(self, workers: int = WORKERS, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._queue: queue.SimpleQueue[int] = queue.SimpleQueue()
        self._jobs: dict[int, Job] = {}   # until delivered; emptied on the interface thread only
        self._lock = threading.Lock()
        self._numbers = itertools.count()
        self._latest: dict[Hashable, Job] = {}
        self._result.connect(self._deliver)
        for k in range(workers):
            threading.Thread(target=self._work, name=f"miles-worker-{k}", daemon=True).start()

    def submit(self, func: Callable[..., Any], *args: Any, done: Callable[[Any], Any] | None = None,
               failed: Callable[[str], Any] | None = None, owner: QObject | None = None,
               key: Hashable | None = None, pass_token: bool = False) -> Job:
        """Run ``func(*args)`` on a worker thread; ``done(result)`` or ``failed(message)`` follow.

        Parameters
        ----------
        func : callable
            The slow function. It must not touch Qt objects or the
            documents: give it copies of what it needs.
        *args
            Its arguments.
        done, failed : callable, optional
            Called on the interface thread with the result, or with the
            error message when ``func`` raised.
        owner : QObject, optional
            Nothing is delivered once it has been deleted.
        key : hashable, optional
            Cancel the previous job with this key: only the newest request
            is answered (a panel recomputing after each change).
        pass_token : bool, optional
            Call ``func(*args, cancel=token)`` so it can stop early.

        Returns
        -------
        Job
            The job; ``cancel()`` drops its result.
        """
        job = Job(func, args, CancelToken(), done, failed, weakref.ref(owner) if owner is not None else None,
                  key, pass_token)
        if key is not None:
            previous = self._latest.get(key)
            if previous is not None:
                previous.cancel()
            self._latest[key] = job
        number = next(self._numbers)
        with self._lock:
            self._jobs[number] = job
        self._queue.put(number)
        return job

    @property
    def pending(self) -> int:
        """Jobs submitted whose result has not been delivered yet."""
        with self._lock:
            return len(self._jobs)

    def cancel(self, key: Hashable) -> None:
        """Cancel the job last submitted with a key, if it has not ended."""
        job = self._latest.pop(key, None)
        if job is not None:
            job.cancel()

    def _work(self) -> None:
        while True:
            number = self._queue.get()
            with self._lock:
                job = self._jobs[number]
            func, args, token, pass_token = job.func, job.args, job.token, job.pass_token
            del job   # the table keeps it: it is let go on the interface thread
            result: Any = None
            error: str | None = None
            if not token.cancelled:
                try:
                    result = func(*args, **({"cancel": token} if pass_token else {}))
                except Exception as exc:  # reported on the interface thread instead of lost here
                    error = str(exc) or type(exc).__name__
            del func, args, token
            with contextlib.suppress(RuntimeError):   # the application is gone: nobody waits any more
                self._result.emit(number, result, error)
            del result

    def _deliver(self, number: int, result: Any, error: str | None) -> None:
        with self._lock:
            job = self._jobs.pop(number)
        if job.key is not None and self._latest.get(job.key) is job:
            del self._latest[job.key]
        if job.wanted():
            if error is None:
                if job.done is not None:
                    job.done(result)
            elif job.failed is not None:
                job.failed(error)
        for ended in job.ended:
            ended(job)


_executor: Executor | None = None


def executor() -> Executor:
    """Return the application's executor, made on first use (on the interface thread)."""
    global _executor
    if _executor is None or not shiboken6.isValid(_executor):
        _executor = Executor(parent=QCoreApplication.instance())
    return _executor


class Task(QObject):
    """A command's job, shown in the busy bar until it ends.

    Created by :meth:`TaskRunner.start` (scripts and plugins use
    :meth:`miles.api.App.run_in_background`).

    Attributes
    ----------
    message : str
        What the busy bar says while the task runs, e.g. "Searching for
        conformers".
    cancelled : bool
        Whether :meth:`cancel` was called; the result is then dropped.
    running : bool
        Whether the job has not ended yet.
    """

    ended = Signal()

    def __init__(self, message: str) -> None:
        super().__init__()
        self.message = message
        self.cancelled = False
        self.running = True
        self.job: Job | None = None

    def cancel(self) -> None:
        """Drop the result and ask the job to stop; a job that takes the token stops at its next check."""
        if self.job is not None:
            self.job.cancel()
        self.cancelled = True
        self._end()

    def _end(self) -> None:
        if self.running:
            self.running = False
            self.ended.emit()


class TaskRunner(QObject):
    """The commands' tasks of a window that are still running.

    Signals
    -------
    changed()
        A task started, ended or was cancelled.
    """

    changed = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.tasks: list[Task] = []

    def start(self, message: str, job: Callable[..., Any], *args: Any,
              done: Callable[[Any], Any] | None = None, failed: Callable[[str], Any] | None = None,
              pass_token: bool = False) -> Task:
        """Run ``job(*args)`` on a worker thread and call ``done(result)`` or ``failed(message)`` afterwards.

        Parameters
        ----------
        message : str
            What the busy bar says meanwhile.
        job : callable
            The slow function. It must not touch the interface or the
            documents: give it copies of what it needs.
        *args
            Its arguments.
        done : callable, optional
            Called with the result on the interface thread, unless the task
            was cancelled.
        failed : callable, optional
            Called with the error message when the function raised.
        pass_token : bool, optional
            Call ``job(*args, cancel=token)`` with a :class:`CancelToken`
            that Cancel sets, so the job can stop early.

        Returns
        -------
        Task
            The running task.
        """
        task = Task(message)
        task.ended.connect(lambda: self._ended(task))
        self.tasks.append(task)
        task.job = executor().submit(job, *args, done=done, failed=failed, owner=task, pass_token=pass_token)
        task.job.ended.append(lambda _job: task._end())
        self.changed.emit()
        return task

    def _ended(self, task: Task) -> None:
        if task in self.tasks:
            self.tasks.remove(task)
            self.changed.emit()

    def cancel_last(self) -> None:
        """Cancel the task started last."""
        if self.tasks:
            self.tasks[-1].cancel()
