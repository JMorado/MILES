"""Slow work on the application's executor: results on the interface thread, dropped when no longer wanted."""

import gc
import threading
import time
import weakref

import shiboken6
from PySide6.QtCore import QObject
from PySide6.QtWidgets import QApplication

from miles.ui.background import Background
from miles.ui.tasks import TaskRunner, executor


def wait_for_executor(app: QApplication, timeout: float = 5.0) -> None:
    """Process events until every submitted job has been delivered."""
    end = time.monotonic() + timeout
    while executor().pending and time.monotonic() < end:
        time.sleep(0.005)
        app.processEvents()
    app.processEvents()
    assert not executor().pending


def _run_and_delete_owner(app: QApplication,
                          func: object) -> tuple[list[object], list[threading.ExceptHookArgs]]:
    errors: list[threading.ExceptHookArgs] = []
    delivered: list[object] = []
    previous = threading.excepthook
    threading.excepthook = errors.append
    try:
        started = threading.Event()
        release = threading.Event()

        def work() -> object:
            started.set()
            release.wait(5)
            return func() if callable(func) else func

        worker = Background()
        worker.finished.connect(delivered.append)
        worker.failed.connect(delivered.append)
        worker.run(work)
        assert started.wait(5)
        shiboken6.delete(worker)   # the panel or window that owned it is gone
        release.set()
        wait_for_executor(app)
    finally:
        threading.excepthook = previous
    return delivered, errors


def test_result_after_owner_deleted_is_dropped(app: QApplication) -> None:
    assert _run_and_delete_owner(app, 42) == ([], [])


def test_error_after_owner_deleted_is_dropped(app: QApplication) -> None:
    def fail() -> None:
        raise ValueError("no grid")

    assert _run_and_delete_owner(app, fail) == ([], [])


def test_result_arrives_on_the_interface_thread(app: QApplication) -> None:
    worker = Background()
    results: list[tuple[object, str]] = []
    worker.finished.connect(lambda result: results.append((result, threading.current_thread().name)))
    worker.run(lambda: 7)
    wait_for_executor(app)
    assert results == [(7, threading.main_thread().name)]


def test_a_newer_request_replaces_an_older_one(app: QApplication) -> None:
    worker = Background()
    results: list[int] = []
    worker.finished.connect(results.append)
    release = threading.Event()
    worker.run(lambda: (release.wait(5), 1)[1])
    worker.run(lambda: 2)
    release.set()
    wait_for_executor(app)
    assert results == [2]


def test_an_owner_dropped_while_its_job_runs_goes_on_the_interface_thread(app: QApplication) -> None:
    finalised: list[str] = []

    class Probe(Background):
        def __del__(self) -> None:
            finalised.append(threading.current_thread().name)

    release = threading.Event()
    probe = Probe()
    probe.run(lambda: (release.wait(5), 42)[1])
    alive = weakref.ref(probe)
    del probe
    gc.collect()
    release.set()
    wait_for_executor(app)
    gc.collect()
    assert alive() is None and finalised == [threading.main_thread().name]


def test_cancel_stops_a_job_that_checks_its_token(app: QApplication) -> None:
    runner = TaskRunner()
    steps: list[int] = []
    results: list[object] = []
    started = threading.Event()

    def job(cancel: object) -> str:
        started.set()
        for k in range(500):
            if cancel.cancelled:   # type: ignore[attr-defined]
                return "stopped"
            steps.append(k)
            time.sleep(0.002)
        return "finished"

    task = runner.start("Counting", job, done=results.append, pass_token=True)
    assert started.wait(5) and runner.tasks == [task]
    task.cancel()
    assert not runner.tasks and not task.running   # the busy bar hides at once
    wait_for_executor(app)
    assert results == [] and len(steps) < 500   # stopped early, and nothing delivered


def test_jobs_run_on_the_worker_threads(app: QApplication) -> None:
    owner = QObject()
    seen: list[str] = []
    executor().submit(lambda: threading.current_thread().name, owner=owner, done=seen.append)
    wait_for_executor(app)
    assert seen and seen[0].startswith("miles-worker-")
