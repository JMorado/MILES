"""Opening the window from Python: miles.launch in plain Python and in a Jupyter kernel."""

import contextlib
import json
import os
import sys
import textwrap
import time
from pathlib import Path
from typing import Any

import pytest
from PySide6.QtCore import QSettings, QTimer
from PySide6.QtWidgets import QApplication

import miles
from miles import interactive
from miles.api import App, Picture
from miles.files.common import FileFormatError
from miles.ui.mainwindow import MainWindow

EXAMPLES = Path(__file__).parent.parent / "examples" / "qm"
FORMALDEHYDE = str(EXAMPLES / "formaldehyde.molden")
NITROGEN_DIOXIDE = str(EXAMPLES / "nitrogen_dioxide.molden")


def close_window(window: MainWindow) -> None:
    for canvas in window.canvases:
        canvas.doc.modified = False
    window.close()


@pytest.fixture
def fresh(app, tmp_path, monkeypatch):
    """No MILES window open; the exception hook and quit setting are restored afterwards."""
    for widget in QApplication.topLevelWidgets():
        if isinstance(widget, MainWindow) and widget.isVisible():
            close_window(widget)
    monkeypatch.setattr(sys, "excepthook", sys.excepthook)
    quits = app.quitOnLastWindowClosed()
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    yield settings
    if interactive._launched is not None:
        close_window(interactive._launched)
    interactive._launched = None
    app.setQuitOnLastWindowClosed(quits)


def test_launch_opens_files_and_returns_the_console_app(fresh):
    app = miles.launch(FORMALDEHYDE, settings=fresh, appearance="dark")
    assert isinstance(app, App) and app.window.isVisible()
    assert interactive.running_app() is app
    assert app.current is not None and app.current.name == "formaldehyde.molden"
    assert app.window.appearance == "dark" and fresh.value("appearance") is None   # not remembered
    assert not QApplication.instance().quitOnLastWindowClosed()

    again = miles.launch(NITROGEN_DIOXIDE)
    assert again is app and len(app.documents) == 2


def test_launch_checks_paths_and_names_first(fresh):
    with pytest.raises(FileNotFoundError):
        miles.launch("no-such-file.xyz", settings=fresh)
    with pytest.raises(ValueError):
        miles.launch(appearance="sepia", settings=fresh)
    assert interactive.running_app() is None   # nothing was opened


def test_open_raises_instead_of_showing_a_dialog(fresh, tmp_path):
    app = miles.launch(settings=fresh)
    bad = tmp_path / "notes.unknownformat"
    bad.write_text("not a structure\n")
    with pytest.raises(FileFormatError):
        app.open(str(bad))
    assert QApplication.activeModalWidget() is None


def test_image_is_shown_inline_at_half_its_pixel_size(fresh, tmp_path):
    app = miles.launch(FORMALDEHYDE, settings=fresh)
    doc = app.current
    assert doc is not None
    picture = doc.image()
    assert isinstance(picture, Picture)
    data, size = picture._repr_png_()
    assert data.startswith(b"\x89PNG\r\n\x1a\n")
    assert size == {"width": round(picture.width / 2), "height": round(picture.height / 2)}
    whole = doc.image(crop=False)
    assert whole.width > picture.width and whole.height >= picture.height   # cropping trims the border
    path = tmp_path / "view.png"
    picture.save(str(path))
    assert path.read_bytes().startswith(b"\x89PNG")


def test_wait_until_closed_returns_when_the_window_closes(fresh):
    app = miles.launch(settings=fresh)
    QTimer.singleShot(100, lambda: close_window(app.window))
    start = time.monotonic()
    app.wait_until_closed()
    assert not app.window.isVisible() and time.monotonic() - start < 10
    assert interactive.running_app() is None

    newer = miles.launch(settings=fresh)   # a closed window is replaced
    assert newer is not app and newer.window.isVisible()


def test_plain_prompt_gets_a_note_and_errors_go_to_the_dialog(fresh, monkeypatch, capsys):
    from miles.ui.errors import ErrorReporter

    monkeypatch.setattr(sys, "ps1", ">>> ", raising=False)
    miles.launch(settings=fresh)
    assert "wait_until_closed" in capsys.readouterr().err
    assert isinstance(sys.excepthook, ErrorReporter)


# -- a real Jupyter kernel ------------------------------------------------------------------------

def _kernel(tmp_path: Path) -> Any:
    """Start an IPython kernel of this Python, offscreen, with settings and folders of its own."""
    from jupyter_client.kernelspec import KernelSpecManager
    from jupyter_client.manager import KernelManager

    spec = tmp_path / "kernels" / "miles-test"
    spec.mkdir(parents=True)
    spec.joinpath("kernel.json").write_text(json.dumps({
        "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "MILES test", "language": "python"}))
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", MILES_SOFTWARE="1",
               MILES_PLUGIN_DIR=str(tmp_path / "no-plugins"),
               MILES_FRAGMENT_DIR=str(tmp_path / "no-fragments"),
               MILES_AUTOSAVE_DIR=str(tmp_path / "autosave"))
    manager = KernelManager(kernel_name="miles-test",
                            kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(spec.parent)]))
    manager.start_kernel(env=env, cwd=str(tmp_path))
    client = manager.client()
    client.start_channels()
    client.wait_for_ready(timeout=60)
    return manager, client


def _run(client: Any, code: str, timeout: float = 60) -> dict[str, Any]:
    """Run a cell; return its printed output, result, rich outputs and error."""
    out: dict[str, Any] = {"stdout": "", "stderr": "", "data": {}, "metadata": {}, "error": None}

    def collect(msg: dict[str, Any]) -> None:
        kind, content = msg["msg_type"], msg["content"]
        if kind == "stream":
            out[content["name"]] += content["text"]
        elif kind in ("execute_result", "display_data"):
            out["data"].update(content["data"])
            out["metadata"].update(content.get("metadata", {}))
        elif kind == "error":
            out["error"] = "\n".join(content["traceback"])

    reply = client.execute_interactive(textwrap.dedent(code), output_hook=collect, timeout=timeout)
    assert reply["content"]["status"] == "ok", out["error"]
    return out


def _later_output(client: Any, wanted: list[str], timeout: float = 20) -> str:
    """Printed output sent between cells (by the window's event loop), until every wanted text came."""
    import queue

    text = ""
    end = time.monotonic() + timeout
    while time.monotonic() < end and not all(w in text for w in wanted):
        try:
            msg = client.get_iopub_msg(timeout=0.5)
        except queue.Empty:
            continue
        if msg["msg_type"] == "stream":
            text += msg["content"]["text"]
    return text


def test_launch_in_a_jupyter_kernel(tmp_path):
    pytest.importorskip("ipykernel")
    pytest.importorskip("jupyter_client")
    manager, client = _kernel(tmp_path)
    try:
        settings = str(tmp_path / "settings.ini")
        out = _run(client, f"""
            from PySide6.QtCore import QSettings
            import miles
            app = miles.launch({FORMALDEHYDE!r},
                                    settings=QSettings({settings!r}, QSettings.Format.IniFormat))
            print(app.current.name, len(app.documents), get_ipython().active_eventloop)
            """)
        assert out["stdout"].split() == ["formaldehyde.molden", "1", "qt"]

        # the event loop runs between cells: a timer set in one cell fires before the next
        _run(client, """
            from PySide6.QtCore import QTimer
            ticks = []
            QTimer.singleShot(50, lambda: ticks.append(1))
            """)
        time.sleep(0.5)
        assert _run(client, "print(ticks, app.window.isVisible())")["stdout"].strip() == "[1] True"

        out = _run(client, """
            doc = app.current
            oxygen = doc.select(smarts="[#8]")
            before = len(doc.molecule)
            doc.edit(lambda m: m.remove_atoms(oxygen))
            print(oxygen, before, len(doc.molecule))
            doc.undo()
            print(len(doc.molecule), doc.selection)
            """)
        assert out["stdout"].split("\n")[:2] == ["[1] 4 3", "4 [1]"]   # undo brings the selection back too

        out = _run(client, "doc.image()")
        assert out["data"]["image/png"] and out["metadata"]["image/png"]["width"] > 10

        out = _run(client, f"""
            again = miles.launch({NITROGEN_DIOXIDE!r})
            print(again is app, len(app.documents), app.current.name)
            """)
        assert out["stdout"].split() == ["True", "2", "nitrogen_dioxide.molden"]

        # an error inside the window between cells: traceback printed and the error dialog shown
        _run(client, """
            from PySide6.QtWidgets import QApplication

            def broken():
                raise ZeroDivisionError("raised inside the window")

            def answer_dialog():
                box = QApplication.activeModalWidget()
                print("dialog:", box.text() if box is not None else None)
                if box is not None:
                    box.done(0)

            QTimer.singleShot(100, broken)
            QTimer.singleShot(800, answer_dialog)
            """)
        later = _later_output(client, ["raised inside the window", "dialog: ZeroDivisionError"])
        assert "Traceback" in later and "dialog: ZeroDivisionError: raised inside the window" in later

        # while a cell runs, a failing command shows in the cell output as IPython's own error
        failed = client.execute_interactive("app.run('No Such Command')", timeout=60)
        assert failed["content"]["status"] == "error" and failed["content"]["ename"] == "KeyError"
        assert _run(client, "print(QApplication.activeModalWidget())")["stdout"].strip() == "None"
    finally:
        with contextlib.suppress(Exception):   # the kernel may be gone already; shutting down still applies
            _run(client, "[c.doc.__setattr__('modified', False) for c in app.window.canvases]", timeout=20)
        client.stop_channels()
        manager.shutdown_kernel(now=False)
