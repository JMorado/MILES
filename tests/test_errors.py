"""Unexpected errors are shown in a dialog and the application keeps running."""

import sys

from PySide6.QtWidgets import QMessageBox

from miles.ui.errors import ErrorReporter, install


def raise_and_report(reporter, error):
    try:
        raise error
    except type(error):
        reporter(*sys.exc_info())


def test_errors_are_shown_once_at_a_time(app, monkeypatch, capsys):
    shown = []
    monkeypatch.setattr(QMessageBox, "exec", lambda box: shown.append((box.text(), box.detailedText())) or 0)
    reporter = ErrorReporter(lambda: None)
    raise_and_report(reporter, ValueError("bad value"))
    assert reporter.shown == 1 and shown[0][0] == "ValueError: bad value"
    assert "Traceback" in shown[0][1] and "bad value" in shown[0][1]
    assert "ValueError: bad value" in capsys.readouterr().err   # still printed for the terminal

    reporter._showing = True   # while a dialog is open, further errors are only printed
    raise_and_report(reporter, RuntimeError("again"))
    assert reporter.shown == 1
    reporter._showing = False
    raise_and_report(reporter, KeyboardInterrupt())
    assert reporter.shown == 1


def test_install_sets_the_exception_hook(app, monkeypatch):
    monkeypatch.setattr(sys, "excepthook", sys.excepthook)
    reporter = install(lambda: None)
    assert sys.excepthook is reporter


def test_a_repeating_error_is_shown_once(app, monkeypatch):
    monkeypatch.setattr(QMessageBox, "exec", lambda box: 0)
    reporter = ErrorReporter(lambda: None)

    def paint():
        raise RuntimeError("drawing failed")

    for _ in range(3):   # the same error on every repaint
        try:
            paint()
        except RuntimeError:
            reporter(*sys.exc_info())
    assert reporter.shown == 1
