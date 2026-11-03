"""Shared fixtures. Qt runs offscreen and the canvas draws with QPainter."""

import atexit
import os
import shutil
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["MILES_SOFTWARE"] = "1"
# A plugin folder that does not exist, so plugins installed for daily use do not run in the tests
os.environ["MILES_PLUGIN_DIR"] = os.path.join(tempfile.gettempdir(), "miles-tests-without-plugins")
os.environ["MILES_FRAGMENT_DIR"] = os.path.join(tempfile.gettempdir(), "miles-tests-no-fragments")
# A fresh autosave folder for each run, so the tests neither see nor leave documents to restore
os.environ["MILES_AUTOSAVE_DIR"] = tempfile.mkdtemp(prefix="miles-tests-autosave-")

import pytest
from PySide6.QtWidgets import QApplication

from miles.testing import miles_app, qt_application  # noqa: F401  (the fixtures plugins test with)

atexit.register(shutil.rmtree, os.environ["MILES_AUTOSAVE_DIR"], True)


@pytest.fixture(scope="session")
def app() -> QApplication:
    instance = QApplication.instance()
    return instance if isinstance(instance, QApplication) else QApplication([])
