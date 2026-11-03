"""The plugin template in examples/plugins/template works, and its own tests pass."""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from miles.testing import document_from_smiles, load_plugin
from miles.ui.mainwindow import MainWindow

TEMPLATE = Path(__file__).resolve().parent.parent / "examples" / "plugins" / "template"


def test_the_template_tests_pass() -> None:
    """Run them as a plugin author would, from the template's folder."""
    environment = {k: v for k, v in os.environ.items() if not k.startswith("MILES_")}
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"], cwd=TEMPLATE,
                            env=environment, capture_output=True, text=True, timeout=300, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_template_registers_like_an_installed_plugin(app: QApplication, tmp_path: Path) -> None:
    spec = importlib.util.spec_from_file_location("miles_atom_counter",
                                                  TEMPLATE / "miles_atom_counter" / "__init__.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    window = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    try:
        info = load_plugin(window.api, module)
        assert (info.name, info.version) == ("Atom Counter", "0.1.0")
        doc = document_from_smiles(window.api, "CC")
        assert module.element_counts(doc).startswith("C 2")
        assert "Tools > Count Atoms" in window.api.commands()
    finally:
        for canvas in window.canvases:
            canvas.doc.modified = False
        window.close()
        window.deleteLater()
        app.processEvents()
