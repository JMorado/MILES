"""The data model, documents and file readers work without Qt."""

import subprocess
import sys
import textwrap


def test_core_and_files_do_not_import_qt():
    code = textwrap.dedent("""
        import importlib, pkgutil, sys
        import miles.core, miles.files
        for package in (miles.core, miles.files):
            for module in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
                importlib.import_module(module.name)
        loaded = sorted(name for name in sys.modules if name.startswith("PySide6"))
        print(",".join(loaded))
    """)
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == ""


def test_a_document_edits_and_undoes_without_a_window():
    code = textwrap.dedent("""
        import sys
        from miles.core.document import Change, Document
        from miles.files.formats import parse_text
        doc = Document(parse_text("3\\nwater\\nO 0 0 0\\nH 0.96 0 0\\nH -0.24 0.93 0\\n"))
        events = []
        doc.listeners.append(lambda event: events.append(event.kind))
        doc.edit(lambda mol: mol.remove_atoms([2]))
        doc.set_representation("tube", [0])
        assert len(doc.mol) == 2 and doc.display.rep[0] == "tube"
        doc.undo()
        doc.undo()
        assert len(doc.mol) == 3 and doc.display.rep[0] == ""
        # undo announces what its step changed
        assert events == [Change.STRUCTURE, Change.STYLE, Change.STYLE, Change.STRUCTURE]
        assert not any(name.startswith("PySide6") for name in sys.modules)
    """)
    subprocess.run([sys.executable, "-c", code], check=True)
