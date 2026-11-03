import os
import sys

import miles


def _prefix(tmp_path, *parts):
    folder = tmp_path.joinpath("share", "openbabel", *parts)
    folder.mkdir(parents=True)
    (folder / "mmff94.ff").write_text("", encoding="utf-8")
    return folder


def test_data_folder_is_found_when_the_environment_does_not_set_it(tmp_path, monkeypatch):
    windows = _prefix(tmp_path / "windows")                 # conda on Windows: share/openbabel
    unix = _prefix(tmp_path / "unix", "3.2.1")              # conda elsewhere: share/openbabel/<version>
    for prefix, expected in ((tmp_path / "windows", windows), (tmp_path / "unix", unix)):
        monkeypatch.delenv("BABEL_DATADIR", raising=False)
        monkeypatch.setattr(sys, "prefix", str(prefix))
        miles._find_openbabel_data()
        assert os.environ["BABEL_DATADIR"] == str(expected)


def test_a_data_folder_set_by_the_environment_is_kept(tmp_path, monkeypatch):
    _prefix(tmp_path)
    monkeypatch.setenv("BABEL_DATADIR", "/chosen/by/the/user")
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    miles._find_openbabel_data()
    assert os.environ["BABEL_DATADIR"] == "/chosen/by/the/user"


def test_nothing_is_set_without_open_babel_data(tmp_path, monkeypatch):
    monkeypatch.delenv("BABEL_DATADIR", raising=False)
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    miles._find_openbabel_data()
    assert "BABEL_DATADIR" not in os.environ
