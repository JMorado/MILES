"""MILES as an application of its own: the macOS bundle, the Linux desktop entry, the Windows icon."""

import configparser
import os
import plistlib
import shutil
import struct
import sys

import pytest
from PySide6.QtWidgets import QApplication

from miles.ui import desktop
from miles.ui.app import parse_args


def test_install_option_is_parsed() -> None:
    assert parse_args(["--install-app"]).install_app == ""
    assert parse_args(["--install-app", "/tmp/apps"]).install_app == "/tmp/apps"
    assert parse_args([]).install_app is None


def test_process_is_named_for_the_desktop(app: QApplication) -> None:
    desktop.identify_process(app)
    assert app.desktopFileName() == desktop.SLUG


def test_logo_png_has_the_size_asked_for(app: QApplication) -> None:
    data = desktop.logo_png(48)
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert struct.unpack(">II", data[16:24]) == (48, 48)


def test_windows_icon_holds_every_size(app: QApplication, tmp_path) -> None:
    path = tmp_path / "miles.ico"
    desktop.write_ico(path)
    data = path.read_bytes()
    reserved, kind, count = struct.unpack("<HHH", data[:6])
    assert (reserved, kind, count) == (0, 1, len(desktop.ICO_SIZES))
    for k, size in enumerate(desktop.ICO_SIZES):
        entry = data[6 + 16 * k:22 + 16 * k]
        width, height, _, _, planes, bits, length, offset = struct.unpack("<BBBBHHII", entry)
        assert width == height == (0 if size >= 256 else size)
        assert (planes, bits) == (1, 32)
        assert data[offset:offset + 8] == b"\x89PNG\r\n\x1a\n"
        assert struct.unpack(">II", data[offset + 16:offset + 24]) == (size, size)
        assert offset + length <= len(data)


def test_linux_desktop_entry(app: QApplication, tmp_path) -> None:
    python = "/opt/envs/my miles/bin/python"
    entry, icon = desktop.install_linux(tmp_path, python)
    assert icon == tmp_path / "icons" / "hicolor" / "scalable" / "apps" / "miles.svg"
    assert icon.read_text(encoding="utf-8").startswith("<svg")
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str  # type: ignore[assignment,method-assign]   # keys are case-sensitive
    parser.read(entry)
    fields = parser["Desktop Entry"]
    assert fields["Name"] == "MILES" and fields["Icon"] == "miles" and fields["Type"] == "Application"
    # the path with a space is quoted
    assert fields["Exec"] == '"/opt/envs/my miles/bin/python" -P -m miles %F'
    assert fields["StartupWMClass"] == desktop.SLUG


@pytest.mark.skipif(sys.platform != "darwin" or shutil.which("iconutil") is None, reason="macOS only")
def test_macos_bundle(app: QApplication, tmp_path) -> None:
    bundle = desktop.install_macos(tmp_path, sys.executable)
    assert bundle == tmp_path / "MILES.app"
    contents = bundle / "Contents"
    info = plistlib.loads((contents / "Info.plist").read_bytes())
    assert info["CFBundleName"] == "MILES" and info["CFBundleExecutable"] == "MILES"
    assert info["CFBundleIdentifier"] == desktop.IDENTIFIER
    launcher = contents / "MacOS" / "MILES"
    assert os.access(launcher, os.X_OK)
    assert '"$here/python" -P -m miles "$@"' in launcher.read_text(encoding="utf-8")
    assert os.readlink(contents / "MacOS" / "python") == sys.executable
    assert (contents / "Resources" / "MILES.icns").read_bytes()[:4] == b"icns"
    again = desktop.install_macos(tmp_path, sys.executable)   # making it again replaces the old one
    assert again == bundle and [p.name for p in tmp_path.iterdir()] == ["MILES.app"]
