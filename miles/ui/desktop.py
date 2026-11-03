"""Making MILES appear as an application of its own on macOS, Linux and Windows.

A program started with Python is named after the interpreter: the Dock,
the taskbar and the application switcher show "python" and its icon,
because each system takes an application's name and icon from somewhere
Python does not provide:

* macOS reads them from an application bundle. :func:`install` makes
  ``MILES.app`` with the logo and a launcher; ``Contents/MacOS/python`` is a
  symbolic link to the Python of this environment, so the running program
  belongs to the bundle and is called MILES, also when the launcher is
  started from a terminal.
* Linux desktops read them from a ``.desktop`` file and match the window to
  it by its application id, which :func:`identify_process` sets.
* Windows groups taskbar buttons by an application id, which
  :func:`identify_process` sets so the window's own icon is used; the Start
  menu entry made by :func:`install` gives the name and the icon.

What :func:`install` makes points at the Python it was run with; run it again
after moving or recreating that environment. Conda packages can make the
same entries when they are installed (see the release notes of the recipe).
"""

import os
import plistlib
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from miles import __version__

if TYPE_CHECKING:
    from PySide6.QtWidgets import QApplication

NAME = "MILES"
SLUG = "miles"
IDENTIFIER = "com.miles.MILES"   # also where Qt keeps the settings on macOS
SUMMARY = "Molecular Improvisation & Ligand Exploration Studio"
LAUNCH = ("-P", "-m", "miles")    # -P: the current folder does not go on the module path

ICNS_SIZES = (16, 32, 128, 256, 512)       # each also at twice the size
ICO_SIZES = (16, 24, 32, 48, 64, 128, 256)


def identify_process(app: "QApplication") -> None:
    """Give the running program the identity of MILES, before any window is shown."""
    app.setDesktopFileName(SLUG)   # Linux: the window belongs to miles.desktop
    if sys.platform == "win32":    # Windows: a taskbar group of its own, with the window's icon
        import ctypes

        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(IDENTIFIER)


def install(folder: Path | None = None, python: str | None = None) -> list[Path]:
    """Add MILES to the system's applications; returns the files and folders made.

    Parameters
    ----------
    folder : Path, optional
        Where to put the entry: the folder for ``MILES.app`` on macOS
        (``~/Applications``), the data folder on Linux (``$XDG_DATA_HOME``
        or ``~/.local/share``), the Start menu folder on Windows.
    python : str, optional
        The Python that runs MILES (the one running now).
    """
    python = python or sys.executable
    if sys.platform == "darwin":
        return [install_macos(folder or Path.home() / "Applications", python)]
    if sys.platform == "win32":
        appdata = Path(os.environ.get("APPDATA", Path.home()))
        start_menu = appdata / "Microsoft" / "Windows" / "Start Menu" / "Programs"
        return install_windows(folder or start_menu, python)
    data = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    return install_linux(folder or data, python)


# -- macOS ------------------------------------------------------------------------------------

def install_macos(folder: Path, python: str) -> Path:
    """Make ``MILES.app`` in a folder, replacing an older one; returns its path."""
    folder.mkdir(parents=True, exist_ok=True)
    bundle = folder / f"{NAME}.app"
    work = Path(tempfile.mkdtemp(prefix=".miles-", dir=folder)) / f"{NAME}.app"
    contents = work / "Contents"
    (contents / "MacOS").mkdir(parents=True)
    (contents / "Resources").mkdir()
    (contents / "Info.plist").write_bytes(plistlib.dumps({
        "CFBundleName": NAME,
        "CFBundleDisplayName": NAME,
        "CFBundleIdentifier": IDENTIFIER,
        "CFBundleExecutable": NAME,
        "CFBundleIconFile": NAME,
        "CFBundlePackageType": "APPL",
        "CFBundleShortVersionString": __version__,
        "CFBundleVersion": __version__,
        "NSHighResolutionCapable": True,
    }))
    launcher = contents / "MacOS" / NAME
    launcher.write_text(
        "#!/bin/sh\n"
        "# Start MILES with the Python next to this file, a link to the environment it was made for:\n"
        "# running it from inside the bundle is what makes macOS call the program MILES.\n"
        'here="$(dirname "$(readlink -f "$0")")"\n'
        f'exec "$here/python" {" ".join(LAUNCH)} "$@"\n', encoding="utf-8")
    launcher.chmod(0o755)
    (contents / "MacOS" / "python").symlink_to(python)
    write_icns(contents / "Resources" / f"{NAME}.icns")
    if bundle.exists():
        shutil.rmtree(bundle)
    work.rename(bundle)
    shutil.rmtree(work.parent, ignore_errors=True)
    os.utime(bundle)   # Finder and the Dock notice the new bundle
    return bundle


def write_icns(path: Path) -> None:
    """Write the logo as a macOS icon file (needs ``iconutil``, part of macOS)."""
    with tempfile.TemporaryDirectory(suffix=".iconset") as iconset:
        for size in ICNS_SIZES:
            Path(iconset, f"icon_{size}x{size}.png").write_bytes(logo_png(size))
            Path(iconset, f"icon_{size}x{size}@2x.png").write_bytes(logo_png(2 * size))
        subprocess.run(["iconutil", "-c", "icns", iconset, "-o", str(path)], check=True)


# -- Linux ------------------------------------------------------------------------------------

def install_linux(data: Path, python: str) -> list[Path]:
    """Write ``miles.desktop`` and the icon into a data folder; returns both paths."""
    from miles.ui.icons import APP_ICON

    icon = data / "icons" / "hicolor" / "scalable" / "apps" / f"{SLUG}.svg"
    icon.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(APP_ICON, icon)
    entry = data / "applications" / f"{SLUG}.desktop"
    entry.parent.mkdir(parents=True, exist_ok=True)
    command = " ".join(_desktop_quote(part) for part in (python, *LAUNCH))
    entry.write_text(
        "[Desktop Entry]\n"
        "Type=Application\n"
        f"Name={NAME}\n"
        "GenericName=Molecule Builder and Viewer\n"
        f"Comment={SUMMARY}\n"
        f"Exec={command} %F\n"
        f"Icon={SLUG}\n"
        "Terminal=false\n"
        "Categories=Science;Chemistry;Education;\n"
        f"StartupWMClass={SLUG}\n", encoding="utf-8")
    return [entry, icon]


def _desktop_quote(arg: str) -> str:
    """Quote one argument of a desktop entry's Exec line (freedesktop rules)."""
    if not any(c in arg for c in ' \t\n"\'\\><~|&;$*?#()`'):
        return arg
    escaped = "".join("\\" + c if c in '"`$\\' else c for c in arg)
    return f'"{escaped}"'


# -- Windows ----------------------------------------------------------------------------------

def install_windows(start_menu: Path, python: str) -> list[Path]:
    """Make a Start menu shortcut with the logo; returns the shortcut and the icon file."""
    icon = Path(os.environ.get("LOCALAPPDATA", Path.home())) / NAME / f"{SLUG}.ico"
    icon.parent.mkdir(parents=True, exist_ok=True)
    write_ico(icon)
    windowed = Path(python).with_name("pythonw.exe")   # no console window
    target = windowed if windowed.exists() else Path(python)
    start_menu.mkdir(parents=True, exist_ok=True)
    shortcut = start_menu / f"{NAME}.lnk"
    script = ("$s = (New-Object -ComObject WScript.Shell).CreateShortcut($env:MILES_LINK); "
              "$s.TargetPath = $env:MILES_TARGET; $s.Arguments = $env:MILES_ARGS; "
              "$s.IconLocation = $env:MILES_ICON; $s.Description = $env:MILES_SUMMARY; $s.Save()")
    env = {**os.environ, "MILES_LINK": str(shortcut), "MILES_TARGET": str(target),
           "MILES_ARGS": " ".join(LAUNCH), "MILES_ICON": str(icon), "MILES_SUMMARY": SUMMARY}
    subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", script], check=True, env=env)
    return [shortcut, icon]


def write_ico(path: Path) -> None:
    """Write the logo as a Windows icon file, one PNG image per size."""
    images = [logo_png(size) for size in ICO_SIZES]
    header = struct.pack("<HHH", 0, 1, len(images))
    offset = len(header) + 16 * len(images)
    directory = b""
    for size, data in zip(ICO_SIZES, images, strict=True):
        side = 0 if size >= 256 else size   # 0 means 256 in the directory
        directory += struct.pack("<BBBBHHII", side, side, 0, 0, 1, 32, len(data), offset)
        offset += len(data)
    path.write_bytes(header + directory + b"".join(images))


# -- the logo as an image ---------------------------------------------------------------------

def logo_png(size: int) -> bytes:
    """Draw the logo at a size in pixels and return it as PNG data."""
    from PySide6.QtCore import QBuffer, QByteArray, QIODevice, QRectF
    from PySide6.QtGui import QImage, QImageWriter, QPainter
    from PySide6.QtSvg import QSvgRenderer

    from miles.ui.icons import APP_ICON

    image = QImage(size, size, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    QSvgRenderer(str(APP_ICON)).render(painter, QRectF(0, 0, size, size))
    painter.end()
    data = QByteArray()
    buffer = QBuffer(data)
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    QImageWriter(buffer, QByteArray(b"png")).write(image)
    return bytes(data.data())
