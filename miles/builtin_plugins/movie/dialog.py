"""Exporting movies: a turntable rotation or the frames of a trajectory.

Frames are rendered with the same code as exported pictures. MP4 and GIF
files are encoded with ``ffmpeg`` when it is installed; otherwise the frames
are written as numbered PNG files that any video tool can assemble.
"""

import math
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from enum import StrEnum
from pathlib import Path
from typing import override

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QGridLayout, QLabel, QProgressDialog, QSpinBox

from miles import rotation_matrix
from miles.plugin import App, DocumentView, PluginDialog, Severity
from miles.plugin.widgets import Segmented, muted

FPS = 30


class MovieKind(StrEnum):
    """What the movie shows."""

    TURNTABLE = "turntable"     # one full turn about the vertical axis
    TRAJECTORY = "trajectory"   # the frames of the document


class MovieFormat(StrEnum):
    """File the frames are written to."""

    MP4 = "mp4"
    GIF = "gif"
    PNG = "png"   # a folder of numbered PNG files


def ffmpeg() -> str | None:
    """Path of the ffmpeg program, if installed."""
    return shutil.which("ffmpeg")


def render_frames(doc: DocumentView, kind: MovieKind, count: int, folder: Path, scale: float,
                  background: str, ambient_occlusion: bool,
                  progress: Callable[[int], bool] | None = None) -> int:
    """Render the frames of a movie as PNG files.

    Parameters
    ----------
    doc : DocumentView
        View to film; its orientation and frame are restored afterwards.
    kind : MovieKind
        Turntable or trajectory.
    count : int
        Number of frames for a turntable (a trajectory uses all its frames).
    folder : pathlib.Path
        Where ``frame_0000.png`` and so on are written.
    scale, background, ambient_occlusion
        As for :meth:`~miles.api.DocumentView.render`.
    progress : callable, optional
        Called with the frame number; return False to stop.

    Returns
    -------
    int
        Number of frames written.
    """
    rotation, frame = doc.rotation, doc.frame
    total = doc.frame_count if kind == MovieKind.TRAJECTORY else count
    written = 0
    try:
        for k in range(total):
            if progress is not None and not progress(k):
                break
            if kind == MovieKind.TRAJECTORY:
                doc.frame = k
            else:
                doc.rotation = rotation_matrix([0, 1, 0], 2 * math.pi * k / total) @ rotation
            image = doc.render(scale, background, ambient_occlusion)
            image.save(str(folder / f"frame_{k:04d}.png"))
            written += 1
    finally:
        doc.rotation = rotation
        if kind == MovieKind.TRAJECTORY:
            doc.frame = frame
    return written


def encode(folder: Path, path: Path, movie_format: MovieFormat, fps: int = FPS) -> None:
    """Join numbered PNG frames into an MP4 or GIF file with ffmpeg.

    Raises
    ------
    RuntimeError
        If ffmpeg is missing or fails.
    """
    program = ffmpeg()
    if program is None:
        raise RuntimeError("ffmpeg is not installed")
    frames = str(folder / "frame_%04d.png")
    if movie_format == MovieFormat.MP4:
        args = ["-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2,format=yuv420p", "-c:v", "libx264", "-crf", "18"]
    else:
        args = ["-vf", "split[a][b];[a]palettegen[p];[b][p]paletteuse", "-loop", "0"]
    command = [program, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", frames, *args, str(path)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "ffmpeg failed")


class MovieDialog(PluginDialog):
    """Options for a movie of the current view; Save asks where, then renders and writes it.

    After :meth:`run` accepts, :attr:`saved_path` is where the movie went.
    """

    def __init__(self, app: App, doc: DocumentView) -> None:
        super().__init__(app, "Export Movie")
        self.doc = doc
        self.app = app
        self.saved_path = ""
        self.body.setSpacing(12)
        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(10)

        self.kind = Segmented(["Turntable", "Trajectory"])
        has_frames = doc.frame_count > 1
        trajectory = self.kind.group.button(1)
        if trajectory is not None:
            trajectory.setEnabled(has_frames)
        if has_frames:
            self.kind.set_index(1)
        self.seconds = QSpinBox()
        self.seconds.setRange(2, 60)
        self.seconds.setValue(6)
        self.seconds.setSuffix(" s")
        self.scale = Segmented(["1x", "2x"])
        self.background = Segmented(["Theme", "White"])
        self.background.set_index(1)
        self.ao = QCheckBox("Ambient occlusion")
        self.format = Segmented(["MP4", "GIF", "PNG frames"])
        if ffmpeg() is None:
            for n in (0, 1):
                b = self.format.group.button(n)
                if b is not None:
                    b.setEnabled(False)
                    b.setToolTip("Install ffmpeg to export videos")
            self.format.set_index(2)

        rows = [("Movie", self.kind), ("Length", self.seconds), ("Size", self.scale),
                ("Background", self.background), ("Options", self.ao), ("Format", self.format)]
        for n, (name, widget) in enumerate(rows):
            grid.addWidget(muted(QLabel(name)), n, 0)
            grid.addWidget(widget, n, 1)
        grid.setColumnStretch(2, 1)
        self.body.addLayout(grid)
        if ffmpeg() is None:
            hint = "MP4 and GIF need ffmpeg, e.g. micromamba install -c conda-forge ffmpeg."
            self.body.addWidget(muted(QLabel(hint)))
        self.buttons("Save...")

    @override
    def accept(self) -> None:
        """Ask where the movie goes, render and write it, then close (or stay open when cancelled)."""
        self._save()

    def _save(self) -> None:
        kind = MovieKind.TRAJECTORY if self.kind.index() == 1 else MovieKind.TURNTABLE
        movie_format = [MovieFormat.MP4, MovieFormat.GIF, MovieFormat.PNG][self.format.index()]
        if movie_format == MovieFormat.PNG:
            target = self.app.ask_folder("Folder for the Frames")
        else:
            file_filter = f"{movie_format.upper()} (*.{movie_format})"
            target = self.app.ask_save_path("Export Movie", f"movie.{movie_format}", file_filter)
        if not target:
            return
        count = self.seconds.value() * FPS
        total = self.doc.frame_count if kind == MovieKind.TRAJECTORY else count
        dialog = QProgressDialog("Rendering frames...", "Cancel", 0, total, self)
        dialog.setWindowModality(Qt.WindowModality.WindowModal)
        dialog.setMinimumDuration(0)

        def progress(k: int) -> bool:
            dialog.setValue(k)
            return not dialog.wasCanceled()

        scale = 2.0 if self.scale.index() == 1 else 1.0
        background = "white" if self.background.index() == 1 else "theme"
        try:
            if movie_format == MovieFormat.PNG:
                written = render_frames(self.doc, kind, count, Path(target), scale, background,
                                        self.ao.isChecked(), progress)
            else:
                with tempfile.TemporaryDirectory() as tmp:
                    written = render_frames(self.doc, kind, count, Path(tmp), scale, background,
                                            self.ao.isChecked(), progress)
                    if written == total:
                        dialog.setLabelText("Encoding...")
                        encode(Path(tmp), Path(target), movie_format)
        except RuntimeError as exc:
            dialog.close()
            self.app.message(str(exc), severity=Severity.WARNING)
            return
        dialog.close()
        if written == total:
            self.saved_path = target
            super().accept()
