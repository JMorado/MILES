"""Aligning the frames of a trajectory."""

from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings

from miles.core.geometry import kabsch
from miles.core.trajectory import align_frames, rmsd_series
from miles.ui.mainwindow import MainWindow

EXAMPLE = Path(__file__).parent.parent / "examples" / "alanine_dipeptide_md.xyz"


@pytest.fixture
def window(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    win.resize(1000, 700)
    win.show()
    app.processEvents()
    yield win
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def tumbling_frames(count=8, atoms=12, seed=3):
    """A shape that changes a little while it turns and drifts."""
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(atoms, 3)) * 2.0
    frames = []
    for _ in range(count):
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        frames.append((base + rng.normal(size=base.shape) * 0.1) @ q.T + rng.normal(size=3) * 5.0)
    return np.array(frames)


def distances(frames):
    return np.linalg.norm(frames[:, :, None] - frames[:, None], axis=-1)


def test_aligned_frames_keep_their_shape_and_fit_best():
    frames = tumbling_frames()
    aligned = align_frames(frames, reference=2)
    assert np.allclose(aligned[2], frames[2])
    assert np.allclose(distances(aligned), distances(frames))
    plain = np.sqrt(((aligned - aligned[2]) ** 2).sum(axis=-1).mean(axis=-1))
    assert np.allclose(plain, rmsd_series(frames, reference=2), atol=1e-6)
    ref = frames[2] - frames[2].mean(axis=0)
    for frame, moved in zip(frames, aligned, strict=True):   # same as fitting each frame on its own
        p = frame - frame.mean(axis=0)
        expected = p @ kabsch(p, ref).T + frames[2].mean(axis=0)
        assert np.allclose(moved, expected, atol=1e-8)


def test_alignment_on_some_atoms_moves_every_atom():
    frames = tumbling_frames()
    aligned = align_frames(frames, atoms=[0, 1, 2, 3, 4], reference=0)
    assert np.allclose(distances(aligned), distances(frames))
    assert np.allclose(aligned[:, :5].mean(axis=1), frames[0, :5].mean(axis=0))


def test_align_frames_command_and_undo(window):
    window.open_path(str(EXAMPLE))
    canvas = window.current_canvas()
    doc = canvas.doc
    doc.set_frame(50)
    before = doc.frames.copy()
    shown = canvas.mol.coords.copy()
    window.api.run("Tools > Align Frames")
    assert doc.frame == 50
    assert np.allclose(canvas.mol.coords, shown)   # the frame on screen is the reference
    assert np.allclose(distances(doc.frames), distances(before), atol=1e-6)
    centres = doc.frames.mean(axis=1)
    assert np.ptp(centres, axis=0).max() < np.ptp(before.mean(axis=1), axis=0).max()
    window.current_canvas().commands.undo()
    assert np.allclose(doc.frames, before)
    window.current_canvas().commands.redo()
    assert not np.allclose(doc.frames, before)


def test_align_frames_needs_frames(window):
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CCO")
    window.open_structure(mol, window.current_doc().info)
    window.api.run("Align Frames")   # a single structure: a message, nothing changes
    assert window.current_doc().frames is None


def test_undoing_a_conformer_search_restores_the_trajectory(window):
    window.open_path(str(EXAMPLE))
    canvas = window.current_canvas()
    doc = canvas.doc
    before = doc.frames.copy()
    window.api.run("Generate Conformers")
    window.api.wait()   # it works in the background
    assert doc.frames.shape[0] != before.shape[0]
    window.current_canvas().commands.undo()
    assert doc.frames.shape == before.shape and np.allclose(doc.frames, before)
