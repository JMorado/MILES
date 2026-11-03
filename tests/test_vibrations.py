from pathlib import Path

import numpy as np
import pytest

from miles.builtin_plugins.vibrations.modes import ir_spectrum, mode_frames
from miles.core.vibrations import VibrationalMode
from miles.files import session
from miles.files.formats import read_file

DATA = Path(__file__).parent / "data"


def test_gaussian_frequencies_are_read():
    mol, info = read_file(str(DATA / "water_freq.log"))
    assert mol.symbols == ["O", "H", "H"]
    modes = info.vibrations
    assert [round(m.frequency, 4) for m in modes] == [1713.4237, 3727.2973, 3849.0122]
    assert [round(m.intensity, 4) for m in modes] == [75.2735, 1.5891, 17.9811]
    for mode in modes:
        assert mode.displacements.shape == (3, 3)
    assert np.allclose(modes[0].displacements[1], [0.0, -0.42, -0.56])
    assert not any(m.imaginary for m in modes)


MOLDEN_MODES = """[FREQ]
1600.0
3700.0
3800.0
[FR-COORD]
O 0.0 0.0 0.22
H 0.0 1.43 -0.88
H 0.0 -1.43 -0.88
[FR-NORM-COORD]
vibration 1
0.0 0.0 -0.07
0.0 -0.42 0.56
0.0 0.42 0.56
vibration 2
0.0 0.0 0.05
0.0 0.58 -0.40
0.0 -0.58 -0.40
vibration 3
0.0 0.07 0.0
0.0 -0.56 0.43
0.0 -0.56 -0.43
"""
MOLDEN_INTENSITIES = "[INT]\n70.0\n5.0\n40.0\n"


@pytest.mark.parametrize("intensities_first", [True, False])
def test_molden_frequencies_and_intensities_are_read(tmp_path, intensities_first):
    freq, rest = MOLDEN_MODES.split("[FR-COORD]")
    if intensities_first:   # the order Open Babel reads by itself
        text = "[Molden Format]\n" + freq + MOLDEN_INTENSITIES + "[FR-COORD]" + rest
    else:                   # the order many programs write
        text = "[Molden Format]\n" + MOLDEN_MODES + MOLDEN_INTENSITIES
    path = tmp_path / "water.molden"
    path.write_text(text, encoding="utf-8")
    mol, info = read_file(str(path))
    assert mol.symbols == ["O", "H", "H"]
    assert [m.frequency for m in info.vibrations] == [1600.0, 3700.0, 3800.0]
    assert [m.intensity for m in info.vibrations] == [70.0, 5.0, 40.0]


def test_mode_frames_loop_with_the_requested_amplitude():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    mode = VibrationalMode(1000.0, 1.0, np.array([[0.0, 0.0, 0.5], [0.0, 0.0, -1.0]]))
    frames = mode_frames(coords, mode, amplitude=0.3, count=8)
    assert frames.shape == (8, 2, 3)
    assert np.allclose(frames[0], coords)
    displacement = frames - coords
    assert np.isclose(np.abs(displacement[:, 1, 2]).max(), 0.3)
    assert np.isclose(np.abs(displacement[:, 0, 2]).max(), 0.15)
    assert np.allclose(displacement[2], -displacement[6])


def test_mode_frames_of_a_null_mode_stay_still():
    coords = np.ones((2, 3))
    frames = mode_frames(coords, VibrationalMode(10.0, None, np.zeros((2, 3))), count=4)
    assert np.allclose(frames, coords)


def test_ir_spectrum_peaks_at_the_bands():
    still = np.zeros((1, 3))
    modes = [VibrationalMode(1700.0, 100.0, still), VibrationalMode(3000.0, 50.0, still),
             VibrationalMode(-200.0, 500.0, still), VibrationalMode(2500.0, None, still)]
    x = np.arange(0.0, 4000.0, 1.0)
    y = ir_spectrum(modes, x, width=20.0)
    assert y.max() == pytest.approx(1.0)
    assert x[np.argmax(y)] == 1700.0
    assert y[3000] == pytest.approx(0.5, abs=0.01)
    assert y[1710] == pytest.approx(0.5, abs=0.01)   # half maximum at half the width
    assert y[200] < 1e-3 and y[2500] < 0.01


def test_spectrum_without_intensities_is_zero():
    modes = [VibrationalMode(1000.0, None, np.zeros((1, 3)))]
    assert not ir_spectrum(modes, np.linspace(0, 2000, 50)).any()


def test_session_keeps_vibrations(tmp_path):
    mol, info = read_file(str(DATA / "water_freq.log"))
    path = tmp_path / "water.miles"
    session.write(str(path), mol, info)
    _, restored = session.read(str(path))
    assert [m.frequency for m in restored.vibrations] == pytest.approx([m.frequency for m in info.vibrations])
    assert [m.intensity for m in restored.vibrations] == pytest.approx([m.intensity for m in info.vibrations])
    assert np.allclose(restored.vibrations[2].displacements, info.vibrations[2].displacements, atol=1e-6)
