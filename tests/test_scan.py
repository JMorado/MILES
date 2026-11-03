"""Relaxed coordinate scans and the energies that travel with frames."""

from pathlib import Path

import numpy as np
import pytest

from miles.builtin_plugins.scan.relaxed_scan import Coordinate, coordinate_kind, measure, scan
from miles.core.chem import from_smiles
from miles.core.forcefield import clean_geometry
from miles.files import session
from miles.files.common import JobInfo
from miles.files.openbabel import HARTREE, read, xyz_comment_energies

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def butane():
    mol, _ = from_smiles("CCCC")
    clean_geometry(mol)
    return mol


def test_coordinate_kind_and_measure(butane):
    assert coordinate_kind([0, 1]) == Coordinate.DISTANCE
    assert coordinate_kind([0, 1, 2]) == Coordinate.ANGLE
    assert coordinate_kind([0, 1, 2, 3]) == Coordinate.DIHEDRAL
    with pytest.raises(ValueError):
        coordinate_kind([0])
    with pytest.raises(ValueError):
        coordinate_kind([0, 1, 1])
    assert measure(butane.coords, [0, 1]) == pytest.approx(1.52, abs=0.03)
    assert 100 < measure(butane.coords, [0, 1, 2]) < 120


def test_dihedral_scan_of_butane(butane):
    before = butane.coords.copy()
    result = scan(butane, [0, 1, 2, 3], -180.0, 180.0, 12)
    assert np.allclose(butane.coords, before)   # the input is not modified
    assert result.force_field == "MMFF94"
    assert result.frames.shape == (13, len(butane), 3)
    assert np.isfinite(result.energies).all()
    targets = np.linspace(-180.0, 180.0, 13)
    offsets = (result.values - targets + 180.0) % 360.0 - 180.0
    assert np.abs(offsets).max() < 1.0
    relative = result.energies - result.energies.min()
    anti, eclipsed = relative[[0, 12]], relative[6]
    assert anti.max() < 0.1   # anti is the minimum
    assert 3.0 < eclipsed < 8.0   # the syn (eclipsed) form is the maximum
    assert relative.argmax() == 6
    gauche = relative[[4, 8]]
    assert np.all((gauche > 0.3) & (gauche < 1.5))


def test_distance_and_angle_scans(butane):
    result = scan(butane, [0, 1], 1.45, 1.65, 4)
    assert np.allclose(result.values, np.linspace(1.45, 1.65, 5), atol=0.005)
    assert result.energies.argmin() in (1, 2)
    result = scan(butane, [0, 1, 2], 105.0, 125.0, 2)
    assert np.allclose(result.values, [105.0, 115.0, 125.0], atol=0.5)


def test_scan_rejects_bad_input(butane):
    with pytest.raises(ValueError):
        scan(butane, [0, 1], 1.4, 1.6, 0)
    with pytest.raises(ValueError):
        scan(butane, [0, 99], 1.4, 1.6, 2)
    with pytest.raises(ValueError):
        scan(butane, [0, 1, 2], 10.0, 190.0, 2)


def test_gaussian_optimisation_energies():
    mol, info = read(str(DATA / "water_opt.log"))
    assert info.frames is not None and len(info.frames) == 3
    assert info.energies is not None
    assert np.allclose(info.energies, np.array([-76.4080123, -76.4088456, -76.4089510]) * HARTREE, atol=1e-3)
    assert np.allclose(mol.coords, info.frames[-1])


def write_xyz(path: Path, comments: list[str]) -> None:
    lines = []
    for k, comment in enumerate(comments):
        lines += ["3", comment, f"O 0.0 0.0 {0.117 + 0.001 * k:.4f}",
                  "H 0.0 0.757 -0.469", "H 0.0 -0.757 -0.469"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.mark.parametrize("comments, first", [
    ([" energy: -5.070544440612 gnorm: 0.000377 xtb: 6.6.1", " energy: -5.071 gnorm: 0.0001 xtb: 6.6.1"],
     -5.070544440612 * HARTREE),
    (["Coordinates from ORCA-job water E -76.326270808", "Coordinates from ORCA-job water E -76.33"],
     -76.326270808 * HARTREE),
    (['Properties=species:S:1:pos:R:3 energy=-14.25 pbc="F F F"',
      "Properties=species:S:1:pos:R:3 energy=-14.3"], -14.25 * 23.060548),
    (["   -22.73812837", "   -22.73800000"], -22.73812837 * HARTREE),
])
def test_xyz_comment_energies(tmp_path, comments, first):
    path = tmp_path / "frames.xyz"
    write_xyz(path, comments)
    _, info = read(str(path))
    assert info.energies is not None and len(info.energies) == 2
    assert info.energies[0] == pytest.approx(first)


def test_xyz_without_energies_or_mixed_comments():
    assert xyz_comment_energies(["frame 1", "frame 2"]) is None
    assert xyz_comment_energies([" energy: -5.0 gnorm: 0.1", "frame 2"]) is None
    assert xyz_comment_energies(["12", "13"]) is None   # bare integers are not energies
    _, info = read(str(Path(__file__).parents[1] / "examples" / "alanine_dipeptide_md.xyz"))
    assert info.frames is not None and info.energies is None


def test_session_keeps_energies(butane, tmp_path):
    frames = np.array([butane.coords, butane.coords + 0.1])
    path = str(tmp_path / "scan.miles")
    session.write(path, butane, JobInfo(title="butane", frames=frames, energies=[1.5, 2.5]))
    _, info = session.read(path)
    assert info.energies == [1.5, 2.5]


def test_document_drops_energies_that_do_not_fit(app, butane):
    from miles.ui.document import Document

    frames = np.array([butane.coords, butane.coords + 0.1])
    doc = Document(butane.copy(), JobInfo(frames=frames, energies=[1.0, 2.0, 3.0]))
    assert doc.frame_energies is None
    doc.set_frames(frames, energies=[1.0, 2.0])
    assert doc.frame_energies == [1.0, 2.0]
    doc.set_frames(frames)
    assert doc.frame_energies is None
