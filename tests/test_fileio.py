import gzip

import numpy as np
import pytest

from miles.core.fragments import group_fragment, ring_fragment
from miles.files import session
from miles.files.common import JobInfo
from miles.files.formats import parse_text, read_file, sniff_format, write_file
from miles.files.openbabel import input_format, read_text


@pytest.mark.parametrize("ext", ["xyz", "sdf", "mol2", "pdb", "cml"])
def test_round_trip(tmp_path, ext):
    m = group_fragment("ester").mol
    path = tmp_path / f"m.{ext}"
    write_file(str(path), m, JobInfo(title="methyl formate"))
    m2, _ = read_file(str(path))
    assert sorted(m2.symbols) == sorted(m.symbols)
    assert m2.formula() == m.formula()
    assert len(m2.bonds) == len(m.bonds)


def test_multi_frame_xyz(tmp_path):
    m = ring_fragment("benzene").mol
    text = ""
    for k in range(4):
        c = m.coords + 0.1 * k
        rows = "".join(f"{a.symbol} {x:.4f} {y:.4f} {z:.4f}\n"
                       for a, (x, y, z) in zip(m.atoms, c, strict=True))
        text += f"{len(m)}\nframe {k}\n{rows}"
    path = tmp_path / "traj.xyz"
    path.write_text(text)
    _, info = read_file(str(path))
    assert info.frames.shape == (4, len(m), 3)
    assert np.allclose(info.frames[2], m.coords + 0.2, atol=1e-3)


PDB = """\
ATOM      1  N   ALA A   1      -0.966   0.493   1.500  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.257   0.418   0.692  1.00  0.00           C
ATOM      3  C   ALA A   1      -0.094   0.017  -0.716  1.00  0.00           C
ATOM      4  O   ALA A   1      -1.056  -0.682  -0.923  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.204  -0.620   1.296  1.00  0.00           C
HETATM    6  O   HOH A 101       5.000   5.000   5.000  1.00  0.00           O
END
"""


def test_pdb_keeps_residues():
    m = read_text(PDB, "pdb")
    ca = m.atoms[1]
    assert (ca.name, ca.resname, ca.resseq, ca.chain) == ("CA", "ALA", 1, "A")
    assert m.atoms[5].hetero and m.atoms[5].resname == "HOH"


@pytest.mark.parametrize("text,fmt", [
    ("3\nwater\nO 0 0 0\nH 0.76 0.59 0\nH -0.76 0.59 0\n", "xyz"),
    ("O 0.0 0.0 0.0\nH 0.76 0.59 0.0\nH -0.76 0.59 0.0", "lines"),
    (PDB, "pdb"),
    ("CC(=O)O", "smiles"),
])
def test_sniff_pasted_text(text, fmt):
    assert sniff_format(text) == fmt


@pytest.mark.parametrize("text", [
    "3\nwater\nO 0 0 0\nH 0.76 0.59 0\nH -0.76 0.59 0\n",
    "O 0.0 0.0 0.0\nH 0.76 0.59 0.0\nH -0.76 0.59 0.0",
    "O, 0, 0.0, 0.0, 0.0\nH, 0, 0.76, 0.59, 0.0\nH, 0, -0.76, 0.59, 0.0",
    "O",
])
def test_parse_pasted_text(text):
    m = parse_text(text)
    assert m.formula() == "H2O"
    assert len(m.bonds) == 2


def test_session_round_trip(tmp_path):
    m = ring_fragment("benzene").mol
    m.measurements = [(0, 1), (0, 1, 2)]
    frames = np.stack([m.coords, m.coords + 1.0])
    path = tmp_path / "s.miles"
    display = {"rep": ["spacefill"] + [""] * (len(m) - 1), "hidden": [False, True] + [False] * (len(m) - 2)}
    info = JobInfo(title="bz", frames=frames, frame_labels=["a", "b"], view={"scale": 40.0}, display=display)
    write_file(str(path), m, info)
    with gzip.open(path, "rt") as fh:
        assert session.FORMAT in fh.read(200)
    m2, info = read_file(str(path))
    assert info.display == display
    assert m2.measurements == [(0, 1), (0, 1, 2)]
    assert info.frames.shape == (2, 12, 3) and info.frame_labels == ["a", "b"]
    assert info.view == {"scale": 40.0}
    assert m2.bonds == m.bonds


def test_sessions_with_styles_on_the_atoms_still_load(tmp_path):
    """Sessions written before the display moved out of the molecule keep rep and hidden with the atoms."""
    m = ring_fragment("benzene").mol
    data = session.to_dict(m, JobInfo(title="old"))
    del data["display"]
    data["atoms"]["rep"] = ["tube"] + [""] * (len(m) - 1)
    data["atoms"]["hidden"] = [True] + [False] * (len(m) - 1)
    _, info = session.from_dict(data)
    assert info.display["rep"][0] == "tube" and info.display["hidden"][0]


def test_qm_outputs_are_recognised(tmp_path):
    for phrase, code in (("Entering Gaussian System", "g09"), ("* O   R   C   A *", "orca")):
        path = tmp_path / f"{code}.out"
        path.write_text(f"header\n {phrase}\n")
        assert input_format(str(path)) == code
