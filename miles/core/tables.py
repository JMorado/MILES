"""Tables of the data behind the plots, for the clipboard and for CSV files.

Each plot of the interface (the IR spectrum, the values along a
trajectory, the Ramachandran plot) can hand its data to a spreadsheet. The
functions here build those tables from the same data the plots draw, so
they can be used and tested without a window. Headers are plain text with
the unit in parentheses; empty cells stand for values that do not exist
(an angle at a chain end, a frame without a label).
"""

import csv
import io
import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from miles.core.trajectory import measurement_series, rmsd_atoms, rmsd_series

UNITS = {2: "angstrom", 3: "deg", 4: "deg"}
"""Unit of a measurement by its number of atoms."""


@dataclass
class Table:
    """Rows of text cells under a header.

    Attributes
    ----------
    header : list of str
        Column titles, with units in parentheses.
    rows : list of list of str
        Cells, already formatted; empty strings for missing values.
    """

    header: list[str]
    rows: list[list[str]] = field(default_factory=list)

    def tsv(self) -> str:
        """Tab-separated text with the header first, as spreadsheets paste it."""
        return "".join("\t".join(cells) + "\n" for cells in [self.header, *self.rows])

    def csv(self) -> str:
        """Comma-separated text with the header first (cells quoted where needed)."""
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(self.header)
        writer.writerows(self.rows)
        return out.getvalue()

    def write_csv(self, path: str) -> None:
        """Write :meth:`csv` to a file.

        Parameters
        ----------
        path : str
            File to write.
        """
        with open(path, "w", newline="", encoding="utf-8") as f:
            f.write(self.csv())


def format_number(value: float | None, decimals: int) -> str:
    """Format a number for a cell; None and NaN become an empty cell.

    Parameters
    ----------
    value : float or None
        The value.
    decimals : int
        Digits after the decimal point.

    Returns
    -------
    str
        The formatted value.
    """
    if value is None or math.isnan(value):
        return ""
    return f"{value:.{decimals}f}"




def frame_table(frames: ArrayLike, symbols: Sequence[str], energies: Sequence[float] | None = None,
                labels: Sequence[str] | None = None,
                measurements: Sequence[tuple[str, Sequence[int]]] = ()) -> Table:
    """Everything the timeline can plot, for every frame of a trajectory.

    Parameters
    ----------
    frames : array_like
        Coordinates, shape ``(n_frames, n_atoms, 3)``.
    symbols : sequence of str
        Element symbols, which choose the atoms of the RMSD.
    energies : sequence of float, optional
        Energy of each frame in kcal/mol.
    labels : sequence of str, optional
        Text of each frame, such as a scan coordinate or conformer energy.
    measurements : sequence of (str, sequence of int)
        Name and atoms (2, 3 or 4) of each pinned measurement.

    Returns
    -------
    Table
        One row per frame: its number from 1, the label (when any frame has
        one), the energy and the energy above the lowest frame, each
        measurement (angstrom or degrees, dihedrals between -180 and 180)
        and the RMSD from the first frame in angstrom.
    """
    xyz = np.asarray(frames, dtype=float)
    count = len(xyz)
    header = ["Frame"]
    columns: list[list[str]] = [[str(k) for k in range(1, count + 1)]]
    if labels is not None and len(labels) == count and any(labels):
        header.append("Label")
        columns.append([str(text) for text in labels])
    if energies is not None and len(energies) == count:
        lowest = min(energies)
        header += ["Energy (kcal/mol)", "Relative energy (kcal/mol)"]
        columns.append([format_number(e, 4) for e in energies])
        columns.append([format_number(e - lowest, 4) for e in energies])
    for name, atoms in measurements:
        if len(atoms) not in UNITS:
            continue
        header.append(f"{name} ({UNITS[len(atoms)]})")
        decimals = 4 if len(atoms) == 2 else 2
        columns.append([format_number(float(v), decimals) for v in measurement_series(xyz, atoms)])
    header.append("RMSD from frame 1 (angstrom)")
    columns.append([format_number(float(v), 4) for v in rmsd_series(xyz, rmsd_atoms(symbols))])
    return Table(header, [list(row) for row in zip(*columns, strict=True)])
