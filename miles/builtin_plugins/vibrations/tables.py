"""Tables of vibrational modes and IR spectra, for copying and exporting."""

from collections.abc import Sequence

import numpy as np

from miles import VibrationalMode
from miles.builtin_plugins.vibrations.modes import ir_spectrum
from miles.plugin.widgets import Table, format_number


def mode_table(modes: Sequence[VibrationalMode]) -> Table:
    """Wavenumber and infrared intensity of every vibrational mode.

    Parameters
    ----------
    modes : sequence of VibrationalMode
        Normal modes, in their order in the file.

    Returns
    -------
    Table
        One row per mode: its number from 1, the wavenumber in cm-1
        (negative for imaginary modes, as most programs print them) and the
        intensity in km/mol (empty when the file gives none).
    """
    rows = [[str(k), format_number(m.frequency, 2), format_number(m.intensity, 4)]
            for k, m in enumerate(modes, start=1)]
    return Table(["Mode", "Wavenumber (cm-1)", "IR intensity (km/mol)"], rows)


def spectrum_table(modes: Sequence[VibrationalMode], low: float = 400.0, high: float = 4000.0,
                   step: float = 1.0, width: float = 20.0) -> Table:
    """Tabulate the broadened infrared spectrum on a regular grid of wavenumbers.

    Parameters
    ----------
    modes : sequence of VibrationalMode
        Normal modes.
    low, high : float, optional
        Range of wavenumbers in cm-1.
    step : float, optional
        Spacing of the grid in cm-1.
    width : float, optional
        Full width at half maximum of each band in cm-1, as in the plot.

    Returns
    -------
    Table
        Wavenumber and absorbance (scaled so the highest point is 1), from
        high to low wavenumber as IR spectra are usually listed.
    """
    grid = np.arange(high, low - step / 2, -step)
    curve = ir_spectrum(modes, grid, width)
    rows = [[format_number(float(w), 1), format_number(float(a), 6)]
            for w, a in zip(grid, curve, strict=True)]
    return Table(["Wavenumber (cm-1)", "Absorbance (relative)"], rows)
