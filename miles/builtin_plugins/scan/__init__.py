"""Scan Coordinate: a relaxed scan of a distance, angle or dihedral, with the energy of each point.

The selected atoms (2, 3 or 4, in order) define the coordinate; the points
become frames of the timeline, which plots their energies. Adds
Tools > Scan Coordinate...
"""

from miles.builtin_plugins.scan.dialog import ScanDialog, value_text
from miles.builtin_plugins.scan.relaxed_scan import ScanResult, scan
from miles.plugin import App, PluginContext, Severity, plural

NAME = "Scan Coordinate"
VERSION = "1.0"
DESCRIPTION = "Scan a distance, angle or dihedral and plot the energy profile"

ATOM_LIMIT = 400   # a force-field scan of larger structures takes too long


def scan_coordinate(app: App) -> None:
    """Scan the coordinate of the selected atoms and show the points as frames."""
    doc = app.current
    if doc is None:
        return
    atoms = doc.selection
    if not 2 <= len(atoms) <= 4:
        app.message("Select 2, 3 or 4 atoms, in order, to scan a distance, angle or dihedral",
                    severity=Severity.WARNING)
        return
    if len(doc.molecule) > ATOM_LIMIT:
        app.message(f"Scans are for structures of up to {ATOM_LIMIT} atoms", severity=Severity.WARNING)
        return
    dialog = ScanDialog(app, doc.molecule, atoms)
    if not dialog.run():
        return
    revision = doc.revision

    def done(result: ScanResult) -> None:
        if doc not in app.documents:
            return
        if doc.revision != revision or result.frames.shape[1:] != doc.molecule.coords.shape:
            app.message("The structure changed during the scan; scan it again", severity=Severity.WARNING)
            return
        doc.set_frames(result.frames, [value_text(result.kind, v) for v in result.values],
                       energies=[float(e) for e in result.energies])
        app.message(f"Scanned {plural(len(result.frames), 'point')} with {result.force_field}; "
                    "the timeline shows the energy of each", severity=Severity.WARNING)

    points = dialog.steps.value() + 1
    app.run_in_background(scan, doc.molecule.copy(), dialog.atoms, dialog.start.value(), dialog.stop.value(),
                          dialog.steps.value(), message=f"Scanning {plural(points, 'point')}", done=done)


def register(context: PluginContext) -> None:
    """Add the command after Tools > Generate Conformers."""
    context.add_command("Tools > Scan Coordinate...", scan_coordinate,
                        before=("Vibrations", "protein tools"))
