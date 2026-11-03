"""Align Frames: superpose every frame of a trajectory onto the frame shown.

Removes overall drift and rotation, fitting the selected atoms (three or
more) or else the heavy atoms. Adds Tools > Align Frames.
"""

import miles
from miles.plugin import App, PluginContext, Severity, plural

NAME = "Align Frames"
VERSION = "1.0"
DESCRIPTION = "Superpose the frames of a trajectory onto the frame shown"


def align_frames(app: App) -> None:
    """Superpose every frame onto the frame shown, as one undoable step."""
    doc = app.current
    if doc is None or not len(doc.molecule):
        return
    frames = doc.frames
    if frames is None or doc.frame_count < 2 or frames.shape[1] != len(doc.molecule):
        app.message("Align Frames needs a structure with several frames, such as a trajectory",
                    severity=Severity.WARNING)
        return
    selected = sorted(doc.selection)
    atoms = selected if len(selected) >= 3 else miles.rmsd_atoms(doc.molecule.symbols)
    frame = doc.frame
    doc.replace_frame_coordinates(miles.align_frames(frames, atoms, reference=frame))
    what = "the selected atoms" if atoms is selected else plural(len(atoms), "atom")
    app.message(f"Aligned {plural(doc.frame_count, 'frame')} onto frame {frame + 1}, fitting {what}")


def register(context: PluginContext) -> None:
    """Add the command after Tools > Scan Coordinate..."""
    context.add_command("Tools > Align Frames", align_frames,
                        before=("Vibrations", "protein tools"))
