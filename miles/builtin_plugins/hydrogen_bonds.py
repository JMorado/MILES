"""Hydrogen Bonds: dashed lines between donors and acceptors, counted in the status line.

Adds View > Hydrogen Bonds, a toggle of each tab. Lines start at the
hydrogen when it is drawn, otherwise at the donor, and end at the acceptor;
an end at a backbone atom shown as a cartoon moves to the ribbon. The
lines follow trajectory frames and are in exported pictures and movies.
"""

import time
from dataclasses import dataclass, field

import numpy as np
from PySide6.QtCore import QTimer

from miles import Representation
from miles.builtin_plugins.shared.hbonds import HBOND_ANGLE, HBOND_MAX, HydrogenBond, HydrogenBondFinder
from miles.plugin import Annotations, DocumentView, DrawContext, Line, PluginContext, plural
from miles.plugin.render import RIBBON_ATOMS

NAME = "Hydrogen Bonds"
VERSION = "1.0"
DESCRIPTION = "Draw hydrogen bonds as dashed lines and count them"

COLOR = "#E63A9BDC"   # blue, 90 % opaque
BUDGET = 0.008        # seconds; slower searches are not repeated for every frame of a playback
STALE_MS = 120        # how long such a search's result is reused while frames change


@dataclass
class _Tab:
    """Hydrogen bonds of one document, searched again only when the geometry changed."""

    finder: tuple[int, HydrogenBondFinder] | None = None   # (layout revision, finder)
    bonds: list[HydrogenBond] = field(default_factory=list)
    key: tuple[int, int] | None = None   # (layout revision, revision) of the bonds
    cost: float = 0.0
    found_at: float = 0.0
    timer: QTimer | None = None


def register(context: PluginContext) -> None:
    """Add the toggle, the lines and the count."""
    app = context.app

    def tab(doc: DocumentView) -> _Tab:
        cache = context.document_cache(doc)
        if "tab" not in cache:
            cache["tab"] = _Tab()
        state: _Tab = cache["tab"]
        return state

    def search(doc: DocumentView, fresh: bool = False) -> list[HydrogenBond]:
        """Hydrogen bonds among the atoms not hidden, for the geometry shown.

        What depends on the bonds is worked out once per structure. For
        large structures whose search takes longer than :data:`BUDGET`,
        frames shown in quick succession reuse the last result for up to
        :data:`STALE_MS` milliseconds, and a redraw follows to catch up;
        ``fresh`` always searches (exported pictures and movies).
        """
        state = tab(doc)
        mol = doc.molecule
        layout = doc.layout_revision   # raised when atoms, bonds or the atoms shown change
        if state.finder is None or state.finder[0] != layout:
            atoms = [i for i, rep in enumerate(doc.atom_reps()) if rep != Representation.HIDDEN]
            state.finder = (layout, HydrogenBondFinder(mol, atoms))
            state.key = None
        key = (layout, doc.revision)
        if key == state.key:
            return state.bonds
        now = time.monotonic()
        if (not fresh and state.key is not None and state.key[0] == layout and state.cost > BUDGET
                and now - state.found_at < STALE_MS / 1000):
            if state.timer is None:
                state.timer = QTimer()
                state.timer.setSingleShot(True)
                state.timer.setInterval(STALE_MS)
                state.timer.timeout.connect(
                    lambda d=doc: lines.update(d) if d in app.documents else None)
            state.timer.start()
            return state.bonds
        state.bonds = state.finder[1].find(mol.coords, max_distance.value, min_angle.value)
        state.key = key
        state.found_at = time.monotonic()
        state.cost = state.found_at - now
        return state.bonds

    def ends(doc: DocumentView, bonds: list[HydrogenBond]) -> np.ndarray:
        """Atoms at the two ends of each line: the hydrogen if drawn, else the donor, and the acceptor."""
        drawn = doc.drawn_atoms()
        pairs = [(b.hydrogen if b.hydrogen is not None and drawn[b.hydrogen] else b.donor, b.acceptor)
                 for b in bonds]
        return np.array(pairs, dtype=np.int64).reshape(-1, 2)

    def draw(doc: DocumentView, drawing: DrawContext) -> Annotations | None:
        if not toggle.value(doc) or not len(doc.molecule):
            return None
        pairs = ends(doc, search(doc, fresh=drawing.export))
        return Annotations(lines=[Line(int(a), int(b), COLOR, 1.6, (3.0, 2.5), round_caps=False, ribbon=True)
                                  for a, b in pairs])

    def count(doc: DocumentView) -> str | None:
        if not toggle.value(doc) or not len(doc.molecule):
            return None
        bonds = search(doc, fresh=True)
        total = len(bonds)
        text = plural(total, "hydrogen bond")
        mol = doc.molecule
        drawn = doc.drawn_atoms()
        inside = sum(1 for pair in ends(doc, bonds)   # side chains or bases inside a cartoon
                     if any(not drawn[i] and mol.atoms[i].name not in RIBBON_ATOMS for i in pair))
        return f"{text} ({total - inside} drawn)" if inside else text

    def switched(doc: DocumentView, on: bool) -> None:
        lines.update(doc)
        status.update(doc)

    def criteria_changed(value: float) -> None:
        for doc in app.documents:
            tab(doc).key = None   # search again with the new criteria
        lines.update()
        status.update()

    max_distance = context.add_setting("max_distance", "Largest donor-acceptor distance", HBOND_MAX,
                                       minimum=2.5, maximum=4.5, step=0.1, unit="angstrom",
                                       changed=criteria_changed)
    min_angle = context.add_setting("min_angle", "Smallest D-H...A angle", HBOND_ANGLE, minimum=90.0,
                                    maximum=180.0, step=5.0, unit="degrees", changed=criteria_changed)
    lines = context.add_annotations(draw)
    status = context.add_status(count)
    toggle = context.add_toggle("View > Hydrogen Bonds", "hydrogen_bonds", switched, shortcut="Ctrl+Shift+B",
                                after="Stereo Labels")
