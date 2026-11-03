"""Panel for volumetric data: which grid or orbital, the isovalue and how it is drawn."""

import math
from typing import Any, override

from PySide6.QtCore import QRectF, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QLinearGradient,
    QPainter,
    QPaintEvent,
)
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QVBoxLayout, QWidget

from miles import HARTREE_EV, OrbitalSet, Volume, mapped_isosurface
from miles.builtin_plugins.grids.combine import Combination, combine_volumes
from miles.builtin_plugins.grids.combine_dialog import CombineDialog
from miles.builtin_plugins.grids.levels import OrbitalLevels
from miles.plugin import DocumentPanel, DocumentView, Theme
from miles.plugin.render import MAPPED_COLORS, VolumeRequest, symmetric_range
from miles.plugin.widgets import LIGHT, Dial, ScrubField, button, label, theme_palette

PLUS_MINUS = "\u00b1"
UPDATE_DELAY_MS = 40   # recompute the isosurface at most this often while dragging
TICK_SPAN = 25         # the dial's tick spacing suits isovalues up to this many times the current one


Entry = tuple[str, int]
"""A row of the grid menu: ("volume", index), ("orbital", index), ("density", 0) or ("spin", 0)."""

UNSAVED = object()
"""What :meth:`VolumePanel.session_choice` returns when a document's session said nothing about grids."""


def request_from_saved(volumes: list[Volume], saved: Any) -> VolumeRequest | None:
    """Rebuild the isosurface a session showed, if it still fits the document's grids.

    Parameters
    ----------
    volumes : list of Volume
        The document's grids.
    saved : object
        What :func:`request_to_saved` gave, None for no isosurface.

    Returns
    -------
    VolumeRequest or None
        The request, or None when nothing was shown or it no longer fits.
    """
    if not isinstance(saved, dict):
        return None
    try:
        index = int(saved["index"])
        mapped = saved.get("mapped")
        mapped = None if mapped is None else int(mapped)
        request = VolumeRequest(index, float(saved["level"]), bool(saved.get("see_through", True)),
                                mapped, float(saved.get("mapped_range", 0.0)))
    except (KeyError, TypeError, ValueError):
        return None
    if not 0 <= request.index < len(volumes):
        return None
    if request.mapped is not None and not 0 <= request.mapped < len(volumes):
        return VolumeRequest(request.index, request.level, request.see_through)
    return request


def request_to_saved(request: VolumeRequest | None) -> dict[str, Any] | None:
    """Describe an isosurface request as JSON-ready data for a session (None for none)."""
    if request is None:
        return None
    return {"index": request.index, "level": request.level, "see_through": request.see_through,
            "mapped": request.mapped, "mapped_range": request.mapped_range}


def _evaluate(orbitals: OrbitalSet, entry: Entry) -> Volume:
    """Worker-thread job: the grid of an orbital or a density."""
    kind, index = entry
    if kind == "density":
        return orbitals.density()
    if kind == "spin":
        return orbitals.spin_density()
    return orbitals.evaluate(index)


def _ticks(top: float) -> tuple[float, int]:
    """Minor tick spacing and decimals for values from 0 to ``top`` (20 to 200 ticks)."""
    minor = 10 ** math.floor(math.log10(top)) / 20
    return minor, max(2, -math.floor(math.log10(minor)) + 1)


class ColorBar(QWidget):
    """Legend of the colours of a mapped grid, from negative on the left to positive."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedSize(72, 10)
        self.setToolTip("Red negative, white zero, blue positive")

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        gradient = QLinearGradient(0, 0, self.width(), 0)
        for stop, rgb in zip((0.0, 0.5, 1.0), MAPPED_COLORS, strict=True):
            gradient.setColorAt(stop, QColor.fromRgbF(*(float(c) for c in rgb)))
        p.setPen(QColor(128, 128, 128, 90))
        p.setBrush(gradient)
        p.drawRoundedRect(QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5), 3, 3)
        p.end()


class VolumePanel(DocumentPanel):
    """Shown when the document has volumetric data or molecular orbitals.

    Grids read from a file are listed by name. Orbitals are listed with
    their HOMO/LUMO label and energy; choosing one evaluates it on a grid
    in the background and adds that grid to the document, so it is computed
    only once and saved with sessions; an energy diagram of the orbitals
    (:class:`OrbitalLevels`) offers them by clicking. A density-like surface
    can be coloured by a second grid, such as the electrostatic potential.

    The panel keeps the isosurface each document shows (:meth:`request`) in
    the document's store and emits ``volumeChanged(document)`` when it
    changes; whoever draws the isosurfaces follows that.
    :meth:`session_choice` gives what a session said a document showed. The
    panel opens by itself for a structure with grids or orbitals.
    """

    addGridRequested = Signal()
    saveGridRequested = Signal()
    volumeChanged = Signal(object)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.wanted = True   # it opens by itself whenever the tab in front has grids or orbitals
        self.theme: Theme = LIGHT
        self._entries: list[Entry] = []                # one per row of the grid menu
        self._mapped_entries: list[int | None] = []    # volume index per "Color by" row
        self._wanted: tuple[DocumentView, Entry] | None = None   # orbital or density computed for display
        self._pending = QTimer(self)
        self._pending.setSingleShot(True)
        self._pending.setInterval(UPDATE_DELAY_MS)
        self._pending.timeout.connect(self._update)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 6, 8, 4)
        outer.setSpacing(2)
        lay = QHBoxLayout()
        lay.setSpacing(10)
        outer.addLayout(lay)
        self.choice = QComboBox()
        self.choice.setMaxVisibleItems(16)
        self.choice.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.choice.currentIndexChanged.connect(self._chosen)
        lay.addWidget(self.chevron)
        lay.addWidget(self.choice)
        lay.addWidget(label("Isovalue", "Muted"))
        self.sign = label(PLUS_MINUS, "Muted")
        lay.addWidget(self.sign)
        self.level = ScrubField()
        self.level.edited.connect(self._level_edited)
        lay.addWidget(self.level)
        self.see_through = button("Seg", "See-Through", checkable=True)
        self.see_through.clicked.connect(self._update)
        lay.addWidget(self.see_through)
        self.visible = button("Seg", "Show", checkable=True)
        self.visible.setChecked(True)
        self.visible.clicked.connect(self._update)
        lay.addWidget(self.visible)
        self.dial = Dial()   # right under the isovalue it sets
        self.dial.changed.connect(self._dialed)
        self.body_layout.setSpacing(2)
        self.body_layout.addWidget(self.dial)
        self.levels = OrbitalLevels()
        self.levels.orbitalClicked.connect(self.choose_orbital)
        self.levels.hide()
        self.body_layout.addWidget(self.levels)

        self.mapping_row = QWidget()
        row = QHBoxLayout(self.mapping_row)
        row.setContentsMargins(0, 2, 0, 2)
        row.setSpacing(10)
        self.mapping_label = label("Color by", "Muted")
        row.addWidget(self.mapping_label)
        self.mapped = QComboBox()
        self.mapped.currentIndexChanged.connect(self._mapped_chosen)
        row.addWidget(self.mapped)
        self.range_sign = label(PLUS_MINUS, "Muted")
        row.addWidget(self.range_sign)
        self.range = ScrubField()
        self.range.setToolTip("Values at minus and plus this get the full colours")
        self.range.edited.connect(lambda _: self._pending.start())
        row.addWidget(self.range)
        self.bar = ColorBar()
        row.addWidget(self.bar)
        row.addStretch(1)
        self.add_grid = button("Link", "Add Grid...", "Read grids for this structure from another file")
        self.add_grid.clicked.connect(lambda: self.addGridRequested.emit())
        row.addWidget(self.add_grid)
        self.combine_link = button("Link", "Combine...", "Make a grid as the difference or sum of two grids")
        self.combine_link.clicked.connect(self.combine_dialog)
        row.addWidget(self.combine_link)
        self.save_grid = button("Link", "Save Grid...", "Write the grid shown to a cube file")
        self.save_grid.clicked.connect(lambda: self.saveGridRequested.emit())
        row.addWidget(self.save_grid)
        self.body_layout.addWidget(self.mapping_row)

        self.status = label("", "Muted")
        self.status.hide()
        self.body_layout.addWidget(self.status)
        outer.addWidget(self.body)

    def set_theme(self, theme: Theme) -> None:
        """Recolour the dial and the orbital diagram."""
        self.theme = theme
        self.set_icon_color(theme.muted)
        self.dial.theme = theme
        self.dial.update()
        self.levels.theme = theme
        self.levels.update()
        for combo in (self.choice, self.mapped):   # readable when the app and system themes differ
            theme_palette(combo, theme)

    def session_choice(self, doc: DocumentView) -> Any:
        """Return what a document's session said it showed, or :data:`UNSAVED`.

        Sessions keep it in the plugin's state; sessions written before the
        grids were a plugin kept it with the document's view.
        """
        state = self.saved(doc)
        return state["volume"] if "volume" in state else doc.info.view.get("volume", UNSAVED)

    def request(self, doc: DocumentView | None = None) -> VolumeRequest | None:
        """Return the isosurface a document shows (default: the current one's), or None.

        A document seen for the first time shows what its session saved,
        else its first grid.
        """
        doc = doc if doc is not None else self.document
        if doc is None:
            return None
        store = self.store(doc)
        if "request" not in store:
            saved = self.session_choice(doc)
            volumes = doc.info.volumes
            if saved is not UNSAVED:
                store["request"] = request_from_saved(volumes, saved)
            else:
                store["request"] = VolumeRequest(0, volumes[0].suggested_level()) if volumes else None
        request: VolumeRequest | None = store["request"]
        return request

    def volume_shown(self, doc: DocumentView | None = None) -> tuple[Volume, VolumeRequest] | None:
        """Return the grid a document shows as an isosurface and how, if any."""
        doc = doc if doc is not None else self.document
        request = self.request(doc)
        if doc is None or request is None or not 0 <= request.index < len(doc.info.volumes):
            return None
        return doc.info.volumes[request.index], request

    def mapped_volume(self, doc: DocumentView | None = None) -> Volume | None:
        """Return the grid colouring a document's isosurface, if any."""
        doc = doc if doc is not None else self.document
        shown = self.volume_shown(doc)
        if doc is None or shown is None:
            return None
        mapped = shown[1].mapped
        if mapped is None or not 0 <= mapped < len(doc.info.volumes):
            return None
        return doc.info.volumes[mapped]

    def set_volume(self, request: VolumeRequest | None, doc: DocumentView | None = None) -> None:
        """Show an isosurface of a document's grids (default: the current document), or none."""
        doc = doc if doc is not None else self.document
        if doc is None:
            return
        self.store(doc)["request"] = request
        self.volumeChanged.emit(doc)

    @override
    def has_content(self, document: DocumentView) -> bool:
        """Whether the document has grids or orbitals."""
        return bool(document.info.volumes or document.info.orbitals)

    @override
    def show_contents(self, document: DocumentView) -> None:
        """Show the controls for the document's volumes and orbitals."""
        doc = document
        orbitals = doc.info.orbitals
        restored = self.session_choice(doc) is not UNSAVED   # a session chose what to show, possibly nothing
        store = self.store(doc)
        if orbitals is not None and self.request(doc) is None and not store.get("started") and not restored:
            store["started"] = True
            self._compute(doc, ("orbital", orbitals.homo()))   # show the HOMO when the file is first opened
        self._fill_choice()
        shown = self.volume_shown(doc)
        row = self.choice.currentIndex()
        index = self._volume_index(row) if row >= 0 else None
        volume = doc.info.volumes[index] if index is not None else None
        top = float(abs(volume.values).max()) if volume is not None else 0.0
        if shown is not None:
            level = shown[1].level
        elif volume is not None:
            level = volume.suggested_level()
        else:
            level = self.level.value or 0.02
        # ticks follow the isovalue rather than the peak, which for densities is far larger
        minor, decimals = _ticks(min(top, TICK_SPAN * level) or top or 1.0)
        self.level.configure(0.0, top or 1.0, minor / 10, decimals)
        self.level.set_value(level)
        self.dial.configure(0.0, top or 1.0, minor / 12, minor, minor * 10, max(0, decimals - 1))
        self.dial.set_value(level)
        self.sign.setVisible(volume.signed if volume is not None else orbitals is not None)
        request = self.request(doc)
        self.visible.setChecked(request is not None or self._wanted is not None)
        self.see_through.setChecked(request.see_through if request is not None
                                    else VolumeRequest(0, 0.0).see_through)
        self._fill_mapping(index, volume)
        self.levels.set_orbitals(orbitals)
        self.levels.setVisible(orbitals is not None)
        kind, number = self._entries[row] if 0 <= row < len(self._entries) else ("", -1)
        self.levels.set_current(number if kind == "orbital" else None)

    def choose_orbital(self, index: int) -> None:
        """Show orbital ``index``, computing it first if needed (as choosing it in the menu does)."""
        if ("orbital", index) not in self._entries:
            return
        self.visible.setChecked(True)
        row = self._entries.index(("orbital", index))
        if row == self.choice.currentIndex():
            self._update()
            self.refresh()
        else:
            self.choice.setCurrentIndex(row)   # emits, which shows or computes it

    def _fill_choice(self) -> None:
        """List the document's grids and orbitals and select the one on screen."""
        assert self.document is not None
        doc = self.document
        volumes, orbitals = doc.info.volumes, doc.info.orbitals
        computed = self._derived()
        self._entries = []
        texts = []
        for k, volume in enumerate(volumes):
            if volume.name not in computed:
                self._entries.append(("volume", k))
                texts.append(volume.name or f"Grid {k + 1}")
        if orbitals is not None:
            self._entries.append(("density", 0))
            texts.append(OrbitalSet.DENSITY_NAME)
            if orbitals.unrestricted:
                self._entries.append(("spin", 0))
                texts.append(OrbitalSet.SPIN_NAME)
            for k, orbital in enumerate(orbitals.orbitals):
                self._entries.append(("orbital", k))
                tag = orbitals.label(k)
                texts.append(f"MO {k + 1}{'  ' + tag if tag else ''}   {orbital.energy * HARTREE_EV:.2f} eV")
        row = 0
        request = self.request(doc)
        if self._wanted is not None and self._wanted[0] is doc:
            row = self._entries.index(self._wanted[1])
        elif request is not None and 0 <= request.index < len(volumes):
            name = volumes[request.index].name
            if orbitals is not None and name in computed:
                row = self._entries.index(computed[name])
            else:
                row = self._entries.index(("volume", request.index))
        elif orbitals is not None:
            row = self._entries.index(("orbital", orbitals.homo()))
        self.choice.blockSignals(True)
        self.choice.clear()
        self.choice.addItems(texts)
        self.choice.setCurrentIndex(row)
        self.choice.setVisible(len(texts) > 1)
        self.choice.blockSignals(False)

    def _fill_mapping(self, index: int | None, volume: Volume | None) -> None:
        """Fill the bottom row: "Color by" for unsigned data such as densities, and the grid links."""
        assert self.document is not None
        request = self.request()
        offered = volume is not None and not volume.signed
        self.mapping_row.setVisible(volume is not None)
        self.add_grid.setVisible(offered)
        self.combine_link.setVisible(len(self.document.info.volumes) >= 2)
        if not offered:
            for widget in (self.mapping_label, self.mapped, self.range_sign, self.range, self.bar):
                widget.setVisible(False)
            return
        computed = self._orbital_names()
        volumes = self.document.info.volumes
        others = [k for k, v in enumerate(volumes) if k != index and v.name not in computed]
        self._mapped_entries = [None, *others]
        mapped = request.mapped if request is not None else None
        self.mapped.blockSignals(True)
        self.mapped.clear()
        self.mapped.addItems(["None"] + [volumes[k].name or f"Grid {k + 1}" for k in others])
        self.mapped.setCurrentIndex(self._mapped_entries.index(mapped) if mapped in others else 0)
        self.mapped.blockSignals(False)
        self.mapping_label.setVisible(bool(others))
        self.mapped.setVisible(bool(others))
        active = mapped in others
        self.range_sign.setVisible(active)
        self.range.setVisible(active)
        self.bar.setVisible(active)
        if active and request is not None and request.mapped_range > 0:
            self._configure_range(request.mapped_range)

    def selected_volume(self) -> Volume | None:
        """Return the grid chosen in the panel, or None for an orbital not computed yet."""
        document = self.document
        row = self.choice.currentIndex()
        if document is None or not 0 <= row < len(self._entries):
            return None
        index = self._volume_index(row)
        return document.info.volumes[index] if index is not None else None

    def _orbital_names(self) -> dict[str, int]:
        """Grid names of the document's orbitals, mapped to the orbital index."""
        assert self.document is not None
        orbitals = self.document.info.orbitals
        if orbitals is None:
            return {}
        return {orbitals.name(k): k for k in range(len(orbitals.orbitals))}

    def _derived(self) -> dict[str, Entry]:
        """Names of the grids computed from the document's orbitals, mapped to their menu entries."""
        assert self.document is not None
        orbitals = self.document.info.orbitals
        if orbitals is None:
            return {}
        names: dict[str, Entry] = {name: ("orbital", k) for name, k in self._orbital_names().items()}
        names[OrbitalSet.DENSITY_NAME] = ("density", 0)
        if orbitals.unrestricted:
            names[OrbitalSet.SPIN_NAME] = ("spin", 0)
        return names

    @staticmethod
    def _entry_name(orbitals: OrbitalSet, entry: Entry) -> str:
        kind, index = entry
        names = {"density": OrbitalSet.DENSITY_NAME, "spin": OrbitalSet.SPIN_NAME}
        return names.get(kind) or orbitals.name(index)

    def _volume_index(self, row: int) -> int | None:
        """Index into the document's volumes for a row; None for an orbital not computed yet."""
        assert self.document is not None
        kind, index = self._entries[row]
        if kind == "volume":
            return index
        orbitals = self.document.info.orbitals
        assert orbitals is not None
        names = [v.name for v in self.document.info.volumes]
        name = self._entry_name(orbitals, self._entries[row])
        return names.index(name) if name in names else None

    def _showing_orbital(self) -> bool:
        """Whether the isosurface on screen is one of the document's orbitals."""
        assert self.document is not None
        shown = self.volume_shown()
        return shown is not None and shown[0].name in self._orbital_names()

    def _mapped_index(self) -> int | None:
        """Volume chosen to colour the surface, if the row is shown and one is chosen."""
        row = self.mapped.currentIndex()
        if self.mapping_row.isHidden() or self.mapped.isHidden() or not 0 <= row < len(self._mapped_entries):
            return None
        return self._mapped_entries[row]

    def _configure_range(self, value: float) -> None:
        minor, decimals = _ticks(value)
        self.range.configure(0.0, value * 20, minor / 10, decimals)
        self.range.set_value(value)

    def _chosen(self, row: int) -> None:
        document = self.document
        if document is None or not 0 <= row < len(self._entries):
            return
        index = self._volume_index(row)
        kind, _ = self._entries[row]
        if index is None:
            self._compute(document, self._entries[row])
            return
        if kind != "orbital" or not self._showing_orbital():   # orbitals keep the isovalue between them
            level = document.info.volumes[index].suggested_level()
            self.level.set_value(level)
            self.dial.set_value(level)
        self._wanted = None
        self.status.hide()
        self._update()
        self.refresh()

    def _mapped_chosen(self, row: int) -> None:
        document = self.document
        shown = self.volume_shown()
        if document is None or shown is None or not 0 <= row < len(self._mapped_entries):
            return
        mapped = self._mapped_entries[row]
        if mapped is not None:
            _, values = mapped_isosurface(shown[0], self.level.value, document.info.volumes[mapped])
            self._configure_range(symmetric_range(values))
            if shown[1].mapped is None:
                self.see_through.setChecked(False)   # colours read best on an opaque surface
        self._update()
        self.refresh()

    def map_grid(self, index: int) -> None:
        """Colour the density-like surface on screen by volume ``index`` (after Add Grid)."""
        shown = self.volume_shown()
        if shown is None or shown[0].signed or index == shown[1].index:
            self.refresh()
            return
        self.refresh()
        if index in self._mapped_entries:
            self.mapped.setCurrentIndex(self._mapped_entries.index(index))   # emits, which maps it

    def combine_dialog(self) -> None:
        """Ask for two grids and an operation, then add and show the result."""
        document = self.document
        if document is None or not document.info.volumes:
            return
        volumes = document.info.volumes
        row = self.choice.currentIndex()
        current = self._volume_index(row) if 0 <= row < len(self._entries) else None
        if self.context is None:
            return
        dialog = CombineDialog(self.context.app, [v.name or f"Grid {k + 1}" for k, v in enumerate(volumes)],
                               current or 0)
        if dialog.run():
            self.combine(dialog.first, dialog.second, dialog.combination, dialog.name)

    def combine(self, first: int, second: int, combination: Combination = Combination.DIFFERENCE,
                name: str = "") -> int:
        """Add the difference or sum of two of the document's grids and show it.

        Parameters
        ----------
        first, second : int
            Indices of grids A and B in the document's volumes.
        combination : Combination, optional
            ``A - B`` or ``A + B``.
        name : str, optional
            Name of the new grid (default from the grids' names); a number is
            added if the document already has a grid of that name.

        Returns
        -------
        int
            Index of the new grid.
        """
        assert self.document is not None
        doc = self.document
        volumes = doc.info.volumes
        result = combine_volumes(volumes[first], volumes[second], combination, name)
        names = {v.name for v in volumes}
        base, n = result.name, 2
        while result.name in names:
            result.name = f"{base} ({n})"
            n += 1
        index = doc.add_volume(result)   # one undoable step
        self.visible.setChecked(True)
        self._wanted = None
        self.set_volume(VolumeRequest(index, result.suggested_level(), self.see_through.isChecked()), doc)
        self.refresh()
        return index

    def _compute(self, doc: DocumentView, entry: Entry) -> None:
        orbitals = doc.info.orbitals
        assert orbitals is not None
        self._wanted = (doc, entry)
        kind, index = entry
        if kind == "orbital":
            tag = orbitals.label(index)
            self.status.setText(f"Computing MO {index + 1}{' (' + tag + ')' if tag else ''} ...")
        else:
            self.status.setText(f"Computing the {self._entry_name(orbitals, entry).lower()} ...")
        self.status.show()
        if self.document is not None and self.document is doc:
            self.levels.set_current(index if kind == "orbital" else None)
        if self.context is None:
            return
        self.context.app.run_in_background(
            _evaluate, orbitals, entry, message=self.status.text().removesuffix(" ..."),
            done=lambda volume: self._computed(doc, entry, volume), failed=self._failed)

    def _computed(self, doc: DocumentView, entry: Entry, volume: Volume) -> None:
        if self.context is None or doc not in self.context.app.documents:
            return
        names = [v.name for v in doc.info.volumes]
        if volume.name not in names:
            doc.keep_volume(volume)   # worked out from the document's orbitals: kept, but not an edit
            names.append(volume.name)
        if self._wanted != (doc, entry) or self.document is None or self.document is not doc:
            return   # the user moved on; keep the grid for later
        self._wanted = None
        self.status.hide()
        request = self.request(doc)
        if entry[0] == "orbital" and request is not None and self._showing_orbital():
            level = request.level
        else:
            level = volume.suggested_level()
        self.set_volume(VolumeRequest(names.index(volume.name), level, self.see_through.isChecked()), doc)
        self.refresh()

    def _failed(self, message: str) -> None:
        self._wanted = None
        self.status.setText(f"Could not compute it: {message}")
        self.status.show()
        self.adjustSize()
        self.relayout.emit()

    def _level_edited(self, value: float) -> None:
        self.dial.set_value(value)
        self._pending.start()

    def _dialed(self, value: float) -> None:
        self.level.set_value(value)
        self._pending.start()

    def _update(self) -> None:
        document = self.document
        if document is None or not self._entries:
            return
        if not self.visible.isChecked():
            self.set_volume(None)
            return
        row = max(0, self.choice.currentIndex())
        index = self._volume_index(row)
        if index is None:
            self._compute(document, self._entries[row])
            return
        mapped = self._mapped_index()
        if mapped == index or document.info.volumes[index].signed:
            mapped = None
        self.set_volume(VolumeRequest(index, self.level.value, self.see_through.isChecked(), mapped,
                                      self.range.value if mapped is not None else 0.0))
