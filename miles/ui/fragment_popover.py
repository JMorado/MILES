"""Fragment popover: periodic table, rings, groups, biomolecules, ligands and search."""

import os
from collections.abc import Callable
from pathlib import Path
from typing import override

from PySide6.QtCore import QPoint, QPointF, QSize, QStandardPaths, Qt, Signal
from PySide6.QtGui import QIcon, QMouseEvent, QPainter, QPaintEvent
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLineEdit,
    QScrollArea,
    QStackedWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from miles.core.elements import ELEMENTS, element
from miles.core.fragments import (
    Category,
    Fragment,
    all_group_fragments,
    all_ring_fragments,
    element_fragment,
    element_types,
    library_fragments,
    molecule_fragment,
)
from miles.core.molecule import Molecule
from miles.files.common import JobInfo
from miles.render.camera import Camera, nice_rotation, thumbnail_camera
from miles.render.painter import paint_scene, render_pixmap
from miles.render.picking import AtomHit, Picker
from miles.render.scene import Highlight, Style, build_scene
from miles.render.theme import Theme
from miles.ui.builder import Builder
from miles.ui.controls import Segmented, button, label, rule
from miles.ui.lookup import StructureLookup, StructureSource, choose_source, search_hint
from miles.ui.text import SEPARATOR

TILE_COLUMNS = 5
PREVIEW_FILL = 0.8
PREVIEW_TURN = 0.012   # radians per pixel

SHORT_LABELS = {
    "Trigonal planar": "Planar",
    "Bent (double bond)": "Imine",
    "Linear (triple bond)": "Nitrile",
    "Tetrahedral (cation)": "Cation",
    "Pyramidal (cation)": "Cation",
    "Double bond": "Carbonyl",
    "Pyramidal (sulfoxide)": "Sulfoxide",
    "Tetrahedral (sulfone)": "Sulfone",
    "Tetrahedral (P=O)": "P=O",
    "Trigonal bipyramidal": "TBP",
    "Square planar": "Square",
    "Tetrahedral (anion)": "Anion",
}
"""Names of bonding types short enough for the buttons under the periodic table."""

TABS: list[tuple[str, Callable[[], list[Fragment]] | None]] = [
    ("Elements", None),
    ("Rings", all_ring_fragments),
    ("Groups", all_group_fragments),
    ("Biomolecules", lambda: library_fragments(Category.BIO)),
    ("Ligands", lambda: library_fragments(Category.LIGAND)),
]
CATEGORY_TABS: dict[str, str] = {Category.BIO: "Biomolecules", Category.LIGAND: "Ligands"}
"""Tabs that fragment libraries add to when they use the bundled categories."""

FOLDER_ENV = "MILES_FRAGMENT_DIR"
"""Environment variable naming another fragment folder (used by the tests)."""


def fragment_folder() -> Path:
    """Folder whose ``.json`` fragment libraries are added at start-up (it may not exist).

    Returns
    -------
    pathlib.Path
        ``$MILES_FRAGMENT_DIR`` if set, otherwise ``fragments`` in the
        application data folder, next to the plugin folder.
    """
    override = os.environ.get(FOLDER_ENV)
    if override:
        return Path(override)
    base = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
    return Path(base or Path.home() / ".miles") / "fragments"


def _first_type(symbol: str) -> Fragment:
    """Element fragment with the element's default bonding type."""
    return element_fragment(symbol, element_types(symbol)[0])


class FragmentPreview(QWidget):
    """Rotatable preview of the current fragment; click a terminal atom to bond through it."""

    def __init__(self, builder: Builder, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.builder = builder
        self.theme: Theme | None = None
        self.camera = Camera()
        self._fragment: Fragment | None = None
        self._press: QPointF | None = None
        self._last = QPointF()
        self._moved = False
        self.setFixedSize(200, 168)
        self.setMouseTracking(True)
        self.setCursor(Qt.CursorShape.OpenHandCursor)

    def sync(self) -> None:
        """Frame a newly chosen fragment, keeping the user's rotation if there is one."""
        frag = self.builder.fragment
        if frag is not self._fragment:
            self._fragment = frag
            rotation = self.builder.preview_rotation
            if rotation is None:
                rotation = nice_rotation(frag.mol.coords, tilt=frag.category != Category.ELEMENT)
                self.builder.preview_rotation = rotation.copy()
            self.camera = thumbnail_camera(frag.mol, self.width(), self.height(), rotation=rotation.copy(),
                                           fill=PREVIEW_FILL)
        self.update()

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None:
            return
        frag = self.builder.fragment
        highlights = {frag.attach: Highlight.ATTACH} if frag.attach is not None else None
        scene = build_scene(frag.mol, Style(depth_fade=False), self.theme, highlights)
        p = QPainter(self)
        paint_scene(p, self.width(), self.height(), scene, self.camera, self.theme, Style(), background=False)
        p.end()

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        self._press = self._last = e.position()
        self._moved = False

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        if self._press is None:
            return
        if not self._moved and (e.position() - self._press).manhattanLength() < 3:
            return
        d = e.position() - self._last
        self._moved = True
        self._last = e.position()
        self.camera.turn(d.y() * PREVIEW_TURN, d.x() * PREVIEW_TURN)
        self.builder.preview_rotation = self.camera.rotation.copy()
        self.update()

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        if self._press is not None and not self._moved:
            mol = self.builder.fragment.mol
            hit = Picker(mol, Style()).pick(self.camera, self.width(), self.height(), e.position().x(),
                                            e.position().y())
            if isinstance(hit, AtomHit):
                rotation = self.builder.preview_rotation
                self.builder.set_attach(hit.index)
                self.builder.preview_rotation = rotation
        self._press = None
        self.update()


class FragmentPopover(QWidget):
    """Elements, rings, groups and search in one popover, with a live preview.

    Signals
    -------
    placeRequested()
        A fragment was built from the search box and should be placed.
    structureLoaded(mol, info)
        A PDB entry was downloaded and should open in its own tab.
    """

    placeRequested = Signal()
    structureLoaded = Signal(object, object)

    def __init__(self, builder: Builder, parent: QWidget | None = None,
                 sources: list[StructureSource] | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.builder = builder
        self.theme: Theme | None = None
        self._thumbnail_theme = ""
        self._tiles: list[tuple[QToolButton, Fragment]] = []
        self._filled: set[int] = set()
        self._tabs = list(TABS)
        self._added: dict[int, list[Fragment]] = {}   # fragments of libraries, by tab
        self.lookup = StructureLookup(sources)
        self.lookup.fragment.connect(self._found_fragment)
        self.lookup.structure.connect(self._found_structure)
        self.lookup.failed.connect(self._failed)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        frame = QFrame()
        frame.setObjectName("Popover")
        outer.addWidget(frame)
        lay = QVBoxLayout(frame)
        lay.setContentsMargins(14, 12, 14, 14)
        lay.setSpacing(10)

        self.search = QLineEdit()
        self.search.setObjectName("Search")
        self.search.setClearButtonEnabled(True)
        self.search.setPlaceholderText(search_hint(self.lookup.sources))
        self.search.returnPressed.connect(self.look_up)
        lay.addWidget(self.search)
        self.tabs = Segmented([name for name, _ in TABS])
        self.tabs.changed.connect(self._show_tab)
        top = QHBoxLayout()
        top.addWidget(self.tabs)
        top.addStretch(1)
        lay.addLayout(top)

        body = QHBoxLayout()
        body.setSpacing(14)
        self.pages = QStackedWidget()
        self.pages.addWidget(self._element_page())
        for _ in self._tabs[1:]:
            self.pages.addWidget(QWidget())
        body.addWidget(self.pages, 1)
        body.addWidget(rule(vertical=True))
        side = QVBoxLayout()
        side.setSpacing(4)
        self.preview = FragmentPreview(builder)
        side.addWidget(self.preview)
        self.name = label("", "Title")
        self.name.setWordWrap(True)
        side.addWidget(self.name)
        self.hint = label("", "Muted")
        self.hint.setWordWrap(True)
        self.hint.setFixedWidth(200)
        side.addWidget(self.hint)
        side.addStretch(1)
        body.addLayout(side)
        lay.addLayout(body)
        builder.fragmentChanged.connect(self.sync)

    # -- pages ------------------------------------------------------------------------------

    def _element_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(10)
        grid = QGridLayout()
        grid.setSpacing(2)
        self.cells: dict[str, QToolButton] = {}
        for symbol, el in ELEMENTS.items():
            b = button("Cell", symbol, f"{el.number}  {el.name}", checkable=True)
            b.setCursor(Qt.CursorShape.ArrowCursor)
            b.setFixedSize(26, 23)
            b.clicked.connect(lambda _=False, s=symbol: self._choose(_first_type(s)))
            grid.addWidget(b, el.row + (1 if el.row >= 8 else 0), el.column)   # gap above the f-block
            self.cells[symbol] = b
        gap = QWidget()
        gap.setFixedHeight(6)
        grid.addWidget(gap, 8, 1)
        lay.addLayout(grid)
        lay.addWidget(rule())
        self.variant_title = label("", "Muted")
        lay.addWidget(self.variant_title)
        self.variant_row = QHBoxLayout()
        self.variant_row.setSpacing(4)
        lay.addLayout(self.variant_row)
        lay.addStretch(1)
        return page

    def _fill_tiles(self, page: QWidget, fragments: list[Fragment]) -> None:
        outer = QVBoxLayout(page)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setObjectName("Tiles")
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setFixedHeight(356)
        scroll.setMinimumWidth(TILE_COLUMNS * 108 + 18)
        inner = QWidget()
        inner.setObjectName("TilesInner")
        scroll.setWidget(inner)
        outer.addWidget(scroll)
        grid = QGridLayout(inner)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(4)
        for n, frag in enumerate(fragments):
            b = self._tile(frag, frag.name.replace(" (chair)", ""), QSize(66, 50), QSize(104, 82))
            b.setToolTip(frag.name)
            grid.addWidget(b, n // TILE_COLUMNS, n % TILE_COLUMNS)
            self._tiles.append((b, frag))
        grid.setRowStretch(len(fragments) // TILE_COLUMNS + 1, 1)

    def _tile(self, frag: Fragment, text: str, icon_size: QSize, size: QSize | None = None) -> QToolButton:
        b = button("Tile", text, checkable=True)
        b.setCursor(Qt.CursorShape.ArrowCursor)
        b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        b.setIconSize(icon_size)
        if size is not None:
            b.setFixedSize(size)
        if self.theme is not None:
            b.setIcon(QIcon(render_pixmap(frag.mol, icon_size.width(), icon_size.height(), self.theme)))
        b.clicked.connect(lambda _=False, f=frag: self._choose(f))
        return b

    def add_fragments(self, category: str, fragments: list[Fragment]) -> None:
        """Add fragments to a tab: a bundled category's (``"bio"``, ``"ligand"``) or one titled ``category``.

        A tab that does not exist yet is added after the others. Tabs
        already shown get the fragments when the popover next opens.
        """
        title = CATEGORY_TABS.get(category, category)
        index = next((k for k, (name, _) in enumerate(self._tabs) if name == title), None)
        if index is None:
            index = len(self._tabs)
            self._tabs.append((title, list))
            self.tabs.add(title)
            self.pages.addWidget(QWidget())
        self._added.setdefault(index, []).extend(fragments)
        self._forget_page(index)   # built again with the new fragments when shown

    def remove_fragments(self, category: str, fragments: list[Fragment]) -> None:
        """Take fragments added by :meth:`add_fragments` out again; a tab left empty goes too."""
        title = CATEGORY_TABS.get(category, category)
        index = next((k for k, (name, _) in enumerate(self._tabs) if name == title), None)
        if index is None or index not in self._added:
            return
        ids = {id(f) for f in fragments}
        self._added[index] = [f for f in self._added[index] if id(f) not in ids]
        self._forget_page(index)
        if not self._added[index] and index == len(self._tabs) - 1 and index >= len(TABS):
            del self._added[index]
            self._tabs.pop()
            self.tabs.remove_last()
            page = self.pages.widget(index)
            if page is not None:
                self.pages.removeWidget(page)
                page.deleteLater()
            if self.pages.currentIndex() >= len(self._tabs):
                self.tabs.set_index(0)
                self._show_tab(0)

    def _forget_page(self, index: int) -> None:
        """Build a tab again the next time it is shown, with the fragments it has then."""
        if index not in self._filled:
            return
        self._filled.discard(index)
        old = self.pages.widget(index)
        if old is not None:
            self._tiles = [(b, f) for b, f in self._tiles if not old.isAncestorOf(b)]
            self.pages.removeWidget(old)
            old.deleteLater()
        self.pages.insertWidget(index, QWidget())

    def _show_tab(self, index: int) -> None:
        make = self._tabs[index][1]
        page = self.pages.widget(index)
        if make is not None and page is not None and index not in self._filled:
            self._filled.add(index)
            self._fill_tiles(page, make() + self._added.get(index, []))
        self._refresh_thumbnails()
        self.pages.setCurrentIndex(index)
        self.sync()

    def _refresh_thumbnails(self) -> None:
        if self.theme is None or self._thumbnail_theme == self.theme.name:
            return
        self._thumbnail_theme = self.theme.name
        self._fill_variants()
        for b, frag in self._tiles:
            b.setIcon(QIcon(render_pixmap(frag.mol, 66, 50, self.theme)))

    def _fill_variants(self) -> None:
        """Buttons for the bonding types of the current element."""
        while self.variant_row.count():
            item = self.variant_row.takeAt(0)
            widget = item.widget() if item is not None else None
            if widget is not None:
                widget.deleteLater()
        frag = self.builder.fragment
        symbol = frag.symbol if frag.category == Category.ELEMENT else self.builder.last_element
        types = element_types(symbol)
        self.variant_title.setText(f"{element(symbol).name} bonding types{SEPARATOR}keys 1-{len(types)}")
        for n, kind in enumerate(types, start=1):
            variant = element_fragment(symbol, kind)
            b = self._tile(variant, f"{n}  {SHORT_LABELS.get(kind.label, kind.label.split(' (')[0])}",
                           QSize(50, 38))
            b.setToolTip(variant.name)
            b.setChecked(frag.key == variant.key)
            self.variant_row.addWidget(b)
        self.variant_row.addStretch(1)

    def _choose(self, frag: Fragment) -> None:
        self.builder.set_fragment(frag)
        self.close()

    # -- search -----------------------------------------------------------------------------

    def look_up(self) -> None:
        """Build or download what is typed in the search box."""
        text = self.search.text().strip()
        if not text:
            return
        source = choose_source(text, self.lookup.sources)
        if source is not None and source.message is not None:
            self.hint.setText(source.message(text))
        else:
            self.hint.setText(f"Building {text} ...")
        self.search.setEnabled(False)
        self.lookup.start(text)

    def search_for(self, text: str) -> None:
        """Put text in the search box and look it up."""
        self.search.setText(text)
        self.look_up()

    def _found_fragment(self, mol: Molecule, attach: int | None, smiles: str, text: str) -> None:
        self.search.setEnabled(True)
        if text != smiles:
            name = text.strip().capitalize()
        else:
            name = smiles if len(smiles) <= 28 else smiles[:26] + "..."
        self.builder.set_fragment(molecule_fragment(mol, name, smiles, attach))
        self.search.clear()
        self.close()
        self.placeRequested.emit()

    def _found_structure(self, mol: Molecule, info: JobInfo) -> None:
        self.search.setEnabled(True)
        self.search.clear()
        self.close()
        self.structureLoaded.emit(mol, info)

    def _failed(self, message: str) -> None:
        self.search.setEnabled(True)
        self.search.setFocus()
        self.hint.setText(message)

    def focus_search(self) -> None:
        """Put the cursor in the search box."""
        if self.search.isEnabled():
            self.search.setFocus()
            self.search.selectAll()

    # -- state ------------------------------------------------------------------------------

    def set_theme(self, theme: Theme) -> None:
        """Redraw thumbnails in a new theme."""
        self.theme = theme
        self.preview.theme = theme
        self._refresh_thumbnails()

    def sync(self) -> None:
        """Show the current fragment: checked cell or tile, bonding types, preview and hint."""
        frag = self.builder.fragment
        element_chosen = frag.category == Category.ELEMENT
        for symbol, b in self.cells.items():
            b.setChecked(element_chosen and symbol == frag.symbol)
        for b, f in self._tiles:
            b.setChecked(f.key == frag.key)
        self._fill_variants()
        self.preview.sync()
        self.name.setText(frag.name)
        if frag.attach is None:
            self.hint.setText("Single atom. Drag to rotate.")
        else:
            anchor = f"{frag.mol.atoms[frag.anchor].symbol}{frag.anchor + 1}"
            self.hint.setText(f"Bonds through {anchor} (green atom is replaced). "
                              "Click another end atom to change; drag to rotate.")

    def show_at(self, pos: QPoint) -> None:
        """Open at a screen position."""
        self.search.setPlaceholderText(search_hint(self.lookup.sources))
        current = self.pages.currentIndex()
        if current in self._added and current not in self._filled:   # a library added to the tab shown
            self._show_tab(current)
        self.sync()
        self.adjustSize()
        self.move(pos)
        self.show()
