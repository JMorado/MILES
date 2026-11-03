"""Mouse, trackpad and keyboard input of a canvas, turned into view moves, tool actions and commands."""

import math
from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtGui import QHelpEvent, QKeyEvent, QMouseEvent, QNativeGestureEvent, QWheelEvent
from PySide6.QtWidgets import QAbstractSpinBox, QApplication, QLineEdit, QToolTip

from miles.core import building
from miles.core.fragments import element_types
from miles.core.geometry import rotation_matrix
from miles.core.internal import rotate_atoms, translate_atoms
from miles.core.molecule import Molecule
from miles.render.camera import MAX_SCALE, MIN_SCALE
from miles.render.picking import AtomHit, BondHit, first_atom
from miles.ui.builder import MouseMode, Tool
from miles.ui.canvas.context_menu import show_context_menu
from miles.ui.canvas.overlay import annotation_tip
from miles.ui.canvas.tools import Drag, Gesture, finish_draw
from miles.ui.document import Change
from miles.ui.text import SEPARATOR

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

DRAG_THRESHOLD = 3          # pixels before a press becomes a drag

ROTATE_SPEED = 0.0095       # radians per pixel

KEY_TURN = math.radians(6)

KEY_ZOOM = 1.15

WHEEL_ZOOM = 0.0016

DRAG_ZOOM = 0.006

MODE_DRAGS = {MouseMode.ROTATE: Drag.ROTATE, MouseMode.TRANSLATE: Drag.PAN, MouseMode.SCALE: Drag.SCALE}


class InputController:
    """What the mouse, the trackpad and the keyboard do on one canvas.

    The canvas passes its events here; presses become a
    :class:`~miles.ui.canvas.tools.Gesture` that the active tool or a view
    move acts on, and keys choose elements, frames and commands.

    Parameters
    ----------
    canvas : Canvas
        The canvas whose input this is.
    """

    def __init__(self, canvas: "Canvas") -> None:
        self.canvas = canvas

    def enter(self) -> None:
        """Take the keyboard when the cursor comes over the canvas, unless the user is typing elsewhere."""
        canvas = self.canvas
        focus = QApplication.focusWidget()
        if canvas.window().isActiveWindow() and not isinstance(focus, (QLineEdit, QAbstractSpinBox)):
            canvas.setFocus()

    def leave(self) -> None:
        """Drop the highlight under the cursor."""
        if self.canvas.hover is not None:
            self.canvas.hover = None
            self.canvas.redraw()

    def press(self, e: QMouseEvent) -> None:
        """Start a gesture: a tool action, a view move or a drag of atoms."""
        canvas = self.canvas
        canvas.setFocus()
        canvas.focused.emit()
        canvas.builder.end_typing()
        canvas.view.stop()
        pos = e.position()
        gesture = canvas.gesture = Gesture(press=pos, last=pos, button=e.button())
        mods = e.modifiers()
        shift = bool(mods & Qt.KeyboardModifier.ShiftModifier)
        if e.button() == Qt.MouseButton.LeftButton:
            hit = canvas.pick(pos)
            if mods & Qt.KeyboardModifier.AltModifier:
                if canvas.selection and canvas.builder.tool == Tool.SELECT:
                    gesture.atoms = list(canvas.selection)
                elif hit is not None:
                    gesture.atoms = sorted(canvas.mol.component(first_atom(hit)))
                gesture.drag = Drag.MOVE_ATOMS if shift else Drag.TURN_ATOMS
            elif shift:
                gesture.drag = Drag.PAN
            elif mods & Qt.KeyboardModifier.ControlModifier:
                gesture.drag = Drag.SPIN
            elif not canvas.tool.press(canvas, gesture, hit):
                gesture.drag = MODE_DRAGS[canvas.builder.mouse_mode]
        elif e.button() == Qt.MouseButton.MiddleButton:
            gesture.drag = Drag.PAN
        elif e.button() == Qt.MouseButton.RightButton:
            gesture.drag = Drag.ZOOM

    def move(self, e: QMouseEvent) -> None:
        """Follow the cursor: the hover highlight, or the drag in progress."""
        canvas = self.canvas
        pos = e.position()
        g = canvas.gesture
        if g is None:
            hit = canvas.pick(pos)
            if hit != canvas.hover:
                canvas.hover = hit
                canvas.redraw()
            return
        if not g.moved:
            if (pos - g.press).manhattanLength() < DRAG_THRESHOLD:
                return
            g.moved = True
            canvas.hover = None
        dx, dy = pos.x() - g.last.x(), pos.y() - g.last.y()
        g.last = pos
        cam = canvas.camera
        match g.drag:
            case Drag.DRAW:
                hit = canvas.pick(pos)
                g.draw_to = hit.index if isinstance(hit, AtomHit) and hit.index != g.draw_from else None
                canvas.hover = None if g.draw_to is None else AtomHit(g.draw_to)
            case Drag.ROTATE:
                cam.turn(dy * ROTATE_SPEED, dx * ROTATE_SPEED)
            case Drag.SPIN:
                cx, cy = canvas.width() / 2 + cam.pan[0], canvas.height() / 2 + cam.pan[1]
                before = math.atan2(pos.y() - dy - cy, pos.x() - dx - cx)
                cam.turn(0, 0, -(math.atan2(pos.y() - cy, pos.x() - cx) - before))
            case Drag.PAN:
                cam.pan += (dx, dy)
            case Drag.ZOOM:
                cam.scale = max(MIN_SCALE, min(MAX_SCALE, cam.scale * math.exp(-dy * DRAG_ZOOM)))
            case Drag.SCALE:
                cam.scale = max(MIN_SCALE, min(MAX_SCALE, cam.scale * math.exp((dx - dy) * DRAG_ZOOM)))
            case Drag.MOVE_ATOMS | Drag.TURN_ATOMS:
                self._drag_atoms(g, dx, dy)
        canvas.redraw()

    def _drag_atoms(self, g: Gesture, dx: float, dy: float) -> None:
        canvas = self.canvas
        cam = canvas.camera
        if g.atoms is None:
            cam.turn(dy * ROTATE_SPEED, dx * ROTATE_SPEED)
            return
        if not g.checkpointed:
            canvas.doc.checkpoint()
            g.checkpointed = True
        if g.drag == Drag.MOVE_ATOMS:
            translate_atoms(canvas.mol, g.atoms, cam.screen_vector_to_world(dx, dy))
        else:
            turn = dx * ROTATE_SPEED, dy * ROTATE_SPEED
            screen = rotation_matrix([0, 1, 0], turn[0]) @ rotation_matrix([1, 0, 0], turn[1])
            rotate_atoms(canvas.mol, g.atoms, cam.rotation.T @ screen @ cam.rotation)
        canvas.doc.notify(Change.COORDS)

    def release(self, e: QMouseEvent) -> None:
        """End a gesture: finish a drag, or act on a click (tool or right-click menu)."""
        canvas = self.canvas
        g = canvas.gesture
        if g is None:
            return
        canvas.gesture = None
        pos = e.position()
        mods = e.modifiers()
        if g.drag == Drag.DRAW and g.moved:
            finish_draw(canvas, g, pos)
        elif not g.moved:
            menu = g.button == Qt.MouseButton.RightButton or (
                g.button == Qt.MouseButton.LeftButton and mods & Qt.KeyboardModifier.MetaModifier)
            if menu:
                show_context_menu(canvas, pos, canvas.pick(pos))
            elif g.button == Qt.MouseButton.LeftButton:
                adding = Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier
                canvas.tool.click(canvas, canvas.pick(pos), pos, additive=bool(mods & adding))
        canvas.hover = canvas.pick(pos)
        canvas.redraw()

    def double_click(self, e: QMouseEvent) -> bool:
        """Let the tool act on a double click; returns whether it starts a new press."""
        canvas = self.canvas
        return e.button() == Qt.MouseButton.LeftButton and bool(canvas.tool.double_click(canvas,
                                                                                          canvas.pick(e.position())))

    def wheel(self, e: QWheelEvent) -> None:
        """Zoom, or with Option over a selection of 2-4 atoms, turn the selection panel's dial."""
        canvas = self.canvas
        delta = e.angleDelta().y() or e.angleDelta().x() or (e.pixelDelta().y() or e.pixelDelta().x()) * 2
        if not delta:
            return
        if e.modifiers() & Qt.KeyboardModifier.AltModifier and canvas.builder.tool == Tool.SELECT \
                and 2 <= len(canvas.selection) <= 4:
            canvas.nudge.emit(delta / 120.0)
            return
        canvas.view.zoom_at(e.position(), math.exp(delta * WHEEL_ZOOM))

    def tooltip(self, e: QHelpEvent) -> None:
        """Show the tip of a plugin's line or label under the cursor."""
        tip = annotation_tip(self.canvas, e.pos().x(), e.pos().y())
        if tip:
            QToolTip.showText(e.globalPos(), tip, self.canvas)
        else:
            QToolTip.hideText()

    def native_gesture(self, gesture: QNativeGestureEvent) -> bool:
        """Pinch to zoom and rotate with two fingers; returns whether the gesture was used."""
        canvas = self.canvas
        if gesture.gestureType() == Qt.NativeGestureType.ZoomNativeGesture:
            canvas.view.zoom_at(gesture.position(), 1.0 + gesture.value())
            return True
        if gesture.gestureType() == Qt.NativeGestureType.RotateNativeGesture:
            canvas.camera.turn(0, 0, -math.radians(gesture.value()))
            canvas.redraw()
            return True
        return False

    def key(self, e: QKeyEvent) -> bool:
        """Act on a key (letters, bonding types, frames, turning, zoom); returns whether it was used."""
        canvas = self.canvas
        key = e.key()
        text = e.text()
        mods = e.modifiers() & ~Qt.KeyboardModifier.KeypadModifier
        plain = mods in (Qt.KeyboardModifier.NoModifier, Qt.KeyboardModifier.ShiftModifier)
        turns = {Qt.Key.Key_Left: (0, -KEY_TURN), Qt.Key.Key_Right: (0, KEY_TURN),
                 Qt.Key.Key_Up: (-KEY_TURN, 0), Qt.Key.Key_Down: (KEY_TURN, 0)}
        doc = canvas.doc
        if plain and len(text) == 1 and text.isalpha():
            self._type_letter(text)
        elif plain and len(text) == 1 and text.isdigit() and text != "0":
            frag = canvas.builder.choose_bonding_type(int(text))
            if frag is not None:
                canvas.toast(frag.name)
        elif key == Qt.Key.Key_Escape:
            if canvas.selection:
                canvas.clear_selection()
            else:
                canvas.builder.set_tool(Tool.SELECT)
        elif key in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self._delete_under_cursor()
        elif key in (Qt.Key.Key_Comma, Qt.Key.Key_Period) and doc.frame_count > 1:
            doc.set_frame(doc.frame + (1 if key == Qt.Key.Key_Period else -1))
        elif key == Qt.Key.Key_Space and doc.frame_count > 1:
            canvas.playRequested.emit()
        elif key in turns:
            canvas.camera.turn(*turns[Qt.Key(key)])
        elif key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal):
            canvas.camera.scale *= KEY_ZOOM
        elif key == Qt.Key.Key_Minus:
            canvas.camera.scale /= KEY_ZOOM
        else:
            return False
        canvas.redraw()
        return True

    def _type_letter(self, letter: str) -> None:
        canvas = self.canvas
        frag = canvas.builder.type_letter(letter)
        if frag is not None:
            canvas.toast(f"{frag.name}{SEPARATOR}1-{len(element_types(frag.symbol))} for other types")

    def _delete_under_cursor(self) -> None:
        canvas = self.canvas
        if canvas.selection:
            canvas.commands.delete_selection()
        elif isinstance(canvas.hover, AtomHit):
            canvas.doc.edit(building.delete_atom, canvas.hover.index)
        elif isinstance(canvas.hover, BondHit):
            canvas.doc.edit(Molecule.remove_bond, canvas.hover.i, canvas.hover.j)
        canvas.hover = None
