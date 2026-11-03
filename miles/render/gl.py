"""OpenGL renderer: atoms and bonds as ray-cast impostors, cartoons as meshes.

Each atom is one instanced quad whose fragment shader solves the sphere
intersection per pixel, and each bond one quad solving the cylinder
intersection. Rotating or zooming only changes two matrices, so frames stay
cheap whatever the molecule size. Pictures for export are drawn into a
multisampled framebuffer, optionally followed by an ambient occlusion pass.
"""

import ctypes
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
import OpenGL

OpenGL.ERROR_CHECKING = False
from OpenGL import GL  # noqa: E402
from OpenGL.GL import shaders as gl_shaders  # noqa: E402

from miles.render import shaders  # noqa: E402
from miles.render.camera import Camera, DepthRange  # noqa: E402
from miles.render.cartoon import Mesh  # noqa: E402
from miles.render.theme import Theme, rgb  # noqa: E402

# Instance attributes: (location, number of floats, byte offset), matching the scene arrays.
# The highlight of each sphere (the last column of the scene array) comes from a buffer of its
# own, one float per sphere, so a hover change does not send every sphere to the GPU again.
SPHERE_ATTRIBUTES = [(1, 3, 0), (2, 1, 12), (3, 3, 16)]
SPHERE_STRIDE = 8 * 4
FLAG_LOCATION = 4
PARTIAL_FLAG_ROWS = 64   # up to this many changed highlights are written row by row
CYLINDER_ATTRIBUTES = [(1, 3, 0), (2, 3, 12), (3, 1, 24), (4, 3, 28), (5, 3, 40), (6, 1, 52)]
CYLINDER_STRIDE = 14 * 4
MESH_STRIDE = 9 * 4  # position, normal, colour

FADE = 0.32
OUTLINE_PX = 1.2
HALO_PX = 2.6
MESH_OUTLINE_PX = 1.3
AO_RADIUS = 1.4   # angstrom
AO_STRENGTH = 0.75
MAX_EXPORT_SAMPLES = 8

_SHARED_UNIFORMS = ("u_view", "u_proj", "u_px", "u_outline", "u_outline_color", "u_bg", "u_fade",
                    "u_zrange", "u_flag_colors", "u_halo", "u_dpr", "u_offset", "u_outline_pass", "u_alpha")


class MeshSlot(StrEnum):
    """Meshes the renderer keeps on the GPU."""

    CARTOON = "cartoon"   # opaque, with an outline
    ARROWS = "arrows"     # dipole and vibration arrows, opaque, with an outline like the cartoon


@dataclass(frozen=True, slots=True)
class Viewport:
    """Size of the drawing surface in logical pixels and its device pixel ratio."""

    width: float
    height: float
    dpr: float = 1.0


class _Program:
    """A linked shader program and the locations of its uniforms."""

    def __init__(self, vertex: str, fragment: str) -> None:
        self.id = gl_shaders.compileProgram(gl_shaders.compileShader(vertex, GL.GL_VERTEX_SHADER),
                                            gl_shaders.compileShader(fragment, GL.GL_FRAGMENT_SHADER),
                                            validate=False)
        self.uniforms = {name: GL.glGetUniformLocation(self.id, name) for name in _SHARED_UNIFORMS}

    def location(self, name: str) -> int:
        """Location of a uniform (-1 if the program does not use it)."""
        if name not in self.uniforms:
            self.uniforms[name] = GL.glGetUniformLocation(self.id, name)
        return self.uniforms[name]


@dataclass
class _MeshBuffers:
    """Vertex array, vertex buffer and index buffer of one mesh."""

    vao: int
    vbo: int
    ebo: int
    count: int = 0
    alpha: float = 1.0


def _mesh_buffers() -> _MeshBuffers:
    vao = GL.glGenVertexArrays(1)
    vbo, ebo = GL.glGenBuffers(2)
    GL.glBindVertexArray(vao)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vbo)
    for location, offset in ((0, 0), (1, 12), (2, 24)):
        GL.glEnableVertexAttribArray(location)
        GL.glVertexAttribPointer(location, 3, GL.GL_FLOAT, GL.GL_FALSE, MESH_STRIDE, ctypes.c_void_p(offset))
    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ebo)
    GL.glBindVertexArray(0)
    return _MeshBuffers(vao, vbo, ebo)


def _instanced_vao(corners: int, vbo: int, attributes: list[tuple[int, int, int]], stride: int) -> int:
    """Vertex array with the shared quad corners plus per-instance attributes."""
    vao = GL.glGenVertexArrays(1)
    GL.glBindVertexArray(vao)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, corners)
    GL.glEnableVertexAttribArray(0)
    GL.glVertexAttribPointer(0, 2, GL.GL_FLOAT, GL.GL_FALSE, 0, None)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vbo)
    for location, size, offset in attributes:
        GL.glEnableVertexAttribArray(location)
        GL.glVertexAttribPointer(location, size, GL.GL_FLOAT, GL.GL_FALSE, stride, ctypes.c_void_p(offset))
        GL.glVertexAttribDivisor(location, 1)
    GL.glBindVertexArray(0)
    return vao


def _check_framebuffer(what: str) -> None:
    """Raise if the bound framebuffer cannot be drawn into (too large, unsupported format)."""
    status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
    if status != GL.GL_FRAMEBUFFER_COMPLETE:
        raise RuntimeError(f"The graphics card cannot make a picture of {what} "
                           f"(framebuffer status {status:#x})")


def _texture(internal_format: int, fmt: int, kind: int, width: int, height: int) -> int:
    tex = GL.glGenTextures(1)
    GL.glBindTexture(GL.GL_TEXTURE_2D, tex)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, internal_format, width, height, 0, fmt, kind, None)
    for param, value in ((GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST),
                         (GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST),
                         (GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE),
                         (GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)):
        GL.glTexParameteri(GL.GL_TEXTURE_2D, param, value)
    return tex


class GLRenderer:
    """Draws scene arrays with OpenGL core profile shaders (GLSL 3.30, in a 4.1 context on macOS).

    Call :meth:`initialize` with a current context first; if it fails,
    :attr:`ok` stays false and :attr:`error` says why, and the caller should
    fall back to QPainter drawing. Call :meth:`release` with the same
    context current before the context goes away.
    """

    def __init__(self) -> None:
        self.ok = False
        self.error = ""
        self.n_spheres = 0
        self.n_cylinders = 0
        self._meshes: dict[MeshSlot, _MeshBuffers] = {}
        self._extra: list[_MeshBuffers] = []   # meshes from plugins, drawn like surfaces after the others
        self._ao: _Program | None = None
        self._quad_vao: int | None = None

    def initialize(self) -> bool:
        """Compile the shaders and create the buffers; returns :attr:`ok`."""
        try:
            self._sphere = _Program(shaders.SPHERE_VS, shaders.SPHERE_FS)
            self._cylinder = _Program(shaders.CYLINDER_VS, shaders.CYLINDER_FS)
            self._mesh = _Program(shaders.MESH_VS, shaders.MESH_FS)
            corners = np.array([-1, -1, 1, -1, -1, 1, 1, 1], np.float32)
            self._corners = GL.glGenBuffers(1)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._corners)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, corners.nbytes, corners, GL.GL_STATIC_DRAW)
            self._sphere_vbo, self._cylinder_vbo, self._flag_vbo = GL.glGenBuffers(3)
            self._sphere_vao = _instanced_vao(self._corners, self._sphere_vbo,
                                              SPHERE_ATTRIBUTES, SPHERE_STRIDE)
            GL.glBindVertexArray(self._sphere_vao)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._flag_vbo)
            GL.glEnableVertexAttribArray(FLAG_LOCATION)
            GL.glVertexAttribPointer(FLAG_LOCATION, 1, GL.GL_FLOAT, GL.GL_FALSE, 4, None)
            GL.glVertexAttribDivisor(FLAG_LOCATION, 1)
            GL.glBindVertexArray(0)
            self._cylinder_vao = _instanced_vao(self._corners, self._cylinder_vbo,
                                                CYLINDER_ATTRIBUTES, CYLINDER_STRIDE)
            self._meshes = {slot: _mesh_buffers() for slot in MeshSlot}
            self.ok = True
        except Exception as exc:  # shader compilation fails or the context is too old
            self.ok = False
            self.error = str(exc)
        return self.ok

    # -- buffers ---------------------------------------------------------------------------

    def upload_spheres(self, data: np.ndarray) -> None:
        """Replace the sphere instances (rows of ``Scene.spheres``), highlights included."""
        self.n_spheres = len(data)
        if self.n_spheres:
            data = np.ascontiguousarray(data, dtype=np.float32)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._sphere_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, data.nbytes, data, GL.GL_DYNAMIC_DRAW)
            self.upload_flags(data[:, 7])

    def upload_flags(self, flags: np.ndarray, rows: np.ndarray | None = None) -> None:
        """Replace the highlight of each sphere, or only of the sphere rows given.

        Parameters
        ----------
        flags : numpy.ndarray
            Highlight value of every sphere (see :class:`~miles.render.scene.Highlight`).
        rows : numpy.ndarray, optional
            The rows that changed since the last upload; few rows are
            written one by one, otherwise the whole buffer is replaced.
        """
        if not self.n_spheres:
            return
        flags = np.ascontiguousarray(flags, dtype=np.float32)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._flag_vbo)
        if rows is not None and len(rows) <= PARTIAL_FLAG_ROWS and len(flags) == self.n_spheres:
            for row in rows.tolist():
                GL.glBufferSubData(GL.GL_ARRAY_BUFFER, 4 * row, 4, flags[row:row + 1])
        else:
            GL.glBufferData(GL.GL_ARRAY_BUFFER, flags.nbytes, flags, GL.GL_DYNAMIC_DRAW)

    def upload_cylinders(self, data: np.ndarray) -> None:
        """Replace the cylinder instances (rows of ``Scene.cylinders``)."""
        self.n_cylinders = len(data)
        if self.n_cylinders:
            data = np.ascontiguousarray(data, dtype=np.float32)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._cylinder_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, data.nbytes, data, GL.GL_DYNAMIC_DRAW)

    def set_layer_alphas(self, alphas: list[float]) -> None:
        """Change the opacity of the plugins' meshes without uploading them again."""
        for buffers, alpha in zip(self._extra, alphas, strict=False):
            buffers.alpha = alpha

    def upload_mesh(self, mesh: Mesh | None, colors: np.ndarray | None = None,
                    slot: MeshSlot = MeshSlot.CARTOON, alpha: float = 1.0) -> None:
        """Replace a mesh, optionally with other vertex colours and an opacity."""
        self._upload(self._meshes[slot], mesh, colors, alpha)

    def upload_layers(self, meshes: list[tuple[Mesh, float]]) -> None:
        """Replace the meshes of plugins (surfaces, isosurfaces), drawn after everything opaque."""
        while len(self._extra) < len(meshes):
            self._extra.append(_mesh_buffers())
        for k, buffers in enumerate(self._extra):
            if k < len(meshes):
                self._upload(buffers, meshes[k][0], None, meshes[k][1])
            else:
                buffers.count = 0

    def _upload(self, buffers: _MeshBuffers, mesh: Mesh | None, colors: np.ndarray | None,
                alpha: float) -> None:
        buffers.alpha = alpha
        if mesh is None or mesh.empty:
            buffers.count = 0
            return
        colors = mesh.colors if colors is None else colors
        data = np.ascontiguousarray(np.hstack([mesh.vertices, mesh.normals, colors]), dtype=np.float32)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, buffers.vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, data.nbytes, data, GL.GL_DYNAMIC_DRAW)
        indices = np.ascontiguousarray(mesh.indices, dtype=np.uint32)
        GL.glBindVertexArray(buffers.vao)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, buffers.ebo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
        GL.glBindVertexArray(0)
        buffers.count = indices.size

    def release(self) -> None:
        """Delete the GPU objects; call with the renderer's context current, before it goes away."""
        if not self.ok:
            return
        self.ok = False
        meshes = [*self._meshes.values(), *self._extra]
        buffers = [self._corners, self._sphere_vbo, self._cylinder_vbo, self._flag_vbo,
                   *(b for m in meshes for b in (m.vbo, m.ebo))]
        GL.glDeleteBuffers(len(buffers), buffers)
        vaos = [self._sphere_vao, self._cylinder_vao, *(m.vao for m in meshes)]
        if self._quad_vao is not None:
            vaos.append(self._quad_vao)
        GL.glDeleteVertexArrays(len(vaos), vaos)
        for program in (self._sphere, self._cylinder, self._mesh, self._ao):
            if program is not None:
                GL.glDeleteProgram(program.id)
        self._meshes, self._extra, self._ao, self._quad_vao = {}, [], None, None
        self.n_spheres = self.n_cylinders = 0

    # -- drawing ---------------------------------------------------------------------------

    def draw(self, camera: Camera, viewport: Viewport, theme: Theme, fade: bool, depth: DepthRange) -> None:
        """Clear the current framebuffer to the theme background and draw the scene."""
        GL.glViewport(0, 0, round(viewport.width * viewport.dpr), round(viewport.height * viewport.dpr))
        GL.glClearColor(*rgb(theme.background), 1.0)
        GL.glClearDepth(1.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        self._draw_scene(camera, viewport, theme, fade, depth)

    def _draw_scene(self, camera: Camera, viewport: Viewport, theme: Theme, fade: bool,
                    depth: DepthRange) -> None:
        opaque = [self._meshes[MeshSlot.CARTOON], self._meshes[MeshSlot.ARROWS]]
        layers = self._extra
        if not (self.n_spheres or self.n_cylinders or any(m.count for m in opaque + layers)):
            return
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthFunc(GL.GL_LESS)
        GL.glDepthMask(GL.GL_TRUE)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_CULL_FACE)
        GL.glEnable(GL.GL_MULTISAMPLE)
        GL.glEnable(GL.GL_SAMPLE_ALPHA_TO_COVERAGE)

        dpr = viewport.dpr
        view = camera.view_matrix().astype(np.float32)
        proj = camera.projection_matrix(viewport.width, viewport.height, depth.extent).astype(np.float32)
        px = 1.0 / (camera.scale * dpr)
        highlight_colors = np.array([rgb(c) for c in theme.highlight_colors()], np.float32)

        def use(program: _Program) -> dict[str, int]:
            u = program.uniforms
            GL.glUseProgram(program.id)
            GL.glUniformMatrix4fv(u["u_view"], 1, GL.GL_TRUE, view)
            GL.glUniformMatrix4fv(u["u_proj"], 1, GL.GL_TRUE, proj)
            GL.glUniform1f(u["u_px"], px)
            GL.glUniform1f(u["u_outline"], OUTLINE_PX * dpr)
            GL.glUniform3f(u["u_outline_color"], *rgb(theme.outline))
            GL.glUniform3f(u["u_bg"], *rgb(theme.background))
            GL.glUniform1f(u["u_fade"], FADE if fade else 0.0)
            GL.glUniform2f(u["u_zrange"], depth.near, depth.far)
            GL.glUniform1f(u["u_dpr"], dpr)
            GL.glUniform3fv(u["u_flag_colors"], len(highlight_colors), highlight_colors)
            GL.glUniform1f(u["u_halo"], HALO_PX * dpr)
            GL.glUniform1f(u["u_alpha"], 1.0)
            return u

        for mesh in opaque:
            if not mesh.count:
                continue
            # back faces pushed outwards along the normals give the outline, then the surface
            u = use(self._mesh)
            GL.glBindVertexArray(mesh.vao)
            GL.glEnable(GL.GL_CULL_FACE)
            GL.glCullFace(GL.GL_FRONT)
            GL.glUniform1f(u["u_outline_pass"], 1.0)
            GL.glUniform1f(u["u_offset"], MESH_OUTLINE_PX * dpr * px)
            GL.glDrawElements(GL.GL_TRIANGLES, mesh.count, GL.GL_UNSIGNED_INT, None)
            GL.glDisable(GL.GL_CULL_FACE)
            GL.glUniform1f(u["u_outline_pass"], 0.0)
            GL.glUniform1f(u["u_offset"], 0.0)
            GL.glDrawElements(GL.GL_TRIANGLES, mesh.count, GL.GL_UNSIGNED_INT, None)

        for program, vao, count in ((self._cylinder, self._cylinder_vao, self.n_cylinders),
                                    (self._sphere, self._sphere_vao, self.n_spheres)):
            if count:
                use(program)
                GL.glBindVertexArray(vao)
                GL.glDrawArraysInstanced(GL.GL_TRIANGLE_STRIP, 0, 4, count)

        for layer in layers:
            if layer.count:
                self._draw_surface(use, layer)

        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_SAMPLE_ALPHA_TO_COVERAGE)
        GL.glDisable(GL.GL_DEPTH_TEST)

    def _draw_surface(self, use: Callable[[_Program], dict[str, int]], surface: _MeshBuffers) -> None:
        """Draw a surface after everything opaque; see-through surfaces are blended."""
        u = use(self._mesh)
        GL.glUniform1f(u["u_outline_pass"], 0.0)
        GL.glUniform1f(u["u_offset"], 0.0)
        GL.glUniform1f(u["u_alpha"], surface.alpha)
        GL.glBindVertexArray(surface.vao)
        if surface.alpha < 1.0:
            # front faces only, blended over what is behind, without hiding it from the depth test
            GL.glDisable(GL.GL_SAMPLE_ALPHA_TO_COVERAGE)
            GL.glEnable(GL.GL_BLEND)
            GL.glBlendFuncSeparate(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA,
                                   GL.GL_ONE, GL.GL_ONE_MINUS_SRC_ALPHA)
            GL.glDepthMask(GL.GL_FALSE)
            GL.glEnable(GL.GL_CULL_FACE)
            GL.glCullFace(GL.GL_BACK)
        GL.glDrawElements(GL.GL_TRIANGLES, surface.count, GL.GL_UNSIGNED_INT, None)
        GL.glDisable(GL.GL_CULL_FACE)
        GL.glDisable(GL.GL_BLEND)
        GL.glDepthMask(GL.GL_TRUE)

    # -- export ----------------------------------------------------------------------------

    def render_image(self, camera: Camera, width: float, height: float, scale: float, theme: Theme,
                     fade: bool, depth: DepthRange, transparent: bool, ambient_occlusion: bool) -> np.ndarray:
        """Draw the scene off screen for export.

        Parameters
        ----------
        camera : Camera
            View to draw.
        width, height : float
            View size in logical pixels.
        scale : float
            Pixels per logical pixel in the picture (capped by the GPU limits).
        theme : Theme
            Colours.
        fade : bool
            Depth fading.
        depth : DepthRange
            Scene depth.
        transparent : bool
            Transparent instead of theme background.
        ambient_occlusion : bool
            Darken creases with a screen-space ambient occlusion pass.

        Returns
        -------
        numpy.ndarray
            RGBA picture with premultiplied alpha, shape ``(H, W, 4)``, top row first.
        """
        limit = int(GL.glGetIntegerv(GL.GL_MAX_RENDERBUFFER_SIZE))
        scale = min(scale, limit / max(width, height))
        W, H = round(width * scale), round(height * scale)
        samples = min(MAX_EXPORT_SAMPLES, int(GL.glGetIntegerv(GL.GL_MAX_SAMPLES)))

        framebuffers: list[int] = []
        renderbuffers: list[int] = []
        textures: list[int] = []
        try:   # the GPU objects are deleted whatever happens while drawing
            multisampled = GL.glGenFramebuffers(1)
            framebuffers.append(multisampled)
            renderbuffers.extend(int(rb) for rb in GL.glGenRenderbuffers(2))
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, multisampled)
            for rb, fmt, attachment in ((renderbuffers[0], GL.GL_RGBA8, GL.GL_COLOR_ATTACHMENT0),
                                        (renderbuffers[1], GL.GL_DEPTH_COMPONENT24, GL.GL_DEPTH_ATTACHMENT)):
                GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, rb)
                GL.glRenderbufferStorageMultisample(GL.GL_RENDERBUFFER, samples, fmt, W, H)
                GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, attachment, GL.GL_RENDERBUFFER, rb)
            _check_framebuffer(f"{W} x {H} with {samples} samples")
            GL.glViewport(0, 0, W, H)
            if transparent:
                GL.glClearColor(0.0, 0.0, 0.0, 0.0)
            else:
                GL.glClearColor(*rgb(theme.background), 1.0)
            GL.glClearDepth(1.0)
            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
            GL.glEnable(GL.GL_SAMPLE_ALPHA_TO_ONE)
            self._draw_scene(camera, Viewport(width, height, scale), theme, fade, depth)
            GL.glDisable(GL.GL_SAMPLE_ALPHA_TO_ONE)

            # resolve the samples into textures the occlusion pass can read
            color = _texture(GL.GL_RGBA8, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, W, H)
            textures.append(color)
            depth_tex = _texture(GL.GL_DEPTH_COMPONENT24, GL.GL_DEPTH_COMPONENT, GL.GL_FLOAT, W, H)
            textures.append(depth_tex)
            resolved = GL.glGenFramebuffers(1)
            framebuffers.append(resolved)
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, resolved)
            texture_2d = GL.GL_TEXTURE_2D
            GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, texture_2d, color, 0)
            GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT, texture_2d, depth_tex, 0)
            _check_framebuffer(f"{W} x {H} for reading back")
            GL.glBindFramebuffer(GL.GL_READ_FRAMEBUFFER, multisampled)
            GL.glBindFramebuffer(GL.GL_DRAW_FRAMEBUFFER, resolved)
            GL.glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL.GL_COLOR_BUFFER_BIT, GL.GL_NEAREST)
            GL.glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL.GL_DEPTH_BUFFER_BIT, GL.GL_NEAREST)

            source = resolved
            if ambient_occlusion:
                out = _texture(GL.GL_RGBA8, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, W, H)
                textures.append(out)
                source = GL.glGenFramebuffers(1)
                framebuffers.append(source)
                GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, source)
                GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, texture_2d, out, 0)
                _check_framebuffer(f"{W} x {H} for ambient occlusion")
                self._occlusion_pass(color, depth_tex, W, H, AO_RADIUS * camera.scale * scale, depth.extent)

            GL.glBindFramebuffer(GL.GL_READ_FRAMEBUFFER, source)
            GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
            raw = GL.glReadPixels(0, 0, W, H, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
            image = np.frombuffer(raw, dtype=np.uint8).reshape(H, W, 4)[::-1].copy()
        finally:
            GL.glDisable(GL.GL_SAMPLE_ALPHA_TO_ONE)
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)
            if framebuffers:
                GL.glDeleteFramebuffers(len(framebuffers), framebuffers)
            if renderbuffers:
                GL.glDeleteRenderbuffers(len(renderbuffers), renderbuffers)
            if textures:
                GL.glDeleteTextures(textures)
        return image

    def _occlusion_pass(self, color: int, depth: int, width: int, height: int, radius_px: float,
                        depth_extent: float) -> None:
        """Draw the colour texture darkened by screen-space ambient occlusion."""
        if self._ao is None:
            self._ao = _Program(shaders.AO_VS, shaders.AO_FS)
            self._quad_vao = GL.glGenVertexArrays(1)
            GL.glBindVertexArray(self._quad_vao)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._corners)
            GL.glEnableVertexAttribArray(0)
            GL.glVertexAttribPointer(0, 2, GL.GL_FLOAT, GL.GL_FALSE, 0, None)
            GL.glBindVertexArray(0)
        ao = self._ao
        GL.glUseProgram(ao.id)
        GL.glActiveTexture(GL.GL_TEXTURE0)
        GL.glBindTexture(GL.GL_TEXTURE_2D, color)
        GL.glActiveTexture(GL.GL_TEXTURE1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, depth)
        GL.glUniform1i(ao.location("u_color"), 0)
        GL.glUniform1i(ao.location("u_depth"), 1)
        GL.glUniform2f(ao.location("u_texel"), 1.0 / width, 1.0 / height)
        GL.glUniform1f(ao.location("u_radius_px"), radius_px)
        GL.glUniform1f(ao.location("u_depth_range"), depth_extent)
        GL.glUniform1f(ao.location("u_world_radius"), AO_RADIUS)
        GL.glUniform1f(ao.location("u_strength"), AO_STRENGTH)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glBindVertexArray(self._quad_vao)
        GL.glDrawArrays(GL.GL_TRIANGLE_STRIP, 0, 4)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glActiveTexture(GL.GL_TEXTURE0)
