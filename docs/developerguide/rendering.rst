Rendering
#########

The ``render`` layer turns a molecule into pixels. It is organised so that
the expensive work, walking every atom and bond, happens only when the
structure or the style changes, while rotating, zooming and hovering only
change a few numbers.

From molecule to scene
**********************

A :class:`~miles.render.scene.Style` holds the display options of one
view: the default representation, whether hydrogens are shown, labels,
depth fading, the cartoon colouring and the charge displays.
:func:`~miles.render.scene.build_scene` combines a molecule, a style
and a :class:`~miles.render.theme.Theme` into a
:class:`~miles.render.scene.Scene`: a few flat float32 arrays that both
renderers consume.

``spheres``, one row of 8 values per drawn atom
   centre (3), radius, colour (3) and highlight. A negative radius is a
   width in pixels, used for wire atoms, which keep the same thickness at any
   zoom.

``cylinders``, one row of 14 values per bond line
   both ends (6), radius, the colours of the two halves (6) and a kind:
   solid, dashed or cartoon trace.

``sphere_atoms`` and ``cylinder_trim``
   the atom behind each sphere row, and the radii of the atoms at both ends
   of each cylinder.

``mesh``
   the cartoon, if any atoms are drawn as one.

Each atom's representation is its own override or the view's default.
Atoms in a cartoon or hidden get no sphere. A bond is drawn when both its
atoms are drawn and neither is in space-fill. Bonds between two
ball-and-stick atoms show their order: a double bond is two thinner
cylinders, a triple bond three, an aromatic bond one solid and one dashed
cylinder. The lines are spread in the plane of a neighbouring atom, and for
aromatic bonds towards the next ring bond, so the dashed line falls inside
the ring. Bonds involving tubes or wires are drawn as one cylinder.

Everything a scene takes from the atoms, bonds, display and theme, but not
from the positions, is worked out once by
:func:`~miles.render.scene.scene_layout`: representations, visibility,
radii, colours, cartoon residues and a
:class:`~miles.render.scene.BondLayout`. The bond layout holds every
cylinder row with its radius, colours and kind, which atoms it joins, how
far it is moved off the bond axis, and for each multiple bond the
neighbours whose plane may spread its lines, in order of preference. The
choice between those neighbours depends on the positions only in one case,
a neighbour lying on the bond axis, so :func:`~miles.render.scene.bond_cylinders`
settles it with array operations for each frame. A trajectory frame or a
drag step therefore costs a few milliseconds whatever the number of double
and aromatic bonds (about 3 ms for 9,600 atoms with 6,400 aromatic bonds).

Atom colours come from the theme's element colours, or from colours passed
to :func:`~miles.render.scene.build_scene` (``colors``), which plugins
give as scene items; the Charges plugin passes
:func:`~miles.render.scene.charge_colors`, the red-white-blue scale of
the charge-coloured surface, symmetric about zero. Bond halves take the
colours of their atoms either way.

The highlight column selects a colour from the theme: selected, hovered,
attachment point or about to be deleted.
:meth:`Scene.flags <miles.render.scene.Scene.flags>` gives the highlight
of each sphere row as an array of its own, which is all that changes while
the mouse moves; :meth:`Scene.highlighted <miles.render.scene.Scene.highlighted>`
returns a copy of the sphere array with that column set, for the QPainter
renderer.

The camera
**********

:class:`~miles.render.camera.Camera` is orthographic. It stores a
rotation matrix, a centre, a scale in pixels per angstrom and a pan offset in
pixels. Orthographic projection suits chemistry: sizes and angles do not
depend on depth, so a bond looks equally long anywhere on the screen.

:meth:`~miles.render.camera.Camera.turn` composes a rotation about the
screen axes with the current rotation and re-orthogonalises the result with
a singular value decomposition, so rounding errors cannot accumulate into a
skewed matrix after thousands of small turns.

Animated camera moves (fit, centre, look along an axis) interpolate between
two cameras with :func:`~miles.render.camera.interpolate`: rotations
along the shortest arc between their quaternions, zoom geometrically and
everything else linearly, with an ease-out curve over a quarter of a second.

The OpenGL renderer
*******************

:class:`~miles.render.gl.GLRenderer` draws the scene with OpenGL 3.3
core profile shaders (GLSL 3.30), in the OpenGL 4.1 core context that
:func:`~miles.ui.canvas.surfaces.configure_gl` asks for, the newest macOS
offers; the shader sources are in ``render/shaders.py``.

Impostors
=========

A sphere is not drawn as a triangle mesh. Each atom is one quad, drawn with
instancing: the vertex shader receives the four corners of a unit square and
the atom's row of the sphere array, and places the quad in front of the atom,
facing the viewer and just large enough to cover it. The fragment shader then
computes, for each pixel, where the view ray meets the sphere: with an
orthographic camera, a pixel at offset :math:`(u, v)` from the centre, in
units of the radius, hits the front of the sphere at depth
:math:`\sqrt{1 - u^2 - v^2}` and misses it where
:math:`u^2 + v^2 > 1`. From the hit point the shader computes the normal for
shading and writes the exact depth.

Bonds work the same way. Each bond is a quad stretched over the projected
cylinder, and the fragment shader solves the ray-cylinder intersection, a
quadratic equation, for each pixel. Dashed bonds discard pixels in every
other segment along the axis.

Because every pixel gets its true depth, spheres and cylinders intersect
correctly where bonds enter atoms, with no gaps or tessellation artefacts at
any zoom. The cost per frame is one quad per atom and per bond line, and
rotating only changes two matrices, so a protein with tens of thousands of
atoms still draws at full frame rate.

Shading and outlines
====================

Shading is a soft wrap-around diffuse term with a small specular highlight,
chosen to read clearly at small sizes rather than to look photographic.
Distant parts can fade towards the background colour over the depth range of
the structure. Silhouettes get a constant-width dark outline: each impostor
knows its edge distance in pixels, so the outline is a simple blend near the
edge. Highlighted atoms get a coloured halo the same way.

Surfaces
========

The renderer keeps two meshes (:class:`~miles.render.gl.MeshSlot`):
the cartoon and the arrows that plugins draw (the dipole, vibrational
modes), and a list of further meshes from plugins, such as the molecular
surface and the isosurfaces of grids and orbitals
(:meth:`~miles.render.gl.GLRenderer.upload_layers`). The cartoon and
the arrows are opaque and get the same outline (back faces pushed out along
the normals); :func:`~miles.render.arrows.arrow_mesh` builds closed
arrows, a capped cylinder and a cone, for that reason. The QPainter
renderer draws the same arrows as a shaft among the bond cylinders and a
cone (a shaded triangle and its elliptical base), sorted by depth with the
atoms. The surface is drawn after everything opaque, and the plugins'
meshes after it; all follow the same rules for see-through drawing. Drawn opaque it is an ordinary mesh; drawn see-through,
only its front faces are drawn, blended over what is behind them with depth
writing turned off, so atoms and cartoons inside remain visible. Its
opacity is a uniform, so switching between opaque and see-through needs no
new upload: when the same meshes are shown with another opacity, the canvas
calls :meth:`~miles.render.gl.GLRenderer.set_layer_alphas` instead of
uploading them again. ``render/surface.py`` colours the surface: by element, by
hydrophobicity (blue polar to orange apolar), by charge (red negative to
blue positive), by the electrostatic potential at each vertex (same colours,
with the range stored in the mesh's ``value_range``) or in one neutral
colour.

The QPainter renderer draws surfaces and isosurfaces as well
(:mod:`miles.render.triangles`): each triangle that faces the viewer
becomes a polygon, flat-shaded with the OpenGL renderer's light, sorted by
depth together with the atoms and bonds and drawn without antialiasing so
neighbouring triangles leave no seams; see-through meshes are drawn with
their opacity. Back faces are skipped, since the meshes are closed. Drawing
costs a few microseconds per triangle, so meshes of more than
:data:`~miles.render.triangles.MAX_TRIANGLES` triangles are first
simplified by vertex clustering (vertices sharing a grid cell merge,
collapsed triangles go); the simplified mesh is cached per mesh. A protein
surface of about 80000 triangles then paints in about 75 ms, a small
molecule's surface or an orbital in about 10 ms.

Surfaces are often computed before the next paint, under a wait cursor,
which clears the scene cache's change flag; the OpenGL canvas therefore also
re-uploads the surface and isosurface meshes whenever the mesh object
differs from the one it uploaded last.

Uploads
=======

The renderer keeps one vertex buffer for spheres, one for cylinders, one for
the cartoon mesh and a small one with the highlight of each sphere, and the
canvas re-uploads only the ones whose contents changed (see :doc:`interface`).
Hovering an atom sends only the highlights that changed: the canvas keeps
the highlights it sent last and writes the differing rows, 4 bytes each
(:meth:`~miles.render.gl.GLRenderer.upload_flags`), or the whole highlight
buffer when many changed. The spheres, the cylinders and the mesh stay on
the GPU.

A window has one OpenGL widget, :class:`~miles.ui.canvas.surfaces.GLView`,
under its canvases; it draws the canvas shown. The buffers, vertex arrays
and shader programs belong to its context, so all tabs share one context,
one set of compiled shaders and one framebuffer the size of the window
(about 22 MB for 1000 by 700 points on a Retina screen), however many tabs
are open. Switching tabs sends the new canvas's arrays to the GPU once.
:meth:`~miles.render.gl.GLRenderer.release` deletes the GPU objects; the
view calls it when it is deleted or its context is about to be destroyed.

Pictures for export
===================

Every framebuffer, renderbuffer and texture that
:meth:`~miles.render.gl.GLRenderer.render_image` makes is deleted when it
returns, also when drawing fails, and each framebuffer is checked before it
is drawn into: a picture larger than the graphics card can hold raises an
error that says so instead of returning an empty image.

:meth:`~miles.render.gl.GLRenderer.render_image` draws into an
off-screen framebuffer with up to eight samples per pixel, at up to four
times the window resolution, then resolves it into colour and depth
textures. With ambient occlusion enabled, a full-screen pass darkens each
pixel by the fraction of 32 nearby samples, taken on a spiral within about
1.4 angstrom, whose depth lies in front of it. This screen-space ambient
occlusion gives creases and pockets soft shadows. The picture is read back
with premultiplied alpha, so transparent backgrounds composite correctly.

Picking
*******

:class:`~miles.render.picking.Picker` finds what is under the cursor
without the GPU. It projects all atoms to the screen with one array
operation and keeps those whose disc contains the cursor; among them the one
nearest the viewer wins. If no atom is hit, the distance from the cursor to
every projected bond segment is computed, again as arrays, and the front-most
bond within a few pixels wins. Picking returns an
:class:`~miles.render.picking.AtomHit`, a
:class:`~miles.render.picking.BondHit` or ``None``.

The picker computes the per-atom pick radii and the bond list once for a
structure and style, and reads the coordinates at each call, so moving atoms
does not require a new picker. In a cartoon only the CA (or P) atom of each
residue can be picked, with a radius that roughly covers the ribbon.

Cartoons
********

:func:`~miles.render.cartoon.cartoon_mesh` builds the ribbon of
proteins and nucleic acids.

1. Residues are grouped into unbroken traces: consecutive residues of the
   same chain whose guide atoms (CA, or P for nucleic acids) are close
   enough to be bonded through the backbone.
2. A Catmull-Rom spline through the guide atoms gives eight points per
   residue with their tangents.
3. A side vector at each residue points from CA towards the carbonyl oxygen,
   which lies in the plane of the peptide; successive side vectors are
   flipped where needed so the ribbon does not twist, and are interpolated
   along the spline.
4. At each spline point a superellipse cross-section is swept, with a size
   and shape that depend on the secondary structure: a thin round tube for
   coils, a wide flat ribbon for helices and a flat band for strands that
   widens into an arrow head before the strand ends. Transitions are
   smoothed over a residue.
5. The rings of points are joined into triangles and the ends capped.

Each vertex remembers its residue, so selecting atoms can tint the residues
they belong to by recolouring the vertices, without rebuilding the mesh.
The mesh is drawn with an outline pass: back faces pushed out along their
normals and drawn in the outline colour, then the surface itself.

The QPainter renderer
*********************

When OpenGL is not available, :func:`~miles.render.painter.paint_scene`
draws the same scene with QPainter. Spheres become discs with a radial
gradient and bonds become shaded quads, one per half bond, all sorted back to
front. In ball-and-stick mode the part of each bond inside its atoms is
trimmed, so bonds appear to enter the spheres. A cartoon is drawn as a thick
backbone trace following the spline. The same function draws the fragment
thumbnails in the fragment popover and the tool rail.

Themes
******

:class:`~miles.render.theme.Theme` holds the colours of the interface
and the molecule. The light and dark themes share an element palette softer
than the usual CPK colours; carbon and hydrogen differ between them so they
stand out against either background. Elements not in the palette get a
desaturated version of their standard colour.
