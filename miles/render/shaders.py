"""GLSL sources for the OpenGL renderer.

Spheres and cylinders are drawn as screen-aligned quads (impostors): the
fragment shader intersects the view ray with the exact surface, shades it and
writes the true depth, so intersections between atoms and bonds are exact.

Uniforms shared by the molecule shaders:

``u_px``
    Size of one framebuffer pixel in angstrom.
``u_outline``, ``u_outline_color``
    Silhouette width in pixels and its colour.
``u_bg``, ``u_fade``, ``u_zrange``
    Depth fading towards the background colour over the scene's depth range.
``u_flag_colors``
    Highlight colours, indexed by ``scene.Highlight``.
"""

_COMMON = """
uniform float u_px;
uniform float u_outline;
uniform vec3 u_outline_color;
uniform vec3 u_bg;
uniform float u_fade;
uniform vec2 u_zrange;
uniform vec3 u_flag_colors[5];

vec3 shade(vec3 base, vec3 n) {
    vec3 L = normalize(vec3(-0.38, 0.52, 0.76));
    float wrap = clamp(dot(n, L) * 0.5 + 0.5, 0.0, 1.0);
    vec3 col = base * (0.46 + 0.58 * wrap * wrap);
    vec3 H = normalize(L + vec3(0.0, 0.0, 1.0));
    col += vec3(0.14) * pow(max(dot(n, H), 0.0), 38.0);
    col *= mix(0.84, 1.0, smoothstep(0.0, 0.55, n.z));
    return col;
}

vec3 fade(vec3 col, float z) {
    float span = max(u_zrange.y - u_zrange.x, 8.0);
    float t = clamp((u_zrange.y - z) / span, 0.0, 1.0);
    return mix(col, u_bg, t * u_fade);
}
"""

SPHERE_VS = """
#version 330 core
layout(location = 0) in vec2 a_corner;
layout(location = 1) in vec3 a_center;
layout(location = 2) in float a_radius;
layout(location = 3) in vec3 a_color;
layout(location = 4) in float a_flag;
uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_px;
uniform float u_halo;
uniform float u_dpr;
out vec2 v_uv;
out vec3 v_center;
out float v_radius;
out vec3 v_color;
flat out int v_flag;
flat out float v_fixed;
void main() {
    v_fixed = a_radius < 0.0 ? 1.0 : 0.0;
    float r = a_radius < 0.0 ? -a_radius * u_px * u_dpr : a_radius;
    vec4 c = u_view * vec4(a_center, 1.0);
    float margin = ((a_flag > 0.5 ? u_halo : 0.0) + 1.5) * u_px;
    v_uv = a_corner * (r + margin) / r;
    v_center = c.xyz;
    v_radius = r;
    v_color = a_color;
    v_flag = int(a_flag + 0.5);
    gl_Position = u_proj * vec4(c.xy + v_uv * r, c.z + r, 1.0);
}
"""

SPHERE_FS = """
#version 330 core
in vec2 v_uv;
in vec3 v_center;
in float v_radius;
in vec3 v_color;
flat in int v_flag;
flat in float v_fixed;
uniform mat4 u_proj;
uniform float u_halo;
out vec4 o_color;
""" + _COMMON + """
void main() {
    float d = length(v_uv);
    float rpx = v_radius / u_px;
    float edge = (1.0 - d) * rpx;
    float halo = v_flag > 0 ? u_halo : 0.0;
    if (edge < -halo - 1.0) discard;
    vec3 col;
    float alpha;
    float z;
    float ow = v_fixed > 0.5 ? 0.0 : min(u_outline, rpx * 0.3);
    if (edge < 0.0) {
        col = u_flag_colors[v_flag];
        alpha = clamp(edge + halo + 0.5, 0.0, 1.0);
        z = v_center.z;
    } else {
        vec3 n = vec3(v_uv, sqrt(max(1.0 - d * d, 0.0)));
        z = v_center.z + n.z * v_radius;
        col = fade(shade(v_color, n), z);
        vec3 oc = u_outline_color;
        if (v_flag > 0) {
            col = mix(col, u_flag_colors[v_flag], 0.18);
            oc = u_flag_colors[v_flag];
            ow = max(ow, 1.0);
        }
        col = mix(oc, col, clamp(edge - ow + 0.5, 0.0, 1.0));
        alpha = v_flag > 0 ? 1.0 : clamp(edge + 0.5, 0.0, 1.0);
    }
    if (alpha <= 0.0) discard;
    vec4 clip = u_proj * vec4(v_center.xy + v_uv * v_radius, z, 1.0);
    gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
    o_color = vec4(col, alpha);
}
"""

CYLINDER_VS = """
#version 330 core
layout(location = 0) in vec2 a_corner;
layout(location = 1) in vec3 a_p0;
layout(location = 2) in vec3 a_p1;
layout(location = 3) in float a_radius;
layout(location = 4) in vec3 a_c0;
layout(location = 5) in vec3 a_c1;
layout(location = 6) in float a_dash;
uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_px;
uniform float u_dpr;
out vec2 v_xy;
flat out float v_fixed;
flat out vec3 v_a;
flat out vec3 v_b;
flat out float v_r;
flat out vec3 v_c0;
flat out vec3 v_c1;
flat out float v_dash;
void main() {
    v_fixed = a_radius < 0.0 ? 1.0 : 0.0;
    float r = a_radius < 0.0 ? -a_radius * u_px * u_dpr : a_radius;
    vec3 a = (u_view * vec4(a_p0, 1.0)).xyz;
    vec3 b = (u_view * vec4(a_p1, 1.0)).xyz;
    vec2 d = b.xy - a.xy;
    float L = length(d);
    vec2 dir = L > 1e-6 ? d / L : vec2(1.0, 0.0);
    vec2 perp = vec2(-dir.y, dir.x);
    float m = r + 1.5 * u_px;
    vec2 base = a_corner.x < 0.0 ? a.xy - dir * m : b.xy + dir * m;
    v_xy = base + perp * a_corner.y * m;
    v_a = a;
    v_b = b;
    v_r = r;
    v_c0 = a_c0;
    v_c1 = a_c1;
    v_dash = a_dash;
    float z = (a_corner.x < 0.0 ? a.z : b.z) + r;
    gl_Position = u_proj * vec4(v_xy, z, 1.0);
}
"""

CYLINDER_FS = """
#version 330 core
in vec2 v_xy;
flat in float v_fixed;
flat in vec3 v_a;
flat in vec3 v_b;
flat in float v_r;
flat in vec3 v_c0;
flat in vec3 v_c1;
flat in float v_dash;
uniform mat4 u_proj;
out vec4 o_color;
""" + _COMMON + """
void main() {
    vec3 axis = v_b - v_a;
    float L = length(axis);
    if (L < 1e-6) discard;
    vec3 ax = axis / L;
    vec2 e2 = v_b.xy - v_a.xy;
    float l2 = length(e2);
    vec2 e = l2 > 1e-6 ? e2 / l2 : vec2(1.0, 0.0);
    vec2 rel = v_xy - v_a.xy;
    float dist = abs(e.x * rel.y - e.y * rel.x);
    float edge = (v_r - dist) / u_px;
    if (edge < -0.5) discard;
    vec3 w = vec3(rel, -v_a.z);
    vec3 Dp = vec3(0.0, 0.0, 1.0) - ax.z * ax;
    vec3 wp = w - dot(w, ax) * ax;
    float A = dot(Dp, Dp);
    if (A < 1e-7) discard;
    float B = 2.0 * dot(Dp, wp);
    float C = dot(wp, wp) - v_r * v_r;
    float disc = max(B * B - 4.0 * A * C, 0.0);
    float t = (-B + sqrt(disc)) / (2.0 * A);
    vec3 P = vec3(v_xy, t);
    float s = dot(P - v_a, ax);
    if (s < 0.0 || s > L) discard;
    if (v_dash > 0.5 && fract(s / 0.24) > 0.52) discard;
    vec3 n = normalize((P - v_a) - s * ax);
    vec3 base = s < 0.5 * L ? v_c0 : v_c1;
    vec3 col = fade(shade(base, n), P.z);
    float ow = (v_fixed > 0.5 || v_dash > 0.5) ? 0.0 : min(u_outline, v_r / u_px * 0.35);
    col = mix(u_outline_color, col, clamp(edge - ow + 0.5, 0.0, 1.0));
    float alpha = clamp(edge + 0.5, 0.0, 1.0);
    vec4 clip = u_proj * vec4(P, 1.0);
    gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
    o_color = vec4(col, alpha);
}
"""

AO_VS = """
#version 330 core
layout(location = 0) in vec2 a_corner;
out vec2 v_uv;
void main() {
    v_uv = a_corner * 0.5 + 0.5;
    gl_Position = vec4(a_corner, 0.0, 1.0);
}
"""

AO_FS = """
#version 330 core
in vec2 v_uv;
uniform sampler2D u_color;
uniform sampler2D u_depth;
uniform vec2 u_texel;
uniform float u_radius_px;
uniform float u_depth_range;
uniform float u_world_radius;
uniform float u_strength;
out vec4 o_color;
float viewz(vec2 uv) {
    return -(2.0 * texture(u_depth, uv).r - 1.0) * u_depth_range;
}
float hash(vec2 p) {
    return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453);
}
void main() {
    vec4 c = texture(u_color, v_uv);
    if (texture(u_depth, v_uv).r >= 0.99999) { o_color = c; return; }
    float z0 = viewz(v_uv);
    float occ = 0.0;
    const int N = 32;
    float a0 = hash(gl_FragCoord.xy) * 6.2831853;
    for (int i = 0; i < N; i++) {
        float t = (float(i) + 0.5) / float(N);
        float a = a0 + float(i) * 2.3999632;
        vec2 uv = v_uv + vec2(cos(a), sin(a)) * sqrt(t) * u_radius_px * u_texel;
        if (texture(u_depth, uv).r >= 0.99999) continue;
        float dz = viewz(uv) - z0;
        float range = 1.0 - smoothstep(u_world_radius, 2.5 * u_world_radius, dz);
        occ += smoothstep(0.03, 0.45 * u_world_radius, dz) * range;
    }
    float ao = 1.0 - u_strength * occ / float(N);
    o_color = vec4(c.rgb * ao, c.a);
}
"""

MESH_VS = """
#version 330 core
layout(location = 0) in vec3 a_pos;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec3 a_color;
uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_offset;
out vec3 v_n;
out vec3 v_color;
out float v_z;
void main() {
    vec4 p = u_view * vec4(a_pos + a_normal * u_offset, 1.0);
    v_n = mat3(u_view) * a_normal;
    v_color = a_color;
    v_z = p.z;
    gl_Position = u_proj * p;
}
"""

MESH_FS = """
#version 330 core
in vec3 v_n;
in vec3 v_color;
in float v_z;
uniform float u_outline_pass;
uniform float u_alpha;
out vec4 o_color;
""" + _COMMON + """
void main() {
    if (u_outline_pass > 0.5) {
        o_color = vec4(u_outline_color, 1.0);
        return;
    }
    vec3 n = normalize(v_n);
    if (!gl_FrontFacing) n = -n;
    o_color = vec4(fade(shade(v_color, n), v_z), u_alpha);
}
"""
