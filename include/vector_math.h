/*
 * Phyzzy: Playful Physics for Imaginative Games
 *
 * Licensed under the GNU General Public License, Version 3.
 * For license details, visit: https://www.gnu.org/licenses/gpl-3.0.html
 *
 * Questions or contributions? Reach out to Simon Devenish:
 * simon.devenish@outlook.com
 */

#pragma once

#include <math.h>
#include <stdbool.h>

#include "math_types.h"

#define FORCEINLINE static inline

FORCEINLINE f32 min(f32 a, f32 b) { return a < b ? a : b; }
FORCEINLINE f32 max(f32 a, f32 b) { return a > b ? a : b; }
FORCEINLINE f32 clamp(f32 value, f32 min_val, f32 max_val) {
    return max(min(value, max_val), min_val);
}

// v2 operations
FORCEINLINE v2 v2_add(const v2* a, const v2* b) {
    v2 result;
    result.x = a->x + b->x;
    result.y = a->y + b->y;
    return result;
}

FORCEINLINE v2 v2_subtract(const v2* a, const v2* b) {
    v2 result;
    result.x = a->x - b->x;
    result.y = a->y - b->y;
    return result;
}

FORCEINLINE v2 v2_scale(const v2* v, f32 s) {
    v2 result;
    result.x = v->x * s;
    result.y = v->y * s;
    return result;
}

FORCEINLINE v2 v2_divide_scalar(const v2* v, f32 s) {
    f32 inv_s = 1.0f / s;
    v2 result;
    result.x = v->x * inv_s;
    result.y = v->y * inv_s;
    return result;
}

FORCEINLINE v2 v2_negate(const v2* v) {
    v2 result;
    result.x = -v->x;
    result.y = -v->y;
    return result;
}

FORCEINLINE void v2_divide_scalar_inplace(v2* v, f32 s) {
    f32 inv_s = 1.0f / s;
    v->x *= inv_s;
    v->y *= inv_s;
}

FORCEINLINE bool v2_equals(const v2* a, const v2* b) {
    return a->x == b->x && a->y == b->y;
}

FORCEINLINE bool v2_not_equals(const v2* a, const v2* b) {
    return a->x != b->x || a->y != b->y;
}

FORCEINLINE f32 v2_magnitude(const v2* v) {
    return sqrtf((v->x * v->x) + (v->y * v->y));
}

FORCEINLINE v2 v2_normalize(const v2* v) {
    f32 mag = v2_magnitude(v);
    if (mag > 0.0f) {
        return v2_scale(v, 1.0f / mag);
    } else {
        return (v2){0.0f, 0.0f};
    }
}

FORCEINLINE v2 v2_abs(const v2* v) {
    v2 result;
    result.x = fabsf(v->x);
    result.y = fabsf(v->y);
    return result;
}

FORCEINLINE f32 v2_dot_product(const v2* a, const v2* b) {
    return (a->x * b->x) + (a->y * b->y);
}

// v3 operations
// v3 operations
FORCEINLINE v3 v3_add(v3 a, v3 b) {
    v3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

FORCEINLINE v3 v3_subtract(v3 a, v3 b) {
    v3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

FORCEINLINE v3 v3_scale(v3 v, f32 s) {
    v3 result;
    result.x = v.x * s;
    result.y = v.y * s;
    result.z = v.z * s;
    return result;
}

FORCEINLINE v3 v3_divide_scalar(const v3* v, f32 s) {
    f32 inv_s = 1.0f / s;
    v3 result;
    result.x = v->x * inv_s;
    result.y = v->y * inv_s;
    result.z = v->z * inv_s;
    return result;
}

FORCEINLINE v3 v3_negate(v3 v) {
    v3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;
    return result;
}

FORCEINLINE void v3_divide_scalar_inplace(v3* v, f32 s) {
    f32 inv_s = 1.0f / s;
    v->x *= inv_s;
    v->y *= inv_s;
    v->z *= inv_s;
}

FORCEINLINE bool v3_equals(const v3* a, const v3* b) {
    return a->x == b->x && a->y == b->y && a->z == b->z;
}

FORCEINLINE bool v3_not_equals(const v3* a, const v3* b) {
    return a->x != b->x || a->y != b->y || a->z != b->z;
}

FORCEINLINE f32 v3_magnitude(v3 v) {
    return sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

FORCEINLINE v3 v3_normalize(v3 v) {
    f32 mag = v3_magnitude(v);
    if (mag > 0.0f) {
        return v3_scale(v, 1.0f / mag);
    } else {
        return (v3){0.0f, 0.0f, 0.0f};
    }
}

FORCEINLINE v3 v3_abs(v3 v) {
    v3 result;
    result.x = fabsf(v.x);
    result.y = fabsf(v.y);
    result.z = fabsf(v.z);
    return result;
}

FORCEINLINE f32 v3_dot_product(v3 a, v3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

FORCEINLINE v3 v3_cross_product(v3 a, v3 b) {
    v3 result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

// v4 operations
FORCEINLINE v4 v4_add(const v4* a, const v4* b) {
    v4 result;
    result.x = a->x + b->x;
    result.y = a->y + b->y;
    result.z = a->z + b->z;
    result.w = a->w + b->w;
    return result;
}

FORCEINLINE v4 v4_subtract(const v4* a, const v4* b) {
    v4 result;
    result.x = a->x - b->x;
    result.y = a->y - b->y;
    result.z = a->z - b->z;
    result.w = a->w - b->w;
    return result;
}

FORCEINLINE v4 v4_scale(const v4* v, f32 s) {
    v4 result;
    result.x = v->x * s;
    result.y = v->y * s;
    result.z = v->z * s;
    result.w = v->w * s;
    return result;
}

FORCEINLINE v4 v4_divide_scalar(const v4* v, f32 s) {
    f32 inv_s = 1.0f / s;
    v4 result;
    result.x = v->x * inv_s;
    result.y = v->y * inv_s;
    result.z = v->z * inv_s;
    result.w = v->w * inv_s;
    return result;
}

FORCEINLINE v4 v4_negate(const v4* v) {
    v4 result;
    result.x = -v->x;
    result.y = -v->y;
    result.z = -v->z;
    result.w = -v->w;
    return result;
}

FORCEINLINE bool v4_equals(const v4* a, const v4* b) {
    return a->x == b->x && a->y == b->y && a->z == b->z && a->w == b->w;
}

FORCEINLINE bool v4_not_equals(const v4* a, const v4* b) {
    return a->x != b->x || a->y != b->y || a->z != b->z || a->w != b->w;
}

// Vector creation functions
FORCEINLINE v2 v2_make(f32 x, f32 y) {
    v2 vec;
    vec.x = x;
    vec.y = y;
    return vec;
}

FORCEINLINE v3 v3_make(f32 x, f32 y, f32 z) {
    v3 vec;
    vec.x = x;
    vec.y = y;
    vec.z = z;
    return vec;
}

FORCEINLINE v4 v4_make(f32 x, f32 y, f32 z, f32 w) {
    v4 vec;
    vec.x = x;
    vec.y = y;
    vec.z = z;
    vec.w = w;
    return vec;
}

FORCEINLINE v3 m3_multiply_v3(m3 matrix, v3 vector) {
    v3 result;
    result.x = matrix.m[0][0] * vector.x + matrix.m[0][1] * vector.y + matrix.m[0][2] * vector.z;
    result.y = matrix.m[1][0] * vector.x + matrix.m[1][1] * vector.y + matrix.m[1][2] * vector.z;
    result.z = matrix.m[2][0] * vector.x + matrix.m[2][1] * vector.y + matrix.m[2][2] * vector.z;
    return result;
}


// Matrix and vector multiplication
FORCEINLINE v2 m4_multiply_v2(const m4* m, const v2* v) {
    v2 result;
    result.x = m->table[0][0] * v->x + m->table[0][1] * v->y + m->table[0][2];
    result.y = m->table[1][0] * v->x + m->table[1][1] * v->y + m->table[1][2];
    return result;
}

FORCEINLINE v3 m4_multiply_v3(const m4* m, const v3* v) {
    v3 result;
    result.x = m->table[0][0] * v->x + m->table[0][1] * v->y + m->table[0][2] * v->z + m->table[0][3];
    result.y = m->table[1][0] * v->x + m->table[1][1] * v->y + m->table[1][2] * v->z + m->table[1][3];
    result.z = m->table[2][0] * v->x + m->table[2][1] * v->y + m->table[2][2] * v->z + m->table[2][3];
    return result;
}

FORCEINLINE v4 m4_multiply_v4(const m4* m, const v4* v) {
    v4 result;
    result.x = m->table[0][0] * v->x + m->table[0][1] * v->y + m->table[0][2] * v->z + m->table[0][3] * v->w;
    result.y = m->table[1][0] * v->x + m->table[1][1] * v->y + m->table[1][2] * v->z + m->table[1][3] * v->w;
    result.z = m->table[2][0] * v->x + m->table[2][1] * v->y + m->table[2][2] * v->z + m->table[2][3] * v->w;
    result.w = m->table[3][0] * v->x + m->table[3][1] * v->y + m->table[3][2] * v->z + m->table[3][3] * v->w;
    return result;
}

// Projection, rejection, reflection, angle between vectors
FORCEINLINE v3 v3_project(v3 a, v3 b) {
    f32 scalar = v3_dot_product(a, b) / v3_dot_product(b, b);
    return v3_scale(b, scalar);
}

FORCEINLINE v3 v3_reject(v3 a, v3 b) {
    v3 projection = v3_project(a, b);
    return v3_subtract(a, projection);
}

FORCEINLINE v3 v3_reflect(v3 incident, v3 normal) {
    f32 scalar = 2.0f * v3_dot_product(incident, normal);
    v3 scaled_normal = v3_scale(normal, scalar);
    return v3_subtract(incident, scaled_normal);
}

FORCEINLINE f32 v3_angle_between(v3 a, v3 b) {
    f32 dot = v3_dot_product(a, b);
    f32 mag_a = sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
    f32 mag_b = sqrtf(b.x * b.x + b.y * b.y + b.z * b.z);
    return acosf(clamp(dot / (mag_a * mag_b), -1.0f, 1.0f));
}

// Matrix multiplication
FORCEINLINE m4 m4_multiply(const m4* a, const m4* b) {
    m4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.table[i][j] = a->table[i][0] * b->table[0][j] +
                                 a->table[i][1] * b->table[1][j] +
                                 a->table[i][2] * b->table[2][j] +
                                 a->table[i][3] * b->table[3][j];
        }
    }
    return result;
}

// Identity matrix
static const m4 identity = {{
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
}};

// Matrix determinant (placeholder)
FORCEINLINE f32 m4_determinant(const m4* m) {
    // Placeholder implementation
    return 0.0f;
}

// Matrix scaling
FORCEINLINE m4 m4_scale(const m4* m, const v3* s) {
    m4 scale = identity;
    scale.table[0][0] = s->x;
    scale.table[1][1] = s->y;
    scale.table[2][2] = s->z;
    return m4_multiply(&scale, m);
}

// Matrix translation
FORCEINLINE m4 m4_translate(const v3* t) {
    m4 trans = identity;
    trans.table[0][3] = t->x;
    trans.table[1][3] = t->y;
    trans.table[2][3] = t->z;
    return trans;
}

// Matrix rotations
FORCEINLINE m4 m4_x_rotation(f32 angle) {
    f32 c = cosf(angle);
    f32 s = sinf(angle);
    m4 rot = identity;
    rot.table[1][1] = c;
    rot.table[1][2] = -s;
    rot.table[2][1] = s;
    rot.table[2][2] = c;
    return rot;
}

FORCEINLINE m4 m4_y_rotation(f32 angle) {
    f32 c = cosf(angle);
    f32 s = sinf(angle);
    m4 rot = identity;
    rot.table[0][0] = c;
    rot.table[0][2] = s;
    rot.table[2][0] = -s;
    rot.table[2][2] = c;
    return rot;
}

FORCEINLINE m4 m4_z_rotation(f32 angle) {
    f32 c = cosf(angle);
    f32 s = sinf(angle);
    m4 rot = identity;
    rot.table[0][0] = c;
    rot.table[0][1] = -s;
    rot.table[1][0] = s;
    rot.table[1][1] = c;
    return rot;
}

// Matrix transpose
FORCEINLINE m4 m4_transpose(const m4* m) {
    m4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.table[i][j] = m->table[j][i];
        }
    }
    return result;
}

// Matrix inversion (placeholder)
FORCEINLINE bool m4_invert(const m4* m, m4* out) {
    // Placeholder: Actual implementation is complex
    return false;
}

// Frustum creation
FORCEINLINE m4 frustum(f32 l, f32 r, f32 b, f32 t, f32 n, f32 f) {
    m4 mat = {0};
    mat.table[0][0] = 2.0f * n / (r - l);
    mat.table[1][1] = 2.0f * n / (t - b);
    mat.table[0][2] = (r + l) / (r - l);
    mat.table[1][2] = (t + b) / (t - b);
    mat.table[2][2] = -(f + n) / (f - n);
    mat.table[2][3] = -2.0f * f * n / (f - n);
    mat.table[3][2] = -1.0f;
    return mat;
}

// Point inside frustum
FORCEINLINE bool point_inside_frustum(const Frustum* frustum, const v3* point) {
    for (int i = 0; i < 6; ++i) {
        const Plane* plane = &frustum->planes[i];
        if (v3_dot_product(plane->normal, *point) + plane->d < 0) {
            return false;
        }
    }
    return true;
}

// AABB inside frustum
FORCEINLINE bool aabb_inside_frustum(const Frustum* frustum, const AABB* aabb) {
    for (int i = 0; i < 6; ++i) {
        const Plane* plane = &frustum->planes[i];
        v3 p = aabb->max;
        if (plane->normal.x < 0) p.x = aabb->min.x;
        if (plane->normal.y < 0) p.y = aabb->min.y;
        if (plane->normal.z < 0) p.z = aabb->min.z;
        if (v3_dot_product(plane->normal, p) + plane->d < 0) {
            return false;
        }
    }
    return true;
}

FORCEINLINE m4 look_at(v3 eye, v3 target, v3 up) {
    v3 forward = v3_normalize(v3_subtract(target, eye));
    v3 side = v3_normalize(v3_cross_product(forward, up));
    v3 up_corrected = v3_cross_product(side, forward);

    m4 view = identity;
    view.table[0][0] = side.x;
    view.table[1][0] = side.y;
    view.table[2][0] = side.z;
    view.table[0][1] = up_corrected.x;
    view.table[1][1] = up_corrected.y;
    view.table[2][1] = up_corrected.z;
    view.table[0][2] = -forward.x;
    view.table[1][2] = -forward.y;
    view.table[2][2] = -forward.z;
    view.table[0][3] = -v3_dot_product(side, eye);
    view.table[1][3] = -v3_dot_product(up_corrected, eye);
    view.table[2][3] = v3_dot_product(forward, eye);
    return view;
}

// Orthographic projection matrix
FORCEINLINE m4 orthographic(f32 left, f32 right, f32 bottom, f32 top, f32 near, f32 far) {
    m4 ortho = identity;
    ortho.table[0][0] = 2.0f / (right - left);
    ortho.table[1][1] = 2.0f / (top - bottom);
    ortho.table[2][2] = -2.0f / (far - near);
    ortho.table[0][3] = -(right + left) / (right - left);
    ortho.table[1][3] = -(top + bottom) / (top - bottom);
    ortho.table[2][3] = -(far + near) / (far - near);
    return ortho;
}

// Perspective projection matrix
FORCEINLINE m4 perspective(f32 fov, f32 aspect, f32 near_plane, f32 far_plane) {
    f32 tangent = tanf(fov * 0.5f * (3.14159265f / 180.0f));
    f32 height = near_plane * tangent;
    f32 width = height * aspect;
    return frustum(-width, width, -height, height, near_plane, far_plane);
}

// Quaternion from axis-angle
FORCEINLINE quat quat_from_axis_angle(const v3* axis, f32 angle) {
    f32 half_angle = angle * 0.5f;
    f32 s = sinf(half_angle);
    quat q;
    q.x = axis->x * s;
    q.y = axis->y * s;
    q.z = axis->z * s;
    q.w = cosf(half_angle);
    return q;
}

// Quaternion to matrix
FORCEINLINE m4 quat_to_matrix(const quat* q) {
    m4 result = identity;
    f32 xx = q->x * q->x;
    f32 yy = q->y * q->y;
    f32 zz = q->z * q->z;
    f32 xy = q->x * q->y;
    f32 xz = q->x * q->z;
    f32 yz = q->y * q->z;
    f32 wx = q->w * q->x;
    f32 wy = q->w * q->y;
    f32 wz = q->w * q->z;

    result.table[0][0] = 1.0f - 2.0f * (yy + zz);
    result.table[0][1] = 2.0f * (xy - wz);
    result.table[0][2] = 2.0f * (xz + wy);

    result.table[1][0] = 2.0f * (xy + wz);
    result.table[1][1] = 1.0f - 2.0f * (xx + zz);
    result.table[1][2] = 2.0f * (yz - wx);

    result.table[2][0] = 2.0f * (xz - wy);
    result.table[2][1] = 2.0f * (yz + wx);
    result.table[2][2] = 1.0f - 2.0f * (xx + yy);

    return result;
}

// Quaternion multiplication
FORCEINLINE quat quat_multiply(quat q1, quat q2) {
    quat result;
    result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return result;
}


// Quaternion normalization
FORCEINLINE quat quat_normalize(const quat* q) {
    f32 mag = sqrtf(q->x * q->x + q->y * q->y + q->z * q->z + q->w * q->w);
    if (mag > 0.0f) {
        f32 inv_mag = 1.0f / mag;
        quat result;
        result.x = q->x * inv_mag;
        result.y = q->y * inv_mag;
        result.z = q->z * inv_mag;
        result.w = q->w * inv_mag;
        return result;
    } else {
        return (quat){0.0f, 0.0f, 0.0f, 1.0f};
    }
}

// Quaternion slerp
FORCEINLINE quat quat_slerp(const quat* q1, const quat* q2, f32 t) {
    f32 dot = q1->x * q2->x + q1->y * q2->y + q1->z * q2->z + q1->w * q2->w;

    quat q2_copy;
    if (dot < 0.0f) {
        dot = -dot;
        q2_copy.x = -q2->x;
        q2_copy.y = -q2->y;
        q2_copy.z = -q2->z;
        q2_copy.w = -q2->w;
    } else {
        q2_copy = *q2;
    }

    const f32 DOT_THRESHOLD = 0.9995f;
    if (dot > DOT_THRESHOLD) {
        quat result;
        result.x = q1->x + t * (q2_copy.x - q1->x);
        result.y = q1->y + t * (q2_copy.y - q1->y);
        result.z = q1->z + t * (q2_copy.z - q1->z);
        result.w = q1->w + t * (q2_copy.w - q1->w);
        return quat_normalize(&result);
    }

    f32 theta_0 = acosf(dot);
    f32 theta = theta_0 * t;

    quat q3;
    q3.x = q2_copy.x - q1->x * dot;
    q3.y = q2_copy.y - q1->y * dot;
    q3.z = q2_copy.z - q1->z * dot;
    q3.w = q2_copy.w - q1->w * dot;
    q3 = quat_normalize(&q3);

    f32 s0 = cosf(theta);
    f32 s1 = sinf(theta);

    quat result;
    result.x = q1->x * s0 + q3.x * s1;
    result.y = q1->y * s0 + q3.y * s1;
    result.z = q1->z * s0 + q3.z * s1;
    result.w = q1->w * s0 + q3.w * s1;
    return result;
}

// Quaternion inverse
FORCEINLINE quat quat_inverse(const quat* q) {
    f32 norm_sq = q->x * q->x + q->y * q->y + q->z * q->z + q->w * q->w;
    if (norm_sq > 0.0f) {
        f32 inv_norm_sq = 1.0f / norm_sq;
        quat result;
        result.x = -q->x * inv_norm_sq;
        result.y = -q->y * inv_norm_sq;
        result.z = -q->z * inv_norm_sq;
        result.w = q->w * inv_norm_sq;
        return result;
    } else {
        return (quat){0.0f, 0.0f, 0.0f, 1.0f};
    }
}

// Quaternion dot product
FORCEINLINE f32 quat_dot(const quat* q1, const quat* q2) {
    return q1->x * q2->x + q1->y * q2->y + q1->z * q2->z + q1->w * q2->w;
}

// Quaternion to axis-angle
FORCEINLINE void quat_to_axis_angle(const quat* q, v3* axis, f32* angle) {
    f32 qw = q->w;
    if (qw > 1.0f) {
        quat q_normalized = quat_normalize(q);
        qw = q_normalized.w;
    }
    *angle = 2.0f * acosf(qw);
    f32 s = sqrtf(1.0f - qw * qw);
    if (s < 0.001f) {
        axis->x = q->x;
        axis->y = q->y;
        axis->z = q->z;
    } else {
        axis->x = q->x / s;
        axis->y = q->y / s;
        axis->z = q->z / s;
    }
}

// AABB intersection
FORCEINLINE bool aabb_intersects(const AABB* a, const AABB* b) {
    return (a->min.x <= b->max.x && a->max.x >= b->min.x) &&
           (a->min.y <= b->max.y && a->max.y >= b->min.y) &&
           (a->min.z <= b->max.z && a->max.z >= b->min.z);
}

// Ray-AABB intersection
FORCEINLINE bool ray_aabb_intersect(const Ray* ray, const AABB* aabb, f32* t) {
    f32 tmin = (aabb->min.x - ray->origin.x) / ray->direction.x;
    f32 tmax = (aabb->max.x - ray->origin.x) / ray->direction.x;

    if (tmin > tmax) { f32 temp = tmin; tmin = tmax; tmax = temp; }

    f32 tymin = (aabb->min.y - ray->origin.y) / ray->direction.y;
    f32 tymax = (aabb->max.y - ray->origin.y) / ray->direction.y;

    if (tymin > tymax) { f32 temp = tymin; tymin = tymax; tymax = temp; }

    if ((tmin > tymax) || (tymin > tmax)) return false;

    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    f32 tzmin = (aabb->min.z - ray->origin.z) / ray->direction.z;
    f32 tzmax = (aabb->max.z - ray->origin.z) / ray->direction.z;

    if (tzmin > tzmax) { f32 temp = tzmin; tzmin = tzmax; tzmax = temp; }

    if ((tmin > tzmax) || (tzmin > tmax)) return false;

    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    *t = tmin;
    return true;
}

// Ray-Sphere intersection
FORCEINLINE bool ray_sphere_intersect(const Ray* ray, const Sphere* sphere, f32* t) {
    v3 oc = v3_subtract(ray->origin, sphere->center);
    f32 a = v3_dot_product(ray->direction, ray->direction);
    f32 b = 2.0f * v3_dot_product(oc, ray->direction);
    f32 c = v3_dot_product(oc, oc) - sphere->radius * sphere->radius;
    f32 discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    *t = (-b - sqrtf(discriminant)) / (2.0f * a);
    return true;
}

// Linear interpolation
FORCEINLINE f32 lerp(f32 a, f32 b, f32 t) {
    return a * (1.0f - t) + b * t;
}

// Vector3 linear interpolation
FORCEINLINE v3 lerp_v3(const v3* a, const v3* b, f32 t) {
    v3 result;
    result.x = lerp(a->x, b->x, t);
    result.y = lerp(a->y, b->y, t);
    result.z = lerp(a->z, b->z, t);
    return result;
}

// Cubic interpolation
FORCEINLINE f32 cubic_interpolate(f32 p0, f32 p1, f32 p2, f32 p3, f32 t) {
    f32 t2 = t * t;
    f32 t3 = t2 * t;
    return 0.5f * ((2.0f * p1) +
                   (-p0 + p2) * t +
                   (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                   (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}