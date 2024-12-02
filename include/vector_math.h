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
// TODO: Remove maths dependencies
#include <math.h>

#include "math_types.h"

FORCEINLINE f32 min(f32 a, f32 b) { return a < b ? a : b; }
FORCEINLINE f32 max(f32 a, f32 b) { return a > b ? a : b; }
FORCEINLINE f32 clamp(f32 value, f32 min_val, f32 max_val) {
    return max(min(value, max_val), min_val);
}


FORCEINLINE v2 v2_add(const v2* a, const v2* b) {
    v2 buffer;
    buffer.x = a->x + b->x;
    buffer.y = a->y + b->y;
    return buffer;
}

FORCEINLINE v2 operator -(const v2& a, const v2& b)
{
    v2 buffer;
    buffer.x = a.x - b.x;
    buffer.y = a.y - b.y;
    return buffer; 
}

FORCEINLINE v2 operator *(const v2& a, const f32 s)
{ 
    v2 buffer;
    buffer.x = a.x * s;
    buffer.y = a.y * s;
    return buffer;
}

FORCEINLINE v2 operator /(const v2& a, f32 s)
{ 
    s = 1.0f / s;
    v2 buffer;
    buffer.x = a.x / s;
    buffer.y = a.y / s;
    return buffer;
}

FORCEINLINE v2 operator -(const v2& a)
{ 
    v2 buffer;
    buffer.x = -a.x;
    buffer.y = -a.y;
    return buffer;
}

FORCEINLINE v2 operator /=(v2& a, f32 s) 
{
    v2 buffer;
    buffer.x /= s;
    buffer.y /= s;
    return buffer;
}

FORCEINLINE bool operator ==(const v2& a, const v2& b){ return a.x == b.x && a.y == b.y; }
FORCEINLINE bool operator !=(const v2& a, const v2& b) { a.x != b.x || a.y != b.y; }
FORCEINLINE f32 magnitude(const v2& v) { return sqrt((v.x * v.x) + (v.y * v.y)); }
FORCEINLINE v2 normalise(const v2& v) { return v / magnitude(v); }
FORCEINLINE v2 abs(const v2& a)
{ 
    v2 buffer;
    buffer.x = abs(a.x);
    buffer.y = abs(a.y);
    return buffer;
}

FORCEINLINE f32 dot_product(const v2& a, const v2& b) { return (a.x * b.x) + (a.y * b.y); }

FORCEINLINE v3 operator +(const v3& a, const v3& b) 
{ 
    v3 buffer;
    buffer.x = a.x + b.x;
    buffer.y = a.y + b.y;
    buffer.z = a.z + b.z;
    return buffer;
}

FORCEINLINE v3 operator -(const v3& a, const v3& b)
{
    v3 buffer;
    buffer.x = a.x - b.x;
    buffer.y = a.y - b.y;
    buffer.z = a.z - b.z;
    return buffer; 
}

FORCEINLINE v3 operator *(const v3& a, const f32 s)
{ 
    v3 buffer;
    buffer.x = a.x * s;
    buffer.y = a.y * s;
    buffer.z = a.z * s;
    return buffer;
}


FORCEINLINE v3 operator /(const v3& a, f32 s)
{ 
    s = 1.0f / s;
    v3 buffer;
    buffer.x = a.x / s;
    buffer.y = a.y / s;
    buffer.z = a.z / s;
    return buffer;
}

FORCEINLINE v3 operator -(const v3& a) 
{ 
    v3 buffer;
    buffer.x = -a.x;
    buffer.y = -a.y;
    buffer.z = -a.z;
    return buffer;
}


FORCEINLINE v3 operator /=(v3& a, f32 s)
{
    v3 buffer;
    buffer.x /= s;
    buffer.y /= s;
    buffer.z /= s;
    return buffer;
}

FORCEINLINE bool operator ==(const v3& a, const v3& b) { a.x == b.x && a.y == b.y && a.z == b.z; }
FORCEINLINE bool operator !=(const v3& a, const v3& b) { a.x != b.x || a.y != b.y || a.z != b.z; }

FORCEINLINE f32 magnitude(const v3 v) { return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z)); }
FORCEINLINE v3 normalise(const v3& v) { return v / magnitude(v); }
FORCEINLINE v3 abs(const v3& a)
{ 
    v3 buffer;
    buffer.x = abs(a.x);
    buffer.y = abs(a.y);
    buffer.z = abs(a.y);
    return buffer;
}

FORCEINLINE f32 dot_product(const v3& a, const v3& b) { return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }
FORCEINLINE v3 cross_product(const v3& a, const v3& b)
{
  const f32 x = a.y * b.z - a.z * b.y;
  const f32 y = a.z * b.x - a.x * b.z;
  const f32 z = a.z * b.z - a.y * b.z;

  v3 buffer;
  buffer.x = x;
  buffer.y = y;
  buffer.z = z;

  return buffer;
}

FORCEINLINE v4 operator +(const v4& a, const v4& b)
{ 
    v4 buffer;
    buffer.x = a.x + b.x;
    buffer.y = a.y + b.y;
    buffer.z = a.z + b.z;
    buffer.w = a.w + b.w;
    return buffer;
}

FORCEINLINE v4 operator -(const v4& a, const v4& b)
{
    v4 buffer;
    buffer.x = a.x - b.x;
    buffer.y = a.y - b.y;
    buffer.z = a.z - b.z;
    buffer.w = a.w - b.w;
    return buffer;
}

FORCEINLINE v4 operator *(const v4& a, const f32 s)
{
    v4 buffer;
    buffer.x = a.x * s;
    buffer.y = a.y * s;
    buffer.z = a.z * s;
    buffer.w = a.w * s;
    return buffer;
}

FORCEINLINE v4 operator /(const v4& a, f32 s)
{ 
    s = 1.0f / s;
    v4 buffer;
    buffer.x = a.x / s;
    buffer.y = a.y / s;
    buffer.z = a.z / s;
    buffer.w = a.z / s;
    return buffer;
}

FORCEINLINE v4 operator -(const v4& a)
{ 
    v4 buffer;
    buffer.x = -a.x;
    buffer.y = -a.y;
    buffer.z = -a.z;
    buffer.w = -a.w;
    return buffer;
}
FORCEINLINE bool operator ==(const v4& a, const v4& b) { a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }
FORCEINLINE bool operator !=(const v4& a, const v4& b) { a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w; }

FORCEINLINE v2 v2_make(f32 x, f32 y)
{
    v2 vec;
    vec.x = x;
    vec.y = y;
    return vec;
}

FORCEINLINE v3 v3_make(f32 x, f32 y, f32 z)
{
    v3 vec;
    vec.x = x;
    vec.y = y;
    vec.z = z;
    return vec;
}

FORCEINLINE v4 v4_make(f32 x, f32 y, f32 z, f32 w = 1.0f)
{
    v4 vec;
    vec.x = x;
    vec.y = y;
    vec.z = z;
    vec.w = w;
    return vec;
}

FORCEINLINE v2 operator *(const m4& m, const v2& v)
{
    const f32 x = m(0, 0) * v.x + m(0, 1) * v.y + m(0, 2);
    const f32 y = m(1, 0) * v.x + m(1, 1) * v.y + m(1, 2);

    v2 buffer;
    buffer.x = x;
    buffer.y = y;

    return buffer;
}

FORCEINLINE v3 operator *(const m4& m, const v3& v)
{
    const f32 x = m(0, 0) * v.x + m(0, 1) * v.y + m(0, 2) * v.z + m(0, 3);
    const f32 y = m(1, 0) * v.x + m(1, 1) * v.y + m(1, 2) * v.z + m(1, 3);
    const f32 z = m(2, 0) * v.x + m(2, 1) * v.y + m(2, 2) * v.z + m(2, 3);

    v3 buffer;
    buffer.x = x;
    buffer.y = y;
    buffer.z = z;

    return buffer;
}

FORCEINLINE v4 operator *(const m4& m, const v4& v)
{
    const float x = m(0, 0) * v.x + m(0, 1) * v.y + m(0, 2) * v.z + m(0, 3) * v.w;
    const float y = m(1, 0) * v.x + m(1, 1) * v.y + m(1, 2) * v.z + m(1, 3) * v.w;
    const float z = m(2, 0) * v.x + m(2, 1) * v.y + m(2, 2) * v.z + m(2, 3) * v.w;
    const float w = m(3, 0) * v.x + m(3, 1) * v.y + m(3, 3) * v.z + m(3, 3) * v.w;

    v4 buffer;
    buffer.x = x;
    buffer.y = y;
    buffer.z = z;
    buffer.w = w;

    return buffer;
}

FORCEINLINE v3 project(const v3& a, const v3& b) {
    return b * (dot_product(a, b) / dot_product(b, b));
}

FORCEINLINE v3 reject(const v3& a, const v3& b) {
    return a - project(a, b);
}

FORCEINLINE v3 reflect(const v3& incident, const v3& normal) {
    return incident - normal * (2.0f * dot_product(incident, normal));
}

FORCEINLINE f32 angle_between(const v3& a, const v3& b) {
    return acos(clamp(dot_product(a, b) / (magnitude(a) * magnitude(b)), -1.0f, 1.0f));
}

FORCEINLINE m4 operator *(m4& a, m4& b)
{
    // Row 1
    const float m00 = a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0) + a(0, 2) * b(2, 0) + a(0, 3) * b(3, 0);
    const float m01 = a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1) + a(0, 2) * b(2, 1) + a(0, 3) * b(3, 1);
    const float m02 = a(0, 0) * b(0, 2) + a(0, 1) * b(1, 2) + a(0, 2) * b(2, 2) + a(0, 3) * b(3, 2);
    const float m03 = a(0, 0) * b(0, 3) + a(0, 1) * b(1, 3) + a(0, 2) * b(2, 3) + a(0, 3) * b(3, 3);

    // Row 2
    const float m10 = a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0) + a(1, 2) * b(2, 0) + a(1, 3) * b(3, 0);
    const float m11 = a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1) + a(1, 3) * b(3, 1);
    const float m12 = a(1, 0) * b(0, 2) + a(1, 1) * b(1, 2) + a(1, 2) * b(2, 2) + a(1, 3) * b(3, 2);
    const float m13 = a(1, 0) * b(0, 3) + a(1, 1) * b(1, 3) + a(1, 2) * b(2, 3) + a(1, 3) * b(3, 3);

    // Row 3
    const float m20 = a(2, 0) * b(0, 0) + a(2, 1) * b(1, 0) + a(2, 2) * b(2, 0) + a(2, 3) * b(3, 0);
    const float m21 = a(2, 0) * b(0, 1) + a(2, 1) * b(1, 1) + a(2, 2) * b(2, 1) + a(2, 3) * b(3, 1);
    const float m22 = a(2, 0) * b(0, 2) + a(2, 1) * b(1, 2) + a(2, 2) * b(2, 2) + a(2, 3) * b(3, 2);
    const float m23 = a(2, 0) * b(0, 3) + a(2, 1) * b(1, 3) + a(2, 2) * b(2, 3) + a(2, 3) * b(3, 3);

    // Row 4
    const float m30 = a(3, 0) * b(0, 0) + a(3, 1) * b(1, 0) + a(3, 2) * b(2, 0) + a(3, 3) * b(3, 0);
    const float m31 = a(3, 0) * b(0, 1) + a(3, 1) * b(1, 1) + a(3, 2) * b(2, 1) + a(3, 3) * b(3, 1);
    const float m32 = a(3, 0) * b(0, 2) + a(3, 1) * b(1, 2) + a(3, 2) * b(2, 2) + a(3, 3) * b(3, 2);
    const float m33 = a(3, 0) * b(0, 3) + a(3, 1) * b(1, 3) + a(3, 2) * b(2, 3) + a(3, 3) * b(3, 3);

    m4 buffer;
    buffer(0, 0) = m00;
    buffer(0, 1) = m01;
    buffer(0, 2) = m02;
    buffer(0, 3) = m03;

    buffer(1, 0) = m10;
    buffer(1, 1) = m11;
    buffer(1, 2) = m12;
    buffer(1, 3) = m13;

    buffer(2, 0) = m20;
    buffer(2, 1) = m21;
    buffer(2, 2) = m22;
    buffer(2, 3) = m23;

    buffer(3, 0) = m30;
    buffer(3, 1) = m31;
    buffer(3, 2) = m32;
    buffer(3, 3) = m33;

    return buffer;
}

const m4 identity = 
{ {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, 1, 0 },
    { 0, 0, 0, 1 }
} };

FORCEINLINE f32 m4_determinant(const m4& m)
{
    return m(0, 0) * m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1) +
        m(0, 1) * m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2) +
        m(2, 0) * m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);
}

FORCEINLINE m4 m4_scale(m4& m, v3& s)
{
    m4 scale =
    { {
        { s.x, 0, 0, 0 },
        { 0, s.y, 0, 0 },
        { 0, 0, s.z, 0 },
        { 0, 0, 0,   1 }
    } };

    m4 res = scale * m;
    return res;
}

FORCEINLINE m4 m4_translate(v3& t)
{
    m4 trans =
    { {
        { 0, 0, 0, t.x },
        { 0, 0, 0, t.y },
        { 0, 0, 0, t.z },
        { 0, 0, 0, 1 }
    } };

    return trans;
}

FORCEINLINE m4 m4_x_rotation(f32 a)
{
    const f32 c = cos(a);
    const f32 s = sin(a);

    m4 rotx =
    { {
        { 1, 0, 0, 0 },
        { 0, c,-s, 0 },
        { 0, s, c, 0 },
        { 0, 0, 0, 1 }
    } };

    return rotx;
}

FORCEINLINE m4 m4_y_rotation(f32 a)
{
    const f32 c = cos(a);
    const f32 s = sin(a);

    m4 roty =
    { {
        { c, 0, s, 0 },
        { 0, 1, 0, 0 },
        {-s, 0, c, 0 },
        { 0, 0, 0, 1 }
    } };

    return roty;
}

FORCEINLINE m4 m4_z_rotation(f32 a)
{
    const f32 c = cos(a);
    const f32 s = sin(a);

    m4 rotz =
    { {
        { c, -s,0, 0 },
        { s, c, 0, 0 },
        { 0, 0, 0, 0 },
        { 0, 0, 0, 1 }
    } };

    return rotz;
}

FORCEINLINE m4 m4_transpose(const m4* m) {
    m4 result;
    for (u32 i = 0; i < 4; i++) {
        for (u32 j = 0; j < 4; j++) {
            result(i, j) = m->table[j][i];
        }
    }
    return result;
}

FORCEINLINE bool m4_invert(const m4* m, m4* out) {
    // Compute matrix inversion (placeholder; use an optimized method)
    // Return true if successful; false if not invertible
    // (Actual implementation is complex and omitted for brevity)
    return false;
}


FORCEINLINE m4 frustum(f32 l, f32 r, f32 bot, f32 t, f32 np, f32 fp)
{
    const f32 a = 2.0f * np / (r - l);
    const f32 b = 2.0f * np / (t - bot);
    const f32 c = (r + l) / (r - l);
    const f32 d = (t + bot) / (t - bot);
    const f32 e = -(fp + np) / (fp - np);
    const f32 f = -l;
    const f32 g = -(2 * fp * np) / (fp - np);
    const f32 h = 0.0f;

    m4 mat =
    { {
        { a, 0, c, 0 },
        { 0, b, d, 0 },
        { 0, 0, e, g },
        { 0, 0, f, h }
    } };

    return mat;
}

FORCEINLINE bool point_inside_frustum(const Frustum* frustum, const v3* point) {
    for (int i = 0; i < 6; ++i) {
        if (dot_product(frustum->planes[i].normal, *point) + frustum->planes[i].d < 0) {
            return false;
        }
    }
    return true;
}

FORCEINLINE bool aabb_inside_frustum(const Frustum* frustum, const AABB* aabb) {
    for (int i = 0; i < 6; ++i) {
        const Plane* plane = &frustum->planes[i];
        v3 p = aabb->max;
        if (plane->normal.x < 0) p.x = aabb->min.x;
        if (plane->normal.y < 0) p.y = aabb->min.y;
        if (plane->normal.z < 0) p.z = aabb->min.z;
        if (dot_product(plane->normal, p) + plane->d < 0) {
            return false;
        }
    }
    return true;
}

FORCEINLINE m4 look_at(const v3& eye, const v3& target, const v3& up) {
    v3 z_axis = normalise(target - eye);
    v3 x_axis = normalise(cross_product(up, z_axis));
    v3 y_axis = cross_product(z_axis, x_axis);

    m4 view = identity;
    view(0, 0) = x_axis.x; view(0, 1) = x_axis.y; view(0, 2) = x_axis.z; view(0, 3) = -dot_product(x_axis, eye);
    view(1, 0) = y_axis.x; view(1, 1) = y_axis.y; view(1, 2) = y_axis.z; view(1, 3) = -dot_product(y_axis, eye);
    view(2, 0) = z_axis.x; view(2, 1) = z_axis.y; view(2, 2) = z_axis.z; view(2, 3) = -dot_product(z_axis, eye);
    return view;
}

FORCEINLINE m4 orthographic(f32 left, f32 right, f32 bottom, f32 top, f32 near, f32 far) {
    m4 ortho = identity;
    ortho(0, 0) = 2.0f / (right - left);
    ortho(1, 1) = 2.0f / (top - bottom);
    ortho(2, 2) = -2.0f / (far - near);
    ortho(0, 3) = -(right + left) / (right - left);
    ortho(1, 3) = -(top + bottom) / (top - bottom);
    ortho(2, 3) = -(far + near) / (far - near);
    return ortho;
}

FORCEINLINE m4 perspective(f32 fov, f32 aspect, f32 f, f32 b)
{
    f32 tangent = tanf(DEG2RAD((fov / 2)));
    f32 height = f * tangent;
    f32 width = height * aspect;

    m4 mat = frustum(-width, width, -height, height, f, b);

    return mat;
}

FORCEINLINE quat quat_from_axis_angle(v3 axis, f32 angle) {
    f32 half_angle = angle * 0.5f;
    f32 s = sin(half_angle);

    quat q = {
        axis.x * s,
        axis.y * s,
        axis.z * s,
        cos(half_angle)
    };

    return q;
}

FORCEINLINE m4 quat_to_matrix(const quat* q) {
    m4 result = identity;
    result(0, 0) = 1 - 2 * (q->y * q->y + q->z * q->z);
    result(0, 1) = 2 * (q->x * q->y - q->w * q->z);
    result(0, 2) = 2 * (q->x * q->z + q->w * q->y);

    result(1, 0) = 2 * (q->x * q->y + q->w * q->z);
    result(1, 1) = 1 - 2 * (q->x * q->x + q->z * q->z);
    result(1, 2) = 2 * (q->y * q->z - q->w * q->x);

    result(2, 0) = 2 * (q->x * q->z - q->w * q->y);
    result(2, 1) = 2 * (q->y * q->z + q->w * q->x);
    result(2, 2) = 1 - 2 * (q->x * q->x + q->y * q->y);

    return result;
}

FORCEINLINE quat quat_multiply(const quat* q1, const quat* q2) {
    quat result;
    result.x = q1->w * q2->x + q1->x * q2->w + q1->y * q2->z - q1->z * q2->y;
    result.y = q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x;
    result.z = q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w;
    result.w = q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
    return result;
}

FORCEINLINE quat quat_normalize(const quat* q) {
    f32 mag = sqrt(q->x * q->x + q->y * q->y + q->z * q->z + q->w * q->w);
    return (quat){ q->x / mag, q->y / mag, q->z / mag, q->w / mag };
}

FORCEINLINE quat quat_slerp(const quat* q1, const quat* q2, f32 t) {
    // Spherical linear interpolation
    f32 dot = q1->x * q2->x + q1->y * q2->y + q1->z * q2->z + q1->w * q2->w;
    const f32 threshold = 0.9995f;

    if (dot > threshold) {
        quat result = {
            q1->x + t * (q2->x - q1->x),
            q1->y + t * (q2->y - q1->y),
            q1->z + t * (q2->z - q1->z),
            q1->w + t * (q2->w - q1->w)
        };
        return quat_normalize(&result);
    }

    dot = fmax(fmin(dot, 1.0f), -1.0f);
    f32 theta_0 = acos(dot);
    f32 theta = theta_0 * t;

    quat q3 = { q2->x - q1->x * dot, q2->y - q1->y * dot, q2->z - q1->z * dot, q2->w - q1->w * dot };
    q3 = quat_normalize(&q3);

    return (quat){
        q1->x * cos(theta) + q3.x * sin(theta),
        q1->y * cos(theta) + q3.y * sin(theta),
        q1->z * cos(theta) + q3.z * sin(theta),
        q1->w * cos(theta) + q3.w * sin(theta)
    };
}

FORCEINLINE quat quat_inverse(const quat* q) {
    f32 norm_sq = q->x * q->x + q->y * q->y + q->z * q->z + q->w * q->w;
    return (quat){ -q->x / norm_sq, -q->y / norm_sq, -q->z / norm_sq, q->w / norm_sq };
}

FORCEINLINE f32 quat_dot(const quat* q1, const quat* q2) {
    return q1->x * q2->x + q1->y * q2->y + q1->z * q2->z + q1->w * q2->w;
}

FORCEINLINE void quat_to_axis_angle(const quat* q, v3* axis, f32* angle) {
    *angle = 2.0f * acos(q->w);
    f32 s = sqrt(1.0f - q->w * q->w);
    if (s > 1e-6f) {
        axis->x = q->x / s;
        axis->y = q->y / s;
        axis->z = q->z / s;
    } else {
        *axis = (v3){ 1.0f, 0.0f, 0.0f };
    }
}

FORCEINLINE bool aabb_intersects(const AABB* a, const AABB* b) {
    return (a->min.x <= b->max.x && a->max.x >= b->min.x) &&
           (a->min.y <= b->max.y && a->max.y >= b->min.y) &&
           (a->min.z <= b->max.z && a->max.z >= b->min.z);
}

FORCEINLINE bool ray_aabb_intersect(const Ray* ray, const AABB* aabb, f32* t) {
    f32 tmin = (aabb->min.x - ray->origin.x) / ray->direction.x;
    f32 tmax = (aabb->max.x - ray->origin.x) / ray->direction.x;

    if (tmin > tmax) {
        f32 temp = tmin;
        tmin = tmax;
        tmax = temp;
    }

    f32 tymin = (aabb->min.y - ray->origin.y) / ray->direction.y;
    f32 tymax = (aabb->max.y - ray->origin.y) / ray->direction.y;

    if (tymin > tymax) {
        f32 temp = tymin;
        tymin = tymax;
        tymax = temp;
    }

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    *t = tmin;
    return true;
}

FORCEINLINE bool ray_sphere_intersect(const Ray* ray, const Sphere* sphere, f32* t) {
    v3 oc = ray->origin - sphere->center;
    f32 a = dot_product(ray->direction, ray->direction);
    f32 b = 2.0f * dot_product(oc, ray->direction);
    f32 c = dot_product(oc, oc) - sphere->radius * sphere->radius;
    f32 discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    *t = (-b - sqrt(discriminant)) / (2.0f * a);
    return true;
}

FORCEINLINE bool plane_point_intersection(const Plane* plane, const v3* point) {
    // Check if the point lies on the plane
    return fabs(dot_product(plane->normal, *point) + plane->d) < 1e-6f;
}

FORCEINLINE f32 plane_distance_to_point(const Plane* plane, const v3* point) {
    // Distance from point to plane
    return dot_product(plane->normal, *point) + plane->d;
}

FORCEINLINE bool circle_contains_point(const Circle* circle, const v2* point) {
    // Check if the point lies within the circle
    v2 diff = *point - circle->center;
    return magnitude(diff) <= circle->radius;
}

FORCEINLINE bool circle_intersects(const Circle* a, const Circle* b) {
    // Check if two circles intersect
    v2 diff = a->center - b->center;
    f32 radii_sum = a->radius + b->radius;
    return magnitude(diff) <= radii_sum;
}

FORCEINLINE f32 triangle_area(const Triangle* tri) {
    // Compute area using cross product
    v3 ab = tri->vertices[1] - tri->vertices[0];
    v3 ac = tri->vertices[2] - tri->vertices[0];
    return 0.5f * magnitude(cross_product(ab, ac));
}

FORCEINLINE bool bounding_capsule_intersects(const BoundingCapsule* a, const BoundingCapsule* b) {
    // Simplified intersection: sphere distance + radii comparison
    v3 dir_a = a->pointB - a->pointA;
    v3 dir_b = b->pointB - b->pointA;
    v3 center_diff = (a->pointA + dir_a * 0.5f) - (b->pointA + dir_b * 0.5f);
    return magnitude(center_diff) <= (a->radius + b->radius);
}

FORCEINLINE bool bounding_ellipsoid_contains_point(const BoundingEllipsoid* ellipsoid, const v3* point) {
    // Check if a point lies inside the ellipsoid
    v3 relative_point = *point - ellipsoid->center;
    f32 scaled_x = relative_point.x / ellipsoid->radii.x;
    f32 scaled_y = relative_point.y / ellipsoid->radii.y;
    f32 scaled_z = relative_point.z / ellipsoid->radii.z;
    return (scaled_x * scaled_x + scaled_y * scaled_y + scaled_z * scaled_z) <= 1.0f;
}

FORCEINLINE bool obb_intersects(const OBB* a, const OBB* b) {
    // Simplified Separating Axis Theorem (SAT) for OBB intersection
    // Placeholder implementation (can be extended)
    return aabb_intersects((AABB*)a, (AABB*)b); // Treat OBB as AABB for simplicity
}

FORCEINLINE v3 bezier_curve_point(const BezierCurve* curve, f32 t) {
    // Calculate a point on the curve using Bernstein polynomials
    f32 one_minus_t = 1.0f - t;
    return curve->p0 * (one_minus_t * one_minus_t * one_minus_t) +
           curve->p1 * (3.0f * one_minus_t * one_minus_t * t) +
           curve->p2 * (3.0f * one_minus_t * t * t) +
           curve->p3 * (t * t * t);
}

FORCEINLINE v3 spline_point(const Spline* spline, f32 t) {
    // Get a point on the spline; assumes linear segments for simplicity
    u32 segment = (u32)(t * (spline->count - 1));
    f32 local_t = (t * (spline->count - 1)) - segment;
    v3 p0 = spline->control_points[segment];
    v3 p1 = spline->control_points[segment + 1];
    return p0 * (1.0f - local_t) + p1 * local_t;
}

FORCEINLINE bool bounding_cylinder_contains_point(const BoundingCylinder* cylinder, const v3* point) {
    // Check if point lies within the cylinder
    v3 relative = *point - cylinder->center;
    f32 height_check = fabs(relative.y) <= (cylinder->height * 0.5f);
    relative.y = 0; // Ignore height axis for radius check
    return height_check && magnitude(relative) <= cylinder->radius;
}

FORCEINLINE u32 grid_space_index(const GridSpace* grid, const v2i* cell) {
    // Calculate a flattened index for the grid space
    return cell->x + cell->y * grid->grid_dimensions.x;
}

FORCEINLINE bool uv_box_contains(const UVBox* box, const UV* uv) {
    // Check if UV coordinates lie within the box
    return (uv->u >= box->min.x && uv->u <= box->max.x &&
            uv->v >= box->min.y && uv->v <= box->max.y);
}

FORCEINLINE f32 line_length(const Line* line) {
    return magnitude(line->end - line->start);
}

FORCEINLINE bool line_intersects_plane(const Line* line, const Plane* plane, v3* intersection) {
    f32 denom = dot_product(plane->normal, line->end - line->start);
    if (fabs(denom) < 1e-6f) return false; // Parallel

    f32 t = -(dot_product(plane->normal, line->start) + plane->d) / denom;
    if (t < 0.0f || t > 1.0f) return false; // Outside segment

    *intersection = line->start + (line->end - line->start) * t;
    return true;
}

FORCEINLINE f32 lerp(f32 a, f32 b, f32 t) {
    return a * (1.0f - t) + b * t;
}

FORCEINLINE v3 lerp_v3(const v3* a, const v3* b, f32 t) {
    return *a * (1.0f - t) + *b * t;
}

FORCEINLINE f32 cubic_interpolate(f32 p0, f32 p1, f32 p2, f32 p3, f32 t) {
    f32 t2 = t * t;
    f32 t3 = t2 * t;
    return 0.5f * ((2.0f * p1) +
                   (-p0 + p2) * t +
                   (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                   (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}