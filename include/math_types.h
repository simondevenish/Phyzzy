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

#include "common.h"

#define PI 3.14159265359f
#define TAU 6.28318530718f
#define DEG2RAD(x) ((x) * (PI / 180.0f))
#define RAD2DEG(x) ((x) * (180.0f / PI))

typedef union v2 v2;
typedef union v3 v3;
typedef union v4 v4;
typedef union m4 m4;

union v2
{
  f32 v[2];
  struct { f32 x, y; };
};

union v3
{
  f32 v[3];
  struct { f32 x, y, z; };
};

union v4
{
  f32 v[4];
  struct { f32 x, y, z, w; };
};

typedef struct {
    f32 m[3][3];
} m3;


union m4
{
  f32 table[4][4];
  v4 v[4];
};

typedef struct {
    v3 normal;  // Normal vector of the plane
    f32 d;      // Distance from the origin
} Plane;

typedef struct {
    v2 center;
    f32 radius;
} Circle;

typedef struct {
    v3 center;
    f32 radius;
} Sphere;

typedef struct {
    v3 start;
    v3 end;
} Line;

typedef struct {
    v3 position;
    quat rotation;  // Quaternion rotation
    v3 scale;       // Scale factor
} Transform;

typedef struct {
    f32 r, g, b, a;
} Color;

typedef struct {
    f32 u, v;
} UV;

typedef struct {
    f32 x, y, z, w;
} quat;

typedef struct {
    v3 min;
    v3 max;
} AABB;

typedef struct {
    v3 origin;
    v3 direction; // Normalized
} Ray;

// Composite Structures
typedef enum {
    BOUNDING_TYPE_AABB,
    BOUNDING_TYPE_SPHERE
} BoundingType;

typedef struct {
    BoundingType type;
    union {
        AABB aabb;
        Sphere sphere;
    };
} BoundingVolume;

typedef struct {
    v3 position;
    v3 normal;
    UV texcoord;
} Vertex;

typedef struct {
    Plane planes[6];  // Left, right, top, bottom, near, far
} Frustum;

typedef struct {
    v3 vertices[3];
} Triangle;

typedef struct {
    v2 min;
    v2 max;
} Rectangle;

typedef struct {
    v3 pointA;
    v3 pointB;
    f32 radius;
} BoundingCapsule;

typedef struct {
    v3 vertices[8];  // 8 vertices of the hexahedron
} Hexahedron;

typedef struct {
    u32 indices[3];
} Index;

typedef struct {
    v3* vertices;
    u32 vertex_count;
} Polygon;

typedef struct {
    v3 position;
    v3 normal;
    UV texcoord;
    v3 tangent;
    v3 bitangent;
} VertexWithTangents;

typedef struct {
    u32 seed;
    f32 frequency;
    f32 amplitude;
} NoiseSettings;

typedef struct {
    v3 center;
    v3 radii;  // Radius along x, y, z axes
} BoundingEllipsoid;

typedef struct {
    v3 position;  // Center position of the cell
    v3 size;      // Dimensions of the cell
} GridCell;

typedef struct {
    v3 p0;  // Start point
    v3 p1;  // Control point 1
    v3 p2;  // Control point 2
    v3 p3;  // End point
} BezierCurve;

typedef struct {
    v3* control_points;
    u32 count;
} Spline;

typedef struct {
    v2 position;  // Position of the viewport (UI space)
    v2 size;      // Size of the viewport (UI space)
    f32 aspect_ratio; // Aspect ratio of the camera
} Viewport;

typedef struct {
    v3 center;
    f32 radius;
    f32 height;
} BoundingCylinder;

typedef struct {
    Vertex* vertices;   // Array of vertices
    u32 vertex_count;
    Index* indices;     // Array of indices
    u32 index_count;
} Mesh;

typedef struct {
    v2i grid_dimensions;  // Dimensions of the grid (rows and columns)
    f32 cell_size;        // Size of each grid cell
} GridSpace;

typedef struct {
    v3 point;   // Point of intersection
    v3 normal;  // Surface normal at intersection
    f32 distance; // Distance from ray origin
    u32 hit_id; // ID of the hit object (if applicable)
} RaycastHit;

typedef struct TransformNode {
    Transform transform;
    struct TransformNode* parent;
    struct TransformNode* children;
    u32 child_count;
} TransformNode;

typedef struct {
    m4* matrices;   // Array of matrices
    u32 count;      // Number of matrices currently in the stack
    u32 capacity;   // Capacity of the stack
} MatrixStack;

typedef struct {
    v3 center;      // Center of the box
    v3 half_extents; // Half-size in each dimension
    m3 rotation;    // Rotation matrix defining the orientation
} OBB;

typedef struct {
    v3 position;  // Position of the camera
    v3 direction; // Forward direction (normalized)
    v3 up;        // Up direction (normalized)
    f32 fov;      // Field of view (degrees)
    f32 aspect_ratio;
    f32 near_clip;
    f32 far_clip;
} Camera;

typedef struct {
    v2 center;
    f32 radius;
} BoundingDisk;

typedef struct {
    v3 start;
    v3 end;
} Edge;

typedef struct {
    quat real;
    quat dual;
} DualQuaternion;

typedef struct {
    m4 view_matrix;
    m4 projection_matrix;
    Viewport viewport;
} Projection;

typedef struct {
    v2 min;
    v2 max;
} UVBox;

typedef struct {
    v2 center;
    v2 size;
} GridCell2D;

typedef struct {
    v3 start;
    v3 control_point; // Middle control point
    v3 end;
} CurveSegment;

typedef struct { i32 x, y; } v2i;
typedef struct { i32 x, y, z; } v3i;
typedef struct { i32 x, y, z, w; } v4i;

