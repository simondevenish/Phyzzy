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

#include <stdbool.h>

#include "math_types.h"

// Physical material properties
typedef struct {
    f32 friction;      // Coefficient of friction (0 = no friction, 1 = high friction)
    f32 restitution;   // Bounciness (0 = no bounce, 1 = perfect bounce)
    f32 density;       // Mass per unit volume
} PhysicsMaterial;

// Rigid body dynamics
typedef struct {
    v3 velocity;       // Linear velocity
    v3 angular_velocity; // Rotational velocity
    f32 mass;          // Mass of the object
    f32 inverse_mass;  // Precomputed inverse mass (0 for immovable objects)
    m3 inertia_tensor; // Rotational inertia tensor
    m3 inverse_inertia_tensor; // Precomputed inverse of the inertia tensor
} RigidBody;

// Collider types
typedef enum {
    COLLIDER_TYPE_SPHERE,
    COLLIDER_TYPE_AABB,
    COLLIDER_TYPE_OBB,
    COLLIDER_TYPE_PLANE,
    COLLIDER_TYPE_CAPSULE,
    COLLIDER_TYPE_MESH
} ColliderType;

typedef struct {
    ColliderType type;
    PhysicsMaterial material;
    Transform transform; // Position, rotation, and scale of the collider
    union {
        Sphere sphere;
        AABB aabb;
        OBB obb;
        Plane plane;
        BoundingCapsule capsule;
        Mesh mesh;
    };
} Collider;

typedef struct {
    bool enabled;      // Is wobble enabled
    f32 amplitude;     // Wobble amplitude
    f32 frequency;     // Wobble frequency
    v3 direction;      // Wobble oscillation direction
    f32 time;          // Internal time tracker for wobble state
} WobbleState;

// Physics body combining RigidBody and Collider
typedef struct {
    RigidBody* rigid_body;   // Pointer to the associated rigid body
    Collider* collider;      // Pointer to the collider
    bool is_static;          // True if the object is immovable
    WobbleState wobble;      // Wobble and oscillation state
} PhysicsBody;

// Force generators (e.g., gravity, springs)
typedef struct {
    v3 force; // Direction and magnitude of the force
} Force;

typedef struct {
    PhysicsBody* body_a; // First body connected by the spring
    PhysicsBody* body_b; // Second body connected by the spring
    f32 rest_length;     // Resting length of the spring
    f32 stiffness;       // How stiff the spring is
    f32 damping;         // Damping factor to reduce oscillations
} Spring;

// Collision data
typedef struct {
    PhysicsBody* body_a;    // First body in the collision
    PhysicsBody* body_b;    // Second body in the collision
    v3 contact_point;       // Point of contact
    v3 contact_normal;      // Normal at the point of contact
    f32 penetration_depth;  // How far the objects are overlapping
} Collision;

// Physics constraints (e.g., fixed joints, hinges)
typedef enum {
    CONSTRAINT_TYPE_FIXED,
    CONSTRAINT_TYPE_HINGE,
    CONSTRAINT_TYPE_SLIDER
} ConstraintType;

typedef struct {
    ConstraintType type;
    PhysicsBody* body_a;
    PhysicsBody* body_b;
    Transform local_anchor_a; // Anchor point relative to body A
    Transform local_anchor_b; // Anchor point relative to body B
} Constraint;

// Environment forces (e.g., gravity, wind)
typedef struct {
    v3 gravity; // Global gravity force (e.g., {0, -9.81, 0})
    v3 wind;    // Wind direction and strength
} EnvironmentForces;

// Physics settings for the simulation
typedef struct {
    f32 time_step;           // Time step for physics updates
    u32 max_iterations;      // Maximum solver iterations for constraints
    bool enable_gravity;     // Whether gravity is enabled
    v3 gravity;              // Gravity vector
    u32 max_bodies;          // Maximum number of physics bodies
    u32 max_colliders;       // Maximum number of colliders
    u32 max_joints;          // Maximum number of joints
    f32 collision_tolerance; // Collision detection tolerance
    bool enable_sleeping;    // Whether inactive bodies can sleep
    f32 sleep_threshold;     // Velocity threshold below which a body can sleep
    f32 linear_damping;      // Default linear damping applied to all bodies
    f32 angular_damping;     // Default angular damping applied to all bodies
    u32 max_zones;           // Maximum number of physics zones
    bool enable_debug_logs;  // Whether debug logging is enabled
} PhysicsSettings;


// Raycast result
typedef struct {
    bool hit;            // True if the ray hit something
    PhysicsBody* body;   // The body that was hit
    v3 hit_point;        // The point of intersection
    v3 hit_normal;       // The normal at the hit point
    f32 hit_distance;    // Distance from ray origin to hit point
} PhysicsRaycastHit;

// Soft body physics
typedef struct {
    v3* vertices;        // Array of vertex positions
    v3* velocities;      // Array of vertex velocities
    f32* masses;         // Mass of each vertex
    u32 vertex_count;    // Number of vertices in the soft body
    u32* links;          // Indices of connected vertices
    u32 link_count;      // Number of links
    f32 stiffness;       // Stiffness of the soft body
    f32 damping;         // Damping factor to reduce oscillations
} SoftBody;

// Fluid simulation
typedef struct {
    v3* particles;       // Positions of fluid particles
    v3* velocities;      // Velocities of fluid particles
    f32* densities;      // Densities of fluid particles
    f32* pressures;      // Pressures of fluid particles
    u32 particle_count;  // Number of fluid particles
    f32 viscosity;       // Viscosity of the fluid
    f32 rest_density;    // Rest density of the fluid
    f32 smoothing_radius; // Radius for calculating forces between particles
} Fluid;

// Particle system (simplified physics)
typedef struct {
    v3* positions;       // Positions of particles
    v3* velocities;      // Velocities of particles
    Color* colors;       // Colors of particles
    f32* lifetimes;      // Remaining lifetime of each particle
    u32 particle_count;  // Number of particles in the system
    f32 emission_rate;   // Rate of particle emission
} ParticleSystem;

typedef struct {
    PhysicsBody* body;     // Associated physics body
    bool on_ground;        // True if the character is on the ground
    f32 max_speed;         // Maximum movement speed
    f32 jump_force;        // Force applied during a jump
} CharacterPhysics;

typedef struct {
    v3 center;       // Center of the interaction zone
    f32 radius;      // Radius for spherical zones
    AABB bounds;     // Optional bounding box for rectangular zones
    bool is_active;  // True if the zone is currently active
} InteractionZone;

typedef struct {
    v3 position;         // Position of the object
    v3 velocity;         // Velocity of the object
    quat rotation;       // Rotation of the object
    u32 object_id;       // Unique ID for the object
    bool is_static;      // True if the object is immovable
} NetworkedPhysicsState;

typedef struct {
    v3 position;       // Current position
    v3 velocity;       // Current velocity
    v3 acceleration;   // Current acceleration
    f32 steering_angle; // Current steering angle
    f32 wheel_friction; // Friction applied to wheels
    f32 mass;          // Mass of the vehicle
} VehiclePhysics;

typedef struct {
    v3* track_points;  // Points defining the track path
    u32 point_count;   // Number of points
    f32 speed;         // Speed of travel along the track
    f32 current_t;     // Current position as a parameter (0 to 1)
} TrackSystem;

typedef struct {
    v3 position;
    v3 velocity;
    bool is_anchor;  // New flag to mark a segment as anchored
} RopeSegment;

typedef struct {
    RopeSegment* segments;    // Array of rope segments
    u32 segment_count;        // Number of segments
    f32 segment_length;       // Length of each segment
    f32 stiffness;            // Stiffness of the rope
    f32 damping;              // Damping to stabilize simulation
    PhysicsBody* attached_body; // Optional attached body
    bool is_dynamic;          // True if the rope can dynamically change segments
} Rope;

typedef struct {
    v3 position;       // Center of the water body
    v2 size;           // Dimensions of the water body
    f32 wave_height;   // Maximum height of waves
    f32 wave_speed;    // Speed of wave propagation
    f32 buoyancy;      // Buoyancy factor for objects in water
} WaterBody;

typedef struct {
    v3 direction;      // Direction of the wind
    f32 strength;      // Wind strength
    v3 bounds_min;     // Minimum bounds for the wind effect area
    v3 bounds_max;     // Maximum bounds for the wind effect area
} WindZone;

typedef struct {
    AABB bounds;       // Bounds of the grid cell
    PhysicsBody** bodies; // Pointers to bodies in this cell
    u32 body_count;    // Number of bodies in this cell
} SpatialGridCell;

typedef struct {
    SpatialGridCell* cells; // Array of grid cells
    u32 cell_count;         // Number of cells
    v3 grid_min;            // Minimum bounds of the grid
    v3 grid_max;            // Maximum bounds of the grid
    v3 cell_size;           // Size of each cell
} SpatialGrid;

typedef struct {
    PhysicsBody** pairs;   // Pairs of potentially colliding bodies
    u32 pair_count;        // Number of pairs
} BroadphasePairCache;

typedef struct {
    v3 position;      // Position on the surface
    v3 normal;        // Surface normal for orientation
    UV uv_coordinates; // UV coordinates on the surface
    u32 texture_id;   // ID of the texture to render the sticker
    f32 scale;        // Scale of the sticker
    Color tint;       // Optional tint to modify the sticker's color
    bool is_permanent; // If true, sticker cannot be removed
} Sticker;

typedef enum {
    PHYSICS_TYPE_RIGID,
    PHYSICS_TYPE_SOFT,
    PHYSICS_TYPE_RUBBERY
} WidgetPhysicsType;

typedef struct {
    WidgetPhysicsType type;    // Type of physical behavior
    PhysicsMaterial material;  // Material properties (friction, restitution, etc.)
    RigidBody rigid_body;      // Optional rigid body for dynamic widgets
    v3 deformation;            // Soft body deformation (applied for soft/rubbery widgets)
} WidgetPhysics;

typedef struct {
    Transform attachment_point; // Position and orientation relative to the player
    WidgetPhysics physics;      // Physics behavior of the widget
    bool follows_movement;      // If true, moves with the player
    bool has_physics_simulation; // Whether the widget participates in physics simulation
    char name[64];              // Name or ID of the widget for identification
} WidgetAttachment;

typedef struct {
    v3 position;      // Surface position
    v3 normal;        // Surface normal
    bool is_attachable; // True if objects can attach to this surface
    bool is_sticky;    // True if stickers or widgets adhere permanently
    u32 surface_id;   // Unique ID for surface identification
} SurfaceProperties;

typedef struct {
    PhysicsBody* body_a;     // First body
    PhysicsBody* body_b;     // Second body
    f32 rest_length;         // Resting length of the elastic constraint
    f32 stiffness;           // Stiffness of the elastic band
    f32 damping;             // Damping factor
    f32 max_stretch_length;  // Maximum allowable stretch
} ElasticConstraint;

typedef struct {
    f32 max_force;     // Maximum force before detachment
    f32 max_torque;    // Maximum torque before detachment
    bool can_reattach; // If true, reattachment is possible
} DetachmentRules;

typedef enum {
    LAYER_DEFAULT,
    LAYER_STICKER,
    LAYER_WIDGET,
    LAYER_PLAYER,
    LAYER_ENVIRONMENT
} PhysicsLayer;

typedef struct {
    PhysicsLayer layer;  // Physics layer of the object
    u32 collision_mask;  // Mask defining which layers this object collides with
} LayerSettings;

typedef struct {
    f32 max_health;    // Maximum health of the object
    f32 current_health; // Current health
    bool is_destructible; // Whether the object can be destroyed
    bool is_repairable;   // Whether the object can be repaired
} DestructibleObject;

typedef struct {
    f32 buoyancy_factor; // How much the object floats
    f32 drag_coefficient; // Resistance while moving through the fluid
} FluidInteraction;

typedef struct {
    v3 position;  // Center of the deformation zone
    f32 radius;   // Radius of influence
    f32 deformation_strength; // How much deformation occurs
} DeformationZone;

typedef struct {
    v3 direction;       // Direction of the wind
    f32 strength;       // Wind strength
    f32 turbulence;     // Variability in wind direction and strength
} AdvancedWindZone;

typedef struct {
    InteractionZone base_zone; // Base interaction zone
    v3 velocity;               // Movement speed and direction
    f32 expansion_rate;        // Rate at which the zone expands or contracts
    bool is_dynamic;           // If true, the zone changes over time
} DynamicInteractionZone;

typedef struct {
    PhysicsBody** bones;      // Array of physics bodies representing bones
    u32 bone_count;           // Number of bones in the ragdoll
    Constraint* joint_constraints; // Constraints defining joint limits
    u32 constraint_count;     // Number of constraints
    bool is_active;           // True if the ragdoll is currently active
} Ragdoll;

typedef struct {
    AABB bounds;          // Bounds of the zone
    v3 gravity_override;  // Override for gravity in this zone
    f32 friction_override; // Override for friction
    bool affects_all;      // True if the zone affects all objects
} PhysicsZone;

typedef struct {
    Sticker base_sticker;  // Sticker details
    v3 deformation_vector; // Deformation applied to the sticker
    bool scales_with_surface; // Whether the sticker scales with surface changes
} DecalPhysics;