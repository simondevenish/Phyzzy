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

#include "physics_types.h"
#include "vector_math.h"
#include "types.h"

// Physics system lifecycle
void physics_initialize(PhysicsSettings* settings);
void physics_cleanup();

// Physics body management
void physics_body_destroy(PhysicsBody* body);  // Properly cleanup and deallocate a physics body
void physics_body_set_mass(RigidBody* body, f32 mass);  // Dynamically change the mass of a body
void physics_body_get_mass(const RigidBody* body, f32* mass);  // Retrieve the current mass of a body
void physics_body_set_velocity(RigidBody* body, const v3* velocity);  // Set the linear velocity
void physics_body_get_velocity(const RigidBody* body, v3* velocity);  // Retrieve the linear velocity
void physics_body_set_angular_velocity(RigidBody* body, const v3* angular_velocity);  // Set rotational velocity
void physics_body_get_angular_velocity(const RigidBody* body, v3* angular_velocity);  // Retrieve rotational velocity
void physics_body_reset_forces(RigidBody* body);  // Reset all accumulated forces and torques
void physics_body_translate(PhysicsBody* body, const v3* translation);  // Move a physics body directly
void physics_body_rotate(PhysicsBody* body, const quat* rotation);  // Rotate a physics body directly
void physics_body_set_position(PhysicsBody* body, const v3* position);  // Explicitly set the position
void physics_body_get_position(const PhysicsBody* body, v3* position);  // Retrieve the current position
void physics_body_set_static(PhysicsBody* body, bool is_static);  // Dynamically toggle static/dynamic state
bool physics_body_is_static(const PhysicsBody* body);  // Check if a body is currently static
void physics_body_apply_damping(RigidBody* body, f32 linear_damping, f32 angular_damping);  // Apply damping to reduce velocity over time
void physics_body_apply_drag(RigidBody* body, f32 drag_coefficient);  // Apply drag proportional to velocity
void physics_body_freeze(PhysicsBody* body);  // Freeze a physics body to stop updates
void physics_body_unfreeze(PhysicsBody* body);  // Unfreeze a physics body to resume updates

// Surface properties
void surface_apply_friction(SurfaceProperties* surface, PhysicsBody* body, f32 delta_time);  // Apply surface-specific friction
void surface_apply_bounce(SurfaceProperties* surface, PhysicsBody* body, const v3* contact_point);  // Handle restitution or bouncing effects
void surface_apply_damping(SurfaceProperties* surface, PhysicsBody* body, f32 delta_time);  // Apply linear/angular damping due to surface interaction
void surface_check_detachment(SurfaceProperties* surface, PhysicsBody* body);  // Determine if an object detaches from a sticky surface
void surface_apply_velocity_constraints(SurfaceProperties* surface, PhysicsBody* body, f32 delta_time);  // Limit velocities on special surfaces (e.g., conveyor belts)
void surface_apply_custom_behavior(SurfaceProperties* surface, PhysicsBody* body, f32 delta_time);  // For game-specific behaviors like ice or quicksand
void surface_register_event(SurfaceProperties* surface, PhysicsBody* body, const char* event_name);  // Notify when specific interactions occur
void surface_reset(SurfaceProperties* surface);  // Reset all properties of the surface to defaults
void surface_query_properties(const SurfaceProperties* surface, PhysicsMaterial* material_out);  // Retrieve surface material properties for debugging or gameplay


// Wobble and oscillation effects
void physics_body_apply_directional_wobble(PhysicsBody* body, const v3* direction, f32 amplitude, f32 frequency);  // Directional wobble effects
void physics_body_sync_wobble_with_animation(PhysicsBody* body, f32 animation_phase);  // Synchronize wobble with an external animation system
void physics_body_apply_impact_wobble(PhysicsBody* body, const v3* impact_force, f32 decay_rate);  // Wobble caused by external impacts
void physics_body_get_wobble_state(const PhysicsBody* body, f32* amplitude_out, f32* frequency_out);  // Retrieve current wobble properties
void physics_body_toggle_wobble(PhysicsBody* body, bool enable);  // Enable or disable wobble dynamically


// Rope physics
void rope_set_material_properties(Rope* rope, f32 stiffness, f32 damping);  // Adjust rope's physical properties dynamically
void rope_cut_segment(Rope* rope, u32 segment_index);  // Cut the rope at a specific segment
void rope_reconnect_segments(Rope* rope, u32 segment_index_a, u32 segment_index_b);  // Reconnect two rope segments
void rope_check_collision_with_body(const Rope* rope, const PhysicsBody* body, u32* colliding_segment);  // Detect collisions with a physics body
void rope_apply_environment_forces(Rope* rope, const EnvironmentForces* forces);  // Apply wind or other environmental forces to the rope
void rope_add_segment(Rope* rope, const v3* position);  // Add a new segment to the rope dynamically
void rope_remove_segment(Rope* rope, u32 segment_index);  // Remove a specific segment from the rope
void rope_get_tension(const Rope* rope, f32* max_tension, u32* max_tension_segment);  // Calculate the maximum tension in the rope
void rope_simulate_swing(Rope* rope, PhysicsBody* anchor_body, const v3* applied_force);  // Simulate swinging motion for the rope
void rope_break_at_tension_limit(Rope* rope, f32 tension_limit);  // Automatically break the rope if tension exceeds the limit

// Collision detection
bool check_collision(const Collider* collider_a, const Collider* collider_b, Collision* collision_data);
bool check_aabb_collision(const AABB* a, const AABB* b);
bool check_obb_collision(const OBB* a, const OBB* b);
bool check_capsule_collision(const BoundingCapsule* a, const BoundingCapsule* b);

// Advanced collision detection
bool check_sphere_collision(const Sphere* a, const Sphere* b);  // Sphere-sphere collision detection
bool check_plane_collision(const Plane* plane, const Collider* collider);  // Plane-collider collision detection
bool check_point_in_collider(const Collider* collider, const v3* point);  // Point containment within a collider
bool check_ray_collision(const Collider* collider, const Ray* ray, PhysicsRaycastHit* hit_data);  // Raycasting against colliders
bool check_mesh_collision(const Mesh* mesh_a, const Mesh* mesh_b);  // Mesh-mesh collision detection

// Collision response
void resolve_collision(Collision* collision_data);  // Resolve collision with appropriate responses
void resolve_collision_with_materials(Collision* collision_data, const PhysicsMaterial* material_a, const PhysicsMaterial* material_b);  // Material-based collision response
void calculate_collision_impulse(const Collision* collision_data, v3* impulse);  // Calculate impulse based on collision

// Swept and continuous collision detection
bool swept_collision_test(const Collider* collider_a, const v3* velocity_a, const Collider* collider_b, const v3* velocity_b, Collision* collision_data);  // Continuous collision detection
bool broadphase_collision_test(const SpatialGrid* grid, const Collider* collider, Collision* collision_results, u32* result_count);  // Broadphase collision testing

// Utility functions
bool check_collision_with_layer(const Collider* collider, PhysicsLayer layer);  // Check collision against specific physics layers
bool is_colliding_with_static_object(const Collider* collider);  // Determine if a collider is in contact with a static object
u32 find_nearest_colliders(const Collider* collider, Collider** nearby_colliders, u32 max_results);  // Find nearby colliders
bool calculate_contact_manifold(const Collider* collider_a, const Collider* collider_b, Collision* collision_data, v3* contact_points, u32* contact_count);  // Generate detailed contact points

// Environment forces
void apply_environment_forces(PhysicsBody* body, const EnvironmentForces* environment);
void apply_physics_zone(const PhysicsZone* zone, PhysicsBody* body);

// General environmental forces
void apply_gravity(PhysicsBody* body, const v3* gravity);  // Apply global or custom gravity to a body
void apply_wind_force(PhysicsBody* body, const EnvironmentForces* wind);  // Apply wind force to a body
void apply_buoyancy_force(PhysicsBody* body, const WaterBody* water_body);  // Apply buoyancy for objects in water
void apply_custom_environment_force(PhysicsBody* body, const v3* force);  // Apply a custom force independent of predefined zones

// Zone-specific effects
void physics_zone_apply_gravity_override(const PhysicsZone* zone, PhysicsBody* body);  // Apply gravity override in a zone
void physics_zone_apply_friction_override(const PhysicsZone* zone, PhysicsBody* body);  // Apply friction override in a zone
void physics_zone_apply_turbulence(const PhysicsZone* zone, PhysicsBody* body, f32 intensity);  // Simulate turbulence within a zone

// Environmental dynamics
void simulate_wind_field(const EnvironmentForces* wind, const SpatialGrid* grid, f32 delta_time);  // Simulate wind effects across a grid
void update_water_body_forces(WaterBody* water_body, PhysicsBody** bodies, u32 body_count, f32 delta_time);  // Update forces for objects in water
void update_environment_effects(const EnvironmentForces* environment, PhysicsBody** bodies, u32 body_count, f32 delta_time);  // Update all active environmental effects

// Noise and variability
void apply_noise_to_environment_force(PhysicsBody* body, const v3* base_force, f32 noise_amplitude);  // Add noise to forces for variability
void apply_turbulence_to_wind(EnvironmentForces* wind, f32 delta_time);  // Introduce turbulence into wind forces

// Utilities
bool is_body_within_physics_zone(const PhysicsBody* body, const PhysicsZone* zone);  // Check if a body is inside a physics zone
bool is_body_affected_by_environment(const PhysicsBody* body, const EnvironmentForces* environment);  // Determine if a body is influenced by an environmental force

// Dynamic floor effects
void physics_body_update_surface_interaction(PhysicsBody* body, const PhysicsZone* floor_zone);
void apply_friction_for_surface(PhysicsBody* body, const PhysicsZone* floor_zone);
void apply_surface_restitution(PhysicsBody* body, const PhysicsZone* floor_zone);
void configure_surface_properties(PhysicsZone* floor_zone, f32 friction, f32 restitution);
void configure_dynamic_surface_effects(PhysicsZone* floor_zone, bool enable_ice, bool enable_stickiness, bool enable_wetness);
void physics_body_adjust_for_wet_surface(PhysicsBody* body, const PhysicsZone* floor_zone);
void physics_body_apply_ice_friction(PhysicsBody* body, const PhysicsZone* floor_zone);

// Advanced floor interaction
void physics_body_adjust_for_slippery_surface(PhysicsBody* body, const PhysicsZone* floor_zone);  // Adjust dynamics for slippery surfaces like polished floors
void physics_body_interact_with_sticky_surface(PhysicsBody* body, const PhysicsZone* floor_zone);  // Simulate stickiness interaction (e.g., glue or tar)
void physics_body_simulate_rough_surface(PhysicsBody* body, const PhysicsZone* floor_zone, f32 roughness_factor);  // Add resistance for rough terrains

// Dynamic property updates
void update_dynamic_surface_properties(PhysicsZone* floor_zone, f32 delta_time);  // Update properties like wetness or stickiness over time
void physics_zone_toggle_surface_effects(PhysicsZone* floor_zone, bool enable_dynamic_friction, bool enable_dynamic_restitution);  // Enable or disable dynamic surface effects

// Environmental effects on surfaces
void apply_temperature_to_surface(PhysicsZone* floor_zone, f32 temperature);  // Simulate temperature changes affecting friction (e.g., ice melting)
void apply_surface_degradation(PhysicsZone* floor_zone, f32 delta_time);  // Simulate wear and tear on the surface over time

// Combined floor effects
void physics_body_apply_combined_surface_effects(PhysicsBody* body, const PhysicsZone* floor_zone);  // Apply all relevant surface effects simultaneously
void physics_body_check_surface_transition(PhysicsBody* body, const PhysicsZone* previous_zone, const PhysicsZone* current_zone);  // Handle transitions between different surface types

// Surface utilities
bool is_surface_wet(const PhysicsZone* floor_zone);  // Check if a surface is wet
bool is_surface_sticky(const PhysicsZone* floor_zone);  // Check if a surface is sticky
bool is_surface_slippery(const PhysicsZone* floor_zone);  // Check if a surface is icy or slippery
bool does_surface_allow_movement(const PhysicsZone* floor_zone, const PhysicsBody* body);  // Determine if movement is allowed on the surface

// Ropes and soft connections
void rope_initialize(Rope* rope, const v3* points, u32 segment_count, f32 segment_length);
void rope_update(Rope* rope, f32 delta_time);
void rope_apply_force(Rope* rope, const v3* force);
void rope_attach_to_body(Rope* rope, PhysicsBody* body, u32 segment_index);
void rope_detach_from_body(Rope* rope, u32 segment_index);
void rope_add_segment(Rope* rope, const v3* position);
void rope_remove_segment(Rope* rope, u32 segment_index);

// Advanced rope dynamics
void rope_simulate_stretch(Rope* rope, f32 stretch_limit, f32 stiffness);  // Simulate stretching and snapping behavior
void rope_apply_damping(Rope* rope, f32 damping_factor);  // Reduce oscillations in the rope over time
void rope_adjust_tension(Rope* rope, f32 tension_factor);  // Adjust rope tension dynamically

// Collision and interaction
void rope_check_collisions(Rope* rope, const PhysicsBody** bodies, u32 body_count);  // Handle rope collisions with other objects
void rope_interact_with_environment(Rope* rope, const EnvironmentForces* forces);  // Simulate wind or other environmental effects on the rope

// Rope anchor
void rope_attach_to_anchor(Rope* rope, const v3* anchor_point);  // Attach rope to a static anchor point in the environment
void rope_detach_from_anchor(Rope* rope);  // Detach the rope from its anchor point

// Utility functions
bool rope_is_stretched(const Rope* rope, f32 threshold);  // Check if the rope is stretched beyond a certain limit
f32 rope_get_total_length(const Rope* rope);  // Calculate the total length of the rope
void rope_recalculate_segments(Rope* rope);  // Recalculate segment positions after a major adjustment

// Soft connection enhancements
void rope_connect_soft_body(Rope* rope, SoftBody* soft_body, u32 soft_body_vertex_index);  // Connect a rope to a soft body vertex
void rope_disconnect_soft_body(Rope* rope, SoftBody* soft_body, u32 soft_body_vertex_index);  // Disconnect a rope from a soft body
void rope_simulate_bounce(Rope* rope, f32 bounce_factor);  // Simulate bouncy behavior for elastic ropes

// Signal-related stickiness
void interaction_zone_sticky_trigger(DynamicInteractionZone* zone, PhysicsBody* body, bool is_attaching);
void interaction_zone_dynamic_update(DynamicInteractionZone* zone, f32 delta_time);

// Advanced interaction zone stickiness behavior
void interaction_zone_set_stickiness(DynamicInteractionZone* zone, f32 stickiness_factor);  // Adjust stickiness dynamically
void interaction_zone_reset_stickiness(DynamicInteractionZone* zone);  // Reset stickiness to default values

// Trigger and event enhancements
void interaction_zone_handle_attachment(DynamicInteractionZone* zone, PhysicsBody* body);  // Manage object attachment behavior
void interaction_zone_handle_detachment(DynamicInteractionZone* zone, PhysicsBody* body);  // Manage object detachment behavior
bool interaction_zone_is_attached(const DynamicInteractionZone* zone, const PhysicsBody* body);  // Check if a body is attached to the zone

// Environmental effects
void interaction_zone_apply_environment_effects(DynamicInteractionZone* zone, const EnvironmentForces* forces, f32 delta_time);  // Apply environmental effects (e.g., wind, gravity)
void interaction_zone_update_sticky_behavior(DynamicInteractionZone* zone, PhysicsBody* body, f32 delta_time);  // Refine stickiness over time
void interaction_zone_set_dynamic_expansion(DynamicInteractionZone* zone, f32 expansion_rate);  // Configure dynamic size changes over time

// Wind and breeze effects
void apply_wind_force_to_object(PhysicsBody* body, const EnvironmentForces* wind, f32 turbulence_factor);
void simulate_tree_bending(PhysicsBody* body, const EnvironmentForces* wind, f32 stiffness);
void wind_update_zone(EnvironmentForces* wind_zone, f32 delta_time);
void wind_zone_apply_to_particles(WindZone* wind_zone, ParticleSystem* particles);
void wind_zone_configure_turbulence(WindZone* wind_zone, f32 turbulence_factor);
void simulate_leaf_movement(PhysicsBody* leaf, const WindZone* wind_zone, f32 delta_time);

// Advanced wind zone behavior
void wind_zone_set_bounds(WindZone* wind_zone, const v3* bounds_min, const v3* bounds_max);  // Adjust wind zone bounds dynamically
void wind_zone_enable(WindZone* wind_zone, bool enabled);  // Enable or disable a wind zone
bool wind_zone_is_enabled(const WindZone* wind_zone);  // Check if a wind zone is active

// Object-specific wind interaction
void apply_wind_force_to_soft_body(SoftBody* body, const WindZone* wind_zone, f32 delta_time);  // Soft bodies react to wind
void apply_wind_force_to_ragdoll(Ragdoll* ragdoll, const WindZone* wind_zone, f32 delta_time);  // Simulate ragdoll interaction with wind
void apply_wind_to_water_body(WindZone* wind_zone, WaterBody* water_body, f32 delta_time);  // Wind effects on water surfaces

// Dynamic effects and turbulence
void wind_zone_add_dynamic_gust(WindZone* wind_zone, const v3* direction, f32 strength, f32 duration);  // Temporary gusts for dynamic wind effects
void wind_zone_update_gusts(WindZone* wind_zone, f32 delta_time);  // Update gust behavior over time

// Vegetation and detailed object simulation
void simulate_grass_waving(PhysicsBody* grass_patch, const WindZone* wind_zone, f32 delta_time, f32 stiffness);
void simulate_cloth_movement(SoftBody* cloth, const WindZone* wind_zone, f32 delta_time);
void simulate_tree_falling(PhysicsBody* tree, const WindZone* wind_zone, f32 threshold_force);  // Handle extreme wind events

// Constraint solvers
void solve_fixed_constraint(const Constraint* constraint);
void solve_hinge_constraint(const Constraint* constraint);
void solve_slider_constraint(const Constraint* constraint);
void solve_constraints(Constraint* constraints, u32 constraint_count, f32 delta_time);

// Advanced constraint solvers
void solve_spring_constraint(const Spring* spring, f32 delta_time);  // Solve spring constraints for soft connections
void solve_elastic_constraint(const ElasticConstraint* constraint, f32 delta_time);  // Solve elastic constraints with stretch limits
void solve_ball_and_socket_constraint(const Constraint* constraint);  // Simulate free rotational movement around a pivot point
void solve_distance_constraint(const Constraint* constraint);  // Maintain a fixed distance between two bodies

// Dynamic constraint adjustments
void constraint_set_stiffness(Constraint* constraint, f32 stiffness);  // Adjust stiffness dynamically
void constraint_set_damping(Constraint* constraint, f32 damping);  // Adjust damping dynamically
void constraint_enable(Constraint* constraint, bool enable);  // Enable or disable a constraint
bool constraint_is_enabled(const Constraint* constraint);  // Check if a constraint is active

// Particle and fluid simulation
void particle_system_update(ParticleSystem* system, f32 delta_time);
void fluid_simulation_update(Fluid* fluid, f32 delta_time);

// Particle system enhancements
void particle_system_initialize(ParticleSystem* system, u32 max_particles, f32 emission_rate);  // Initialize particle system
void particle_system_emit(ParticleSystem* system, u32 count, const v3* position, const v3* velocity);  // Emit particles at a position with a velocity
void particle_system_apply_force(ParticleSystem* system, const v3* force);  // Apply a global force (e.g., wind or gravity)
void particle_system_apply_drag(ParticleSystem* system, f32 drag_coefficient);  // Apply drag to particles
void particle_system_set_lifetime(ParticleSystem* system, f32 lifetime);  // Configure the lifetime of particles
void particle_system_set_color(ParticleSystem* system, const Color* start_color, const Color* end_color);  // Set start and end colors for particles

// Advanced fluid simulation
void fluid_initialize(Fluid* fluid, u32 particle_count, f32 rest_density, f32 viscosity);  // Initialize fluid simulation
void fluid_apply_force(Fluid* fluid, const v3* force);  // Apply a global force to the fluid
void fluid_add_emitter(Fluid* fluid, const v3* position, f32 rate);  // Add a fluid emitter
void fluid_simulate_interactions(Fluid* fluid, f32 delta_time);  // Simulate particle-to-particle interactions
void fluid_calculate_pressure(Fluid* fluid);  // Compute pressures for fluid particles
void fluid_apply_boundary_conditions(Fluid* fluid, const AABB* bounds);  // Constrain fluid within boundaries
void fluid_render(const Fluid* fluid);  // Render fluid particles or a smoothed surface

// Soft body dynamics
void soft_body_update(SoftBody* soft_body, f32 delta_time);
void soft_body_apply_force(SoftBody* soft_body, const v3* force);
void soft_body_solve_constraints(SoftBody* soft_body);
void soft_body_vertex_deformation(SoftBody* body, const v3* deformation_center, f32 influence_radius, f32 deformation_strength);
void soft_body_apply_global_force(SoftBody* body, const v3* force);
void soft_body_add_vertex(SoftBody* soft_body, const v3* position, f32 mass);
void soft_body_remove_vertex(SoftBody* soft_body, u32 vertex_index);
void soft_body_update_lattice(SoftBody* soft_body, const v3* deformation_center, f32 influence_radius, f32 lattice_strength);

// Soft body initialization and setup
void soft_body_initialize(SoftBody* soft_body, const v3* initial_positions, u32 vertex_count, f32 stiffness, f32 damping);  // Initialize a soft body with default properties
void soft_body_set_material_properties(SoftBody* soft_body, f32 stiffness, f32 damping);  // Adjust material properties dynamically

// Force and interaction extensions
void soft_body_apply_local_force(SoftBody* soft_body, u32 vertex_index, const v3* force);  // Apply force to a specific vertex
void soft_body_apply_gravity(SoftBody* soft_body, const v3* gravity);  // Apply gravity to the entire soft body
void soft_body_apply_pressure(SoftBody* soft_body, f32 pressure, const v3* direction);  // Simulate internal/external pressure forces

// Collision handling
void soft_body_handle_collision(SoftBody* soft_body, const Collider* collider);  // Resolve collisions with rigid or static bodies
void soft_body_self_collision(SoftBody* soft_body, f32 repulsion_force);  // Handle self-intersection with repulsion

// Topology manipulation
void soft_body_split(SoftBody* soft_body, u32 vertex_index_a, u32 vertex_index_b);  // Split the soft body between two vertices
void soft_body_merge(SoftBody* soft_body, const SoftBody* other_body, u32 contact_point);  // Merge two soft bodies

// Ragdoll simulation
void ragdoll_initialize(Ragdoll* ragdoll);
void ragdoll_update(Ragdoll* ragdoll, f32 delta_time);
void ragdoll_activate(Ragdoll* ragdoll);
void ragdoll_deactivate(Ragdoll* ragdoll);
void ragdoll_blend_to_animation(Ragdoll* ragdoll, f32 blend_factor);
void ragdoll_apply_joint_constraints(Ragdoll* ragdoll, const Constraint* constraints, u32 count);

// Ragdoll pose manipulation
void ragdoll_set_initial_pose(Ragdoll* ragdoll, const Transform* pose_array, u32 bone_count);  // Set an initial pose for the ragdoll
void ragdoll_save_current_pose(const Ragdoll* ragdoll, Transform* pose_array, u32* bone_count);  // Save the current ragdoll pose
void ragdoll_reset_pose(Ragdoll* ragdoll);  // Reset ragdoll to its default pose

// Advanced physics behavior
void ragdoll_apply_impulse(Ragdoll* ragdoll, const v3* impulse, const v3* contact_point);  // Apply an external impulse to the entire ragdoll
void ragdoll_apply_force(Ragdoll* ragdoll, const v3* force);  // Apply a global force to the ragdoll
void ragdoll_apply_damping(Ragdoll* ragdoll, f32 linear_damping, f32 angular_damping);  // Apply damping to stabilize movement

// Interaction and collision
void ragdoll_interact_with_environment(Ragdoll* ragdoll, const PhysicsZone* zone);  // Handle interaction with environmental zones (e.g., wind, water)
bool ragdoll_check_collision(const Ragdoll* ragdoll, const Collider* collider);  // Check for collisions with other objects
void ragdoll_resolve_self_collision(Ragdoll* ragdoll, f32 repulsion_force);  // Handle self-intersecting bones

// Gameplay-specific extensions
void ragdoll_make_ragdoll_dynamic(Ragdoll* ragdoll, f32 blend_duration);  // Gradually transition from animation to dynamic ragdoll
void ragdoll_make_ragdoll_static(Ragdoll* ragdoll);  // Lock ragdoll in place for gameplay events
void ragdoll_attach_to_object(Ragdoll* ragdoll, PhysicsBody* object, u32 bone_index);  // Attach a specific bone to another object
void ragdoll_detach_from_object(Ragdoll* ragdoll, u32 bone_index);  // Detach a bone from an object

// Spatial grid management
void spatial_grid_initialize(SpatialGrid* grid, const v3* grid_min, const v3* grid_max, const v3* cell_size);
void spatial_grid_insert(SpatialGrid* grid, PhysicsBody* body);
void spatial_grid_remove(SpatialGrid* grid, PhysicsBody* body);

// Dynamic grid updates
void spatial_grid_update(SpatialGrid* grid, PhysicsBody* body, const AABB* previous_bounds);  // Update grid placement after a body moves
void spatial_grid_resize(SpatialGrid* grid, const v3* new_grid_min, const v3* new_grid_max);  // Dynamically resize the grid bounds

// Advanced queries
u32 spatial_grid_query_nearest(const SpatialGrid* grid, const v3* position, PhysicsBody** results, u32 max_results, f32 max_distance);  // Find nearest bodies to a point
u32 spatial_grid_query_ray(const SpatialGrid* grid, const Ray* ray, PhysicsBody** results, u32 max_results);  // Raycast through the grid to find intersecting bodies

// Optimization utilities
void spatial_grid_optimize_cell(SpatialGrid* grid, u32 cell_index);  // Optimize a single cell by removing inactive or static bodies
void spatial_grid_rebuild(SpatialGrid* grid);  // Rebuild the grid structure to account for large changes
void spatial_grid_prune(SpatialGrid* grid);  // Remove inactive or outdated entries from the grid

// Advanced utilities
PhysicsRaycastHit physics_raycast(const v3* origin, const v3* direction, f32 max_distance);
bool raycast_collider(const Collider* collider, const Ray* ray, PhysicsRaycastHit* hit_data);

// Raycasting utilities
u32 physics_raycast_all(const v3* origin, const v3* direction, f32 max_distance, PhysicsRaycastHit* results, u32 max_results);  // Cast a ray and gather all hits along the path
u32 physics_raycast_closest(const v3* origin, const v3* direction, f32 max_distance, PhysicsRaycastHit* closest_hit);  // Cast a ray and return the closest hit
u32 physics_raycast_sorted(const v3* origin, const v3* direction, f32 max_distance, PhysicsRaycastHit* results, u32 max_results);  // Return hits sorted by distance

// Shape casting
bool physics_sphere_cast(const Sphere* sphere, const v3* direction, f32 max_distance, PhysicsRaycastHit* hit_data);  // Cast a sphere along a path
bool physics_aabb_cast(const AABB* aabb, const v3* direction, f32 max_distance, PhysicsRaycastHit* hit_data);  // Cast an AABB along a path
bool physics_obb_cast(const OBB* obb, const v3* direction, f32 max_distance, PhysicsRaycastHit* hit_data);  // Cast an OBB along a path

// Overlap testing
bool physics_overlap_sphere(const Sphere* sphere, PhysicsBody** results, u32 max_results);  // Find all bodies overlapping a sphere
bool physics_overlap_aabb(const AABB* aabb, PhysicsBody** results, u32 max_results);  // Find all bodies overlapping an AABB
bool physics_overlap_obb(const OBB* obb, PhysicsBody** results, u32 max_results);  // Find all bodies overlapping an OBB

// Line and shape intersection
bool physics_line_intersects_collider(const v3* point_a, const v3* point_b, const Collider* collider, PhysicsRaycastHit* hit_data);  // Check if a line segment intersects a collider
bool physics_sweep_test(const Collider* collider, const v3* direction, f32 max_distance, PhysicsRaycastHit* hit_data);  // Sweep a collider along a path and detect intersections

// Interaction Zones
void interaction_zone_initialize(InteractionZone* zone, const v3* center, f32 radius);
void interaction_zone_update(InteractionZone* zone, const PhysicsBody* body, f32 delta_time);
bool interaction_zone_check(const InteractionZone* zone, const PhysicsBody* body);
void interaction_zone_apply_force(const InteractionZone* zone, PhysicsBody* body);
void interaction_zone_trigger_events(const InteractionZone* zone, PhysicsBody* body, bool is_entering);
void interaction_zone_handle_overlap(InteractionZone* zone, PhysicsBody* body, bool is_entering);
void interaction_zone_apply_custom_behavior(InteractionZone* zone, PhysicsBody* body);

// Advanced Interaction Zone Features
void interaction_zone_set_behavior(InteractionZone* zone, void (*custom_behavior)(InteractionZone* zone, PhysicsBody* body)); // Set a custom behavior callback for the zone
void interaction_zone_set_activation_condition(InteractionZone* zone, bool (*condition)(const InteractionZone* zone, const PhysicsBody* body)); // Set a condition for zone activation
void interaction_zone_adjust_properties(InteractionZone* zone, const v3* new_center, f32 new_radius); // Dynamically adjust zone properties
void interaction_zone_apply_gravity_override(InteractionZone* zone, PhysicsBody* body, const v3* gravity_override); // Apply custom gravity within the zone
void interaction_zone_apply_friction_override(InteractionZone* zone, PhysicsBody* body, f32 friction_override); // Apply custom friction within the zone
void interaction_zone_handle_exit(InteractionZone* zone, PhysicsBody* body); // Handle behavior when a body exits the zone
void interaction_zone_check_proximity(InteractionZone* zone, const v3* position, f32 threshold, bool* is_within); // Check if a position is within a certain threshold of the zone
void interaction_zone_enable_dynamic_resizing(InteractionZone* zone, bool enable, f32 resize_rate); // Enable dynamic resizing of the zone
void interaction_zone_apply_force_pattern(InteractionZone* zone, PhysicsBody* body, const v3* force_pattern); // Apply a specific force pattern to bodies in the zone
void interaction_zone_link_with_signal(InteractionZone* zone, const char* signal_name); // Link the interaction zone with a specific signal

// Destructible objects
void destruct_object(const DestructibleObject* object, PhysicsBody** fragments, u32 fragment_count);
void destructible_object_shatter(DestructibleObject* object, PhysicsBody** fragments, u32 fragment_count, const v3* explosion_center, f32 explosion_force);
void destructible_object_apply_damage(DestructibleObject* object, f32 damage_amount);
void destructible_object_repair(DestructibleObject* object, f32 repair_amount);

// Advanced Destructible Object Features
void destructible_object_check_integrity(const DestructibleObject* object, bool* is_intact); // Check if the object is still intact
void destructible_object_set_health(DestructibleObject* object, f32 max_health, f32 initial_health); // Configure health settings for the object
void destructible_object_on_destroy(DestructibleObject* object, void (*callback)(const DestructibleObject* object)); // Set a callback for when the object is destroyed
void destructible_object_apply_force_to_fragments(DestructibleObject* object, PhysicsBody** fragments, u32 fragment_count, const v3* force); // Apply force to all fragments after destruction
void destructible_object_reassemble(DestructibleObject* object, PhysicsBody** fragments, u32 fragment_count); // Reassemble the object from fragments
void destructible_object_set_explosion_properties(DestructibleObject* object, f32 blast_radius, f32 force_multiplier); // Configure explosion behavior
void destructible_object_enable_dynamic_fragments(DestructibleObject* object, bool enable); // Enable/disable dynamic physics for fragments
void destructible_object_handle_environment_effects(DestructibleObject* object, const PhysicsZone* zone); // Handle environmental effects like fire or water

// Force and impulse amplification
void physics_body_apply_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor);
void physics_body_apply_amplified_impulse(RigidBody* body, const v3* impulse, const v3* contact_point, f32 amplification_factor);

// Advanced Force and Impulse Amplification Features
void physics_body_apply_directional_force(RigidBody* body, const v3* force, const v3* direction, f32 amplification_factor); // Apply force in a specific direction with amplification
void physics_body_apply_timed_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor, f32 duration); // Apply amplified force over a specific duration
void physics_body_apply_area_impulse(PhysicsBody* body, const v3* origin, f32 radius, f32 amplification_factor); // Apply an amplified impulse affecting objects within a radius
void physics_body_set_amplification_limits(RigidBody* body, f32 max_force, f32 max_impulse); // Configure maximum amplification limits for forces and impulses
void physics_body_apply_decay_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor, f32 decay_rate); // Apply force with amplification that decays over time

// Advanced math utilities
quat integrate_angular_velocity(const quat* orientation, const v3* angular_velocity, f32 delta_time);

// Friction and rolling resistance
void apply_dynamic_friction(RigidBody* body, f32 dynamic_coefficient);
void apply_rolling_resistance(RigidBody* body, f32 rolling_coefficient);

// Friction Enhancements
void apply_static_friction(RigidBody* body, const v3* external_force, f32 static_coefficient); // Apply static friction to resist initial motion
void adjust_friction_for_surface(RigidBody* body, const PhysicsZone* surface_zone); // Modify friction based on surface properties dynamically

// Rolling Resistance Enhancements
void apply_rolling_resistance_with_velocity(RigidBody* body, f32 rolling_coefficient, const v3* velocity); // Account for velocity in rolling resistance
void apply_torque_based_rolling_resistance(RigidBody* body, f32 rolling_coefficient); // Apply rolling resistance as a torque opposing rotation

// Specialized Friction Models
void apply_combined_friction(RigidBody* body, f32 dynamic_coefficient, f32 static_coefficient); // Combine dynamic and static friction forces
void simulate_anisotropic_friction(RigidBody* body, const v3* friction_direction, f32 coefficient); // Simulate directional friction (e.g., sliding along grooves)

// Environmental noise
v3 generate_wind_force(const EnvironmentForces* environment, const v3* position, f32 time);

// Advanced Environmental Noise
v3 generate_turbulent_wind_force(const EnvironmentForces* environment, const v3* position, f32 time, f32 turbulence_strength); 
// Generate wind force with turbulence effects

void apply_noise_to_force(v3* force, f32 noise_intensity, f32 time); 
// Add random noise to an existing force vector for unpredictable behavior

void simulate_perlin_wind_force(const EnvironmentForces* environment, const v3* position, f32 time, f32 scale, f32 amplitude, v3* output_force); 
// Simulate smooth, natural wind variations using Perlin noise

void update_environmental_noise(EnvironmentForces* environment, f32 delta_time); 
// Update environmental forces dynamically based on time or player actions

// Soft body deformation
void soft_body_vertex_deformation(SoftBody* body, const v3* deformation_center, f32 influence_radius, f32 deformation_strength);
void soft_body_apply_global_force(SoftBody* body, const v3* force);

// Advanced Soft Body Deformation
void soft_body_simulate_internal_pressure(SoftBody* body, f32 pressure_amount);
// Simulate internal pressure forces, such as inflating or deflating soft bodies.

void soft_body_apply_surface_tension(SoftBody* body, f32 tension_factor);
// Apply surface tension forces to maintain shape consistency under deformation.

void soft_body_add_anchor_point(SoftBody* body, u32 vertex_index, const v3* anchor_position);
// Attach a vertex of the soft body to a fixed or moving anchor point.

void soft_body_remove_anchor_point(SoftBody* body, u32 vertex_index);
// Detach a previously anchored vertex from its fixed position.

void soft_body_update_mass_distribution(SoftBody* body, f32 total_mass);
// Recalculate mass distribution across vertices after deformation or changes.

void soft_body_compute_stress_map(SoftBody* body, f32* stress_values);
// Compute a stress map indicating tension levels at each vertex.

void soft_body_apply_external_collision(SoftBody* body, const Collision* collision);
// Handle deformation caused by collisions with other objects.

// Springs
void spring_initialize(Spring* spring, PhysicsBody* body_a, PhysicsBody* body_b, f32 rest_length, f32 stiffness, f32 damping);
void spring_update(Spring* spring, f32 delta_time);
void spring_apply_forces(Spring* spring);

// Advanced Spring Mechanics
void spring_set_length_limits(Spring* spring, f32 min_length, f32 max_length);
// Set minimum and maximum length constraints for the spring.

void spring_break_if_overstretched(Spring* spring, f32 stretch_threshold);
// Break the spring if it is stretched beyond the specified threshold.

void spring_apply_external_force(Spring* spring, const v3* force);
// Apply an external force directly to the spring's endpoints.

void spring_simulate_oscillation(Spring* spring, f32 frequency, f32 amplitude);
// Simulate harmonic oscillations for a spring system.

void spring_reconnect(Spring* spring, PhysicsBody* new_body_a, PhysicsBody* new_body_b);
// Reconnect the spring to new physics bodies, useful for dynamic systems.

void spring_visualize(const Spring* spring, void* renderer_context);
// Render the spring for debugging or visualization purposes.

void spring_adjust_damping(Spring* spring, f32 new_damping_value);
// Dynamically adjust the damping coefficient of the spring.

void spring_compute_tension(const Spring* spring, f32* out_tension_value);
// Compute the current tension force exerted by the spring.

// Waterbody
void water_body_initialize(WaterBody* water_body, const v3* position, const v2* size, f32 wave_height, f32 wave_speed, f32 buoyancy);
void water_body_update(WaterBody* water_body, f32 delta_time);
void water_body_apply_buoyancy(WaterBody* water_body, PhysicsBody* body);
void water_body_simulate_waves(WaterBody* water_body, f32 delta_time);

// Advanced Waterbody Mechanics
void water_body_apply_drag(WaterBody* water_body, PhysicsBody* body, f32 drag_coefficient);
// Apply drag forces to objects moving through the water.

void water_body_generate_ripples(WaterBody* water_body, const v3* impact_point, f32 ripple_intensity, f32 ripple_decay);
// Generate ripple effects at a specific impact point, simulating disturbances in the water.

void water_body_check_submersion(const WaterBody* water_body, const PhysicsBody* body, bool* is_submerged, f32* submersion_depth);
// Check if a physics body is submerged in the water and compute the depth of submersion.

void water_body_interact_with_particles(WaterBody* water_body, ParticleSystem* particles, f32 delta_time);
// Simulate interactions between the waterbody and particle systems, such as splashes or mist.

void water_body_adjust_wave_properties(WaterBody* water_body, f32 new_wave_height, f32 new_wave_speed);
// Dynamically adjust the wave properties of the waterbody for real-time effects.

void water_body_simulate_current(WaterBody* water_body, const v3* flow_direction, f32 flow_speed);
// Simulate directional water currents affecting all objects within the waterbody.

void water_body_render_debug(const WaterBody* water_body, void* renderer_context);
// Visualize the waterbody for debugging purposes, including wave patterns and buoyancy zones.

void water_body_apply_temperature_effects(WaterBody* water_body, PhysicsBody* body, f32 temperature);
// Simulate temperature-related effects, such as freezing or reduced buoyancy in colder water.

void water_body_interact_with_wind(WaterBody* water_body, const WindZone* wind_zone);
// Simulate how wind interacts with the water surface, creating dynamic wave patterns.


// Windzone
void wind_zone_initialize(WindZone* wind_zone, const v3* direction, f32 strength, const v3* bounds_min, const v3* bounds_max);
void wind_zone_apply_force(WindZone* wind_zone, PhysicsBody* body);
void wind_zone_update(WindZone* wind_zone, f32 delta_time);

// Advanced Windzone Mechanics
void wind_zone_adjust_strength(WindZone* wind_zone, f32 new_strength);
// Dynamically adjust the wind strength for environmental changes.

void wind_zone_apply_turbulence(WindZone* wind_zone, PhysicsBody* body, f32 turbulence_intensity, f32 noise_factor);
// Apply turbulence effects to physics bodies for chaotic wind interactions.

void wind_zone_generate_gusts(WindZone* wind_zone, f32 gust_intensity, f32 gust_frequency, f32 delta_time);
// Simulate periodic wind gusts for more dynamic wind behavior.

void wind_zone_interact_with_particles(WindZone* wind_zone, ParticleSystem* particles, f32 delta_time);
// Simulate wind effects on particle systems, such as scattering or directional flow.

void wind_zone_configure_bounds(WindZone* wind_zone, const v3* new_bounds_min, const v3* new_bounds_max);
// Reconfigure the effective area of the wind zone for real-time adjustments.

void wind_zone_render_debug(const WindZone* wind_zone, void* renderer_context);
// Visualize the wind zone for debugging purposes, including bounds and direction.

void wind_zone_simulate_interaction_with_trees(WindZone* wind_zone, PhysicsBody** trees, u32 tree_count, f32 stiffness, f32 damping);
// Apply wind effects specifically tailored for tree bending and movement.

void wind_zone_apply_to_soft_bodies(WindZone* wind_zone, SoftBody* soft_body, f32 delta_time);
// Simulate wind effects on soft bodies, such as cloth or deformable objects.

void wind_zone_apply_to_water(WindZone* wind_zone, WaterBody* water_body, f32 wave_amplification_factor);
// Integrate wind effects with water bodies, creating dynamic wave patterns.

void wind_zone_toggle(WindZone* wind_zone, bool is_active);
// Enable or disable the wind zone dynamically.

// Spatial grid
void spatial_grid_initialize(SpatialGrid* grid, const v3* grid_min, const v3* grid_max, const v3* cell_size);
void spatial_grid_insert(SpatialGrid* grid, PhysicsBody* body);
void spatial_grid_remove(SpatialGrid* grid, PhysicsBody* body);
void spatial_grid_query(const SpatialGrid* grid, const AABB* bounds, PhysicsBody** results, u32* result_count);

// Spatial Grid Utilities
void spatial_grid_optimize(SpatialGrid* grid);
// Optimize the spatial grid layout for better performance (e.g., merge or split cells based on density).

void spatial_grid_clear(SpatialGrid* grid);
// Remove all entries from the spatial grid without deallocating memory.

bool spatial_grid_is_empty(const SpatialGrid* grid);
// Check if the spatial grid is empty.

void spatial_grid_get_cell(const SpatialGrid* grid, const v3* position, SpatialGridCell* out_cell);
// Retrieve the cell containing a specific position.

u32 spatial_grid_get_neighbors(const SpatialGrid* grid, const SpatialGridCell* cell, SpatialGridCell** neighbors, u32 max_neighbors);
// Retrieve neighboring cells for a given cell.

void spatial_grid_render_debug(const SpatialGrid* grid, void* renderer_context);
// Render the spatial grid structure for debugging purposes.

void spatial_grid_query_closest(const SpatialGrid* grid, const v3* position, f32 max_distance, PhysicsBody** results, u32* result_count);
// Query the grid for bodies closest to a given position within a maximum distance.
void spatial_grid_set_body_layer(SpatialGrid* grid, PhysicsBody* body, u32 layer_mask);
// Set the layer mask for a specific body to control which layers it interacts with.

// Broadphase cache
void broadphase_pair_cache_initialize(BroadphasePairCache* cache, u32 initial_capacity);
void broadphase_pair_cache_add(BroadphasePairCache* cache, PhysicsBody* body_a, PhysicsBody* body_b);
void broadphase_pair_cache_clear(BroadphasePairCache* cache);

// Broadphase Pair Cache Utilities
bool broadphase_pair_cache_contains(const BroadphasePairCache* cache, PhysicsBody* body_a, PhysicsBody* body_b);
// Check if a specific pair of bodies exists in the cache.

void broadphase_pair_cache_remove(BroadphasePairCache* cache, PhysicsBody* body_a, PhysicsBody* body_b);
// Remove a specific pair of bodies from the cache.

u32 broadphase_pair_cache_query(const BroadphasePairCache* cache, PhysicsBody* body, PhysicsBody** results, u32 max_results);
// Query the cache for all pairs involving a specific body.

void broadphase_pair_cache_resize(BroadphasePairCache* cache, u32 new_capacity);
// Resize the cache to accommodate more pairs.

void broadphase_pair_cache_iterate(const BroadphasePairCache* cache, void (*callback)(PhysicsBody* body_a, PhysicsBody* body_b, void* user_data), void* user_data);
// Iterate through all pairs in the cache with a callback function.

u32 broadphase_pair_cache_get_count(const BroadphasePairCache* cache);
// Retrieve the number of pairs currently in the cache.

void broadphase_pair_cache_render_debug(const BroadphasePairCache* cache, void* renderer_context);
// Render debug information about the cache, such as pair connections.

// Sticker
void sticker_initialize(Sticker* sticker, const v3* position, const v3* normal, u32 texture_id, f32 scale, Color tint);
void sticker_apply_to_surface(Sticker* sticker, SurfaceProperties* surface);
void sticker_remove(Sticker* sticker);
void sticker_update(Sticker* sticker, f32 delta_time);

// Sticker management
void sticker_set_tint(Sticker* sticker, Color tint);
// Update the tint of a sticker dynamically.

void sticker_scale(Sticker* sticker, f32 scale_factor);
// Adjust the scale of an existing sticker.

bool sticker_is_on_surface(const Sticker* sticker, const SurfaceProperties* surface);
// Check if a sticker is still properly attached to a given surface.

void sticker_attach_to_dynamic_surface(Sticker* sticker, PhysicsBody* dynamic_surface);
// Attach a sticker to a dynamic surface, updating its position and orientation as the surface moves.

void sticker_apply_deformation(Sticker* sticker, const v3* deformation_vector);
// Apply a deformation effect to a sticker, such as stretching or bending.

void sticker_interact_with_physics(Sticker* sticker, PhysicsBody* body, f32 interaction_strength);
// Enable interactions between stickers and nearby physical objects (e.g., "stickers fluttering off" during movement).

void sticker_render_debug(const Sticker* sticker, void* renderer_context);
// Render debug information for the sticker, such as its bounds or attachment points.

// Elastic constraint
void elastic_constraint_initialize(ElasticConstraint* constraint, PhysicsBody* body_a, PhysicsBody* body_b, f32 rest_length, f32 stiffness, f32 damping, f32 max_stretch_length);
void elastic_constraint_update(ElasticConstraint* constraint, f32 delta_time);
void elastic_constraint_apply_forces(ElasticConstraint* constraint);

// Elastic constraint management
void elastic_constraint_set_stiffness(ElasticConstraint* constraint, f32 stiffness);
// Dynamically adjust the stiffness of an elastic constraint.

void elastic_constraint_set_damping(ElasticConstraint* constraint, f32 damping);
// Dynamically adjust the damping factor to control oscillations.

bool elastic_constraint_is_within_limits(const ElasticConstraint* constraint);
// Check if the current stretch of the constraint is within the defined maximum stretch length.

void elastic_constraint_break(ElasticConstraint* constraint);
// Break the constraint when its maximum stretch length is exceeded.

void elastic_constraint_visualize(const ElasticConstraint* constraint, void* renderer_context);
// Visualize the elastic constraint for debugging or gameplay effects.

// Detatchment rules
void detachment_rules_apply(const DetachmentRules* rules, PhysicsBody* body, f32 current_force, f32 current_torque);
bool detachment_rules_can_reattach(const DetachmentRules* rules);

// Detachment rules management
void detachment_rules_set_max_force(DetachmentRules* rules, f32 max_force);
// Dynamically adjust the maximum force threshold for detachment.

void detachment_rules_set_max_torque(DetachmentRules* rules, f32 max_torque);
// Dynamically adjust the maximum torque threshold for detachment.

void detachment_rules_reset(DetachmentRules* rules);
// Reset the detachment rules to their default or initial state.

void detachment_rules_visualize(const DetachmentRules* rules, void* renderer_context);
// Visualize the detachment thresholds and status for debugging or gameplay effects.

bool detachment_rules_is_near_breaking(const DetachmentRules* rules, f32 current_force, f32 current_torque);
// Check if the current force or torque is approaching the detachment thresholds.

// Fluid interaction
void fluid_interaction_apply(FluidInteraction* interaction, PhysicsBody* body);
void fluid_interaction_update(FluidInteraction* interaction, f32 delta_time);

// Enhanced fluid interaction
void fluid_interaction_calculate_drag(FluidInteraction* interaction, PhysicsBody* body, f32* out_drag_force);
// Compute the drag force applied to a body based on its velocity and the fluid's drag coefficient.

void fluid_interaction_calculate_buoyancy(FluidInteraction* interaction, PhysicsBody* body, f32* out_buoyancy_force);
// Compute the buoyancy force applied to a body based on its submerged volume and the fluid's buoyancy factor.

void fluid_interaction_set_properties(FluidInteraction* interaction, f32 buoyancy_factor, f32 drag_coefficient);
// Adjust the fluid interaction properties dynamically during runtime.

bool fluid_interaction_is_submerged(const FluidInteraction* interaction, const PhysicsBody* body);
// Check if a body is fully or partially submerged in the fluid.

void fluid_interaction_apply_turbulence(FluidInteraction* interaction, PhysicsBody* body, const v3* turbulence_vector);
// Apply turbulence forces to the body for dynamic and chaotic fluid behavior.

void fluid_interaction_visualize(const FluidInteraction* interaction, void* renderer_context);
// Visualize fluid interaction effects for debugging or visual feedback in gameplay.

// Decal physics
void decal_physics_initialize(DecalPhysics* decal, const Sticker* base_sticker, const v3* deformation_vector, bool scales_with_surface);
void decal_physics_apply(DecalPhysics* decal, SurfaceProperties* surface);
void decal_physics_update(DecalPhysics* decal, f32 delta_time);
void decal_physics_remove(DecalPhysics* decal);

// Advanced decal physics
void decal_physics_set_properties(DecalPhysics* decal, f32 elasticity, f32 friction);
// Set the material properties for the decal to affect how it interacts with surfaces.

void decal_physics_apply_dynamic_forces(DecalPhysics* decal, const v3* force, const v3* point_of_application);
// Apply dynamic forces to decals, useful for interactive decals like patches on deformable surfaces.

bool decal_physics_check_collision(const DecalPhysics* decal, const SurfaceProperties* surface);
// Check if the decal is in contact with a surface, enabling dynamic behavior like peeling.

void decal_physics_attach_to_body(DecalPhysics* decal, PhysicsBody* body);
// Attach a decal to a physics body for dynamic movement and deformation.

void decal_physics_adjust_for_surface_deformation(DecalPhysics* decal, const v3* deformation_vector);
// Adjust the decal's position and orientation to align with deformed or animated surfaces.

void decal_physics_render_debug(const DecalPhysics* decal, void* renderer_context);
// Render debugging visuals for decal physics interactions.