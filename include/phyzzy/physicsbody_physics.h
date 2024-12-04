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
#include "physics.h"

#ifdef __cplusplus
extern "C" {
#endif

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

// Wobble and oscillation effects
void physics_body_apply_directional_wobble(PhysicsBody* body, const v3* direction, f32 amplitude, f32 frequency);  // Directional wobble effects
void physics_body_sync_wobble_with_animation(PhysicsBody* body, f32 animation_phase);  // Synchronize wobble with an external animation system
void physics_body_apply_impact_wobble(PhysicsBody* body, const v3* impact_force, f32 decay_rate);  // Wobble caused by external impacts
void physics_body_get_wobble_state(const PhysicsBody* body, f32* amplitude_out, f32* frequency_out);  // Retrieve current wobble properties
void physics_body_toggle_wobble(PhysicsBody* body, bool enable);  // Enable or disable wobble dynamically

// Dynamic floor effects
void physics_body_update_surface_interaction(PhysicsBody* body, PhysicsZone* floor_zone);
void apply_friction_for_surface(PhysicsBody* body, const PhysicsZone* floor_zone);
void apply_surface_restitution(PhysicsBody* body, const PhysicsZone* floor_zone);
void physics_body_adjust_for_wet_surface(PhysicsBody* body, PhysicsZone* floor_zone);
void physics_body_apply_ice_friction(PhysicsBody* body, PhysicsZone* floor_zone);

// Advanced floor interaction
void physics_body_adjust_for_slippery_surface(PhysicsBody* body, const PhysicsZone* floor_zone);  // Adjust dynamics for slippery surfaces like polished floors
void physics_body_interact_with_sticky_surface(PhysicsBody* body, const PhysicsZone* floor_zone);  // Simulate stickiness interaction (e.g., glue or tar)
void physics_body_simulate_rough_surface(PhysicsBody* body, const PhysicsZone* floor_zone, f32 roughness_factor);  // Add resistance for rough terrains

// Force and impulse amplification
void physics_body_apply_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor);
void physics_body_apply_amplified_impulse(RigidBody* body, const v3* impulse, const v3* contact_point, f32 amplification_factor);

// Combined floor effects
void physics_body_apply_combined_surface_effects(PhysicsBody* body, PhysicsZone* floor_zone);  // Apply all relevant surface effects simultaneously
void physics_body_check_surface_transition(PhysicsBody* body, const PhysicsZone* previous_zone, const PhysicsZone* current_zone);  // Handle transitions between different surface types

// Advanced Force and Impulse Amplification Features
void physics_body_apply_directional_force(RigidBody* body, const v3* force, const v3* direction, f32 amplification_factor); // Apply force in a specific direction with amplification
void physics_body_apply_timed_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor, f32 duration); // Apply amplified force over a specific duration
void physics_body_apply_area_impulse(PhysicsBody* body, const v3* origin, f32 radius, f32 amplification_factor); // Apply an amplified impulse affecting objects within a radius
void physics_body_set_amplification_limits(RigidBody* body, f32 max_force, f32 max_impulse); // Configure maximum amplification limits for forces and impulses
void physics_body_apply_decay_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor, f32 decay_rate); // Apply force with amplification that decays over time

void physics_body_update_timed_forces(RigidBody* body, f32 delta_time);
void physics_body_update_decay_forces(RigidBody* body, f32 delta_time);

#ifdef __cplusplus
}
#endif
