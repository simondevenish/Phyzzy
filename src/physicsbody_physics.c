/*
 * Phyzzy: Playful Physics for Imaginative Games
 *
 * Licensed under the GNU General Public License, Version 3.
 * For license details, visit: https://www.gnu.org/licenses/gpl-3.0.html
 *
 * Questions or contributions? Reach out to Simon Devenish:
 * simon.devenish@outlook.com
 */

#include "physicsbody_physics.h"
#include "vector_math.h"
#include <stdlib.h>
#include <string.h> // For memset

// Properly cleanup and deallocate a physics body
void physics_body_destroy(PhysicsBody* body) {
    if (body == NULL) {
        return;
    }

    // Destroy the rigid body if it exists
    if (body->rigid_body != NULL) {
        free(body->rigid_body);
        body->rigid_body = NULL;
    }

    // Destroy the collider if it exists
    if (body->collider != NULL) {
        free(body->collider);
        body->collider = NULL;
    }

    // Free the physics body itself
    free(body);
}

// Dynamically change the mass of a body
void physics_body_set_mass(RigidBody* body, f32 mass) {
    if (body == NULL || mass <= 0.0f) {
        return;
    }

    body->mass = mass;
    body->inverse_mass = 1.0f / mass;

    // Recalculate inertia tensor and its inverse if necessary
    // Placeholder for inertia tensor calculation
    // body->inertia_tensor = compute_inertia_tensor(body);
    // body->inverse_inertia_tensor = invert_matrix(body->inertia_tensor);
}

// Retrieve the current mass of a body
void physics_body_get_mass(const RigidBody* body, f32* mass) {
    if (body == NULL || mass == NULL) {
        return;
    }

    *mass = body->mass;
}

// Set the linear velocity
void physics_body_set_velocity(RigidBody* body, const v3* velocity) {
    if (body == NULL || velocity == NULL) {
        return;
    }

    body->velocity = *velocity;
}

// Retrieve the linear velocity
void physics_body_get_velocity(const RigidBody* body, v3* velocity) {
    if (body == NULL || velocity == NULL) {
        return;
    }

    *velocity = body->velocity;
}

// Set rotational velocity
void physics_body_set_angular_velocity(RigidBody* body, const v3* angular_velocity) {
    if (body == NULL || angular_velocity == NULL) {
        return;
    }

    body->angular_velocity = *angular_velocity;
}

// Retrieve rotational velocity
void physics_body_get_angular_velocity(const RigidBody* body, v3* angular_velocity) {
    if (body == NULL || angular_velocity == NULL) {
        return;
    }

    *angular_velocity = body->angular_velocity;
}

// Reset all accumulated forces and torques
void physics_body_reset_forces(RigidBody* body) {
    if (body == NULL) {
        return;
    }

    body->force = (v3){0.0f, 0.0f, 0.0f};
    body->torque = (v3){0.0f, 0.0f, 0.0f};
}

// Move a physics body directly
void physics_body_translate(PhysicsBody* body, const v3* translation) {
    if (body == NULL || translation == NULL) {
        return;
    }

    // Update the body's transform
    body->transform.position = v3_add(body->transform.position, *translation);

    // If rigid body exists, update its position
    if (body->rigid_body != NULL) {
        body->rigid_body->position = body->transform.position;
    }

    // If collider exists, update its transform
    if (body->collider != NULL) {
        body->collider->transform.position = body->transform.position;
    }
}

// Rotate a physics body directly
void physics_body_rotate(PhysicsBody* body, const quat* rotation) {
    if (body == NULL || rotation == NULL) {
        return;
    }

    // Update the body's transform
    body->transform.rotation = quat_multiply(*rotation, body->transform.rotation);

    // If rigid body exists, update its orientation
    if (body->rigid_body != NULL) {
        body->rigid_body->orientation = body->transform.rotation;
    }

    // If collider exists, update its transform
    if (body->collider != NULL) {
        body->collider->transform.rotation = body->transform.rotation;
    }
}

// Explicitly set the position
void physics_body_set_position(PhysicsBody* body, const v3* position) {
    if (body == NULL || position == NULL) {
        return;
    }

    // Update the body's transform
    body->transform.position = *position;

    // If rigid body exists, update its position
    if (body->rigid_body != NULL) {
        body->rigid_body->position = *position;
    }

    // If collider exists, update its transform
    if (body->collider != NULL) {
        body->collider->transform.position = *position;
    }
}

// Retrieve the current position
void physics_body_get_position(const PhysicsBody* body, v3* position) {
    if (body == NULL || position == NULL) {
        return;
    }

    *position = body->transform.position;
}

// Dynamically toggle static/dynamic state
void physics_body_set_static(PhysicsBody* body, bool is_static) {
    if (body == NULL) {
        return;
    }

    body->is_static = is_static;

    if (body->rigid_body != NULL) {
        if (is_static) {
            body->rigid_body->inverse_mass = 0.0f;
            body->rigid_body->mass = 0.0f;
            body->rigid_body->velocity = (v3){0.0f, 0.0f, 0.0f};
            body->rigid_body->angular_velocity = (v3){0.0f, 0.0f, 0.0f};
            // Zero out inverse inertia tensor
            memset(&body->rigid_body->inverse_inertia_tensor, 0, sizeof(m3));
        } else {
            // Set default mass if mass is zero
            if (body->rigid_body->mass == 0.0f) {
                body->rigid_body->mass = 1.0f; // Default mass
            }
            body->rigid_body->inverse_mass = 1.0f / body->rigid_body->mass;
            // Recalculate inverse inertia tensor
            // body->rigid_body->inverse_inertia_tensor = invert_matrix(body->rigid_body->inertia_tensor);
        }
    }
}

// Check if a body is currently static
bool physics_body_is_static(const PhysicsBody* body) {
    if (body == NULL) {
        return false;
    }

    return body->is_static;
}

// Apply damping to reduce velocity over time
void physics_body_apply_damping(RigidBody* body, f32 linear_damping, f32 angular_damping) {
    if (body == NULL) {
        return;
    }

    // Ensure damping factors are between 0 and 1
    linear_damping = clamp(linear_damping, 0.0f, 1.0f);
    angular_damping = clamp(angular_damping, 0.0f, 1.0f);

    body->velocity = v3_scale(body->velocity, 1.0f - linear_damping);
    body->angular_velocity = v3_scale(body->angular_velocity, 1.0f - angular_damping);
}

// Apply drag proportional to velocity
void physics_body_apply_drag(RigidBody* body, f32 drag_coefficient) {
    if (body == NULL) {
        return;
    }

    // Calculate drag force: F_drag = -c * v
    v3 drag_force = v3_scale(body->velocity, -drag_coefficient);

    // Accumulate the drag force
    body->force = v3_add(body->force, drag_force);
}

// Freeze a physics body to stop updates
void physics_body_freeze(PhysicsBody* body) {
    if (body == NULL) {
        return;
    }

    if (body->rigid_body != NULL) {
        body->rigid_body->velocity = (v3){0.0f, 0.0f, 0.0f};
        body->rigid_body->angular_velocity = (v3){0.0f, 0.0f, 0.0f};
        body->rigid_body->inverse_mass = 0.0f;
        memset(&body->rigid_body->inverse_inertia_tensor, 0, sizeof(m3));
    }

    body->is_static = true;
}

// Unfreeze a physics body to resume updates
void physics_body_unfreeze(PhysicsBody* body) {
    if (body == NULL) {
        return;
    }

    if (body->rigid_body != NULL) {
        if (body->rigid_body->mass == 0.0f) {
            body->rigid_body->mass = 1.0f; // Default mass
        }
        body->rigid_body->inverse_mass = 1.0f / body->rigid_body->mass;
        // Recalculate inverse inertia tensor if necessary
        // body->rigid_body->inverse_inertia_tensor = invert_matrix(body->rigid_body->inertia_tensor);
    }

    body->is_static = false;
}

void physics_body_apply_directional_wobble(PhysicsBody* body, const v3* direction, f32 amplitude, f32 frequency) {
    if (body == NULL || direction == NULL) {
        return;
    }

    body->wobble.enabled = true;
    body->wobble.amplitude = amplitude;
    body->wobble.frequency = frequency;
    body->wobble.direction = v3_normalize(*direction);
    body->wobble.time = 0.0f; // Reset time tracker
}

void physics_body_sync_wobble_with_animation(PhysicsBody* body, f32 animation_phase) {
    if (body == NULL) {
        return;
    }

    // Assume animation_phase is between 0 and 1
    body->wobble.time = animation_phase / body->wobble.frequency;
}

void physics_body_apply_impact_wobble(PhysicsBody* body, const v3* impact_force, f32 decay_rate) {
    if (body == NULL || impact_force == NULL) {
        return;
    }

    body->wobble.enabled = true;
    body->wobble.direction = v3_normalize(*impact_force);
    body->wobble.amplitude = v3_magnitude(*impact_force) * decay_rate;
    body->wobble.frequency = decay_rate; // You can adjust this as needed
    body->wobble.time = 0.0f;
}

void physics_body_get_wobble_state(const PhysicsBody* body, f32* amplitude_out, f32* frequency_out) {
    if (body == NULL) {
        return;
    }

    if (amplitude_out != NULL) {
        *amplitude_out = body->wobble.amplitude;
    }

    if (frequency_out != NULL) {
        *frequency_out = body->wobble.frequency;
    }
}

void physics_body_toggle_wobble(PhysicsBody* body, bool enable) {
    if (body == NULL) {
        return;
    }

    body->wobble.enabled = enable;

    if (!enable) {
        // Reset wobble parameters
        body->wobble.amplitude = 0.0f;
        body->wobble.frequency = 0.0f;
        body->wobble.time = 0.0f;
    }
}

void physics_body_update_wobble(PhysicsBody* body, f32 delta_time) {
    if (body == NULL || !body->wobble.enabled) {
        return;
    }

    // Update time
    body->wobble.time += delta_time;

    // Calculate wobble offset
    f32 wobble_angle = body->wobble.amplitude * sinf(body->wobble.frequency * body->wobble.time);

    // Create a rotation quaternion around the wobble direction
    quat wobble_rotation = quat_from_axis_angle(&body->wobble.direction, wobble_angle);

    // Apply the wobble rotation to the body's orientation
    body->transform.rotation = quat_multiply(wobble_rotation, body->transform.rotation);

    // Update the rigid body's orientation
    if (body->rigid_body != NULL) {
        body->rigid_body->orientation = body->transform.rotation;
    }

    // Update the collider's rotation if needed
    if (body->collider != NULL) {
        body->collider->transform.rotation = body->transform.rotation;
    }

    // Decrease amplitude over tim to decay the wobble
    // body->wobble.amplitude *= damping_factor; // where damping_factor < 1.0f
}

void apply_friction_for_surface(PhysicsBody* body, const PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL) {
        return;
    }

    // Get the friction coefficients
    f32 kinetic_friction_coefficient = floor_zone->friction_override;
    f32 static_friction_coefficient = floor_zone->static_friction_coefficient; // Ensure this field exists

    // Calculate the normal force (assuming N = mass * gravity magnitude)
    f32 gravity_magnitude = v3_magnitude(physics_state.settings.gravity);
    f32 normal_force = body->rigid_body->mass * gravity_magnitude;

    // Get the body's current velocity
    v3 velocity = body->rigid_body->velocity;
    f32 speed = v3_magnitude(velocity);

    // Small epsilon to avoid division by zero
    const f32 epsilon = 1e-5f;

    if (speed > epsilon) {
        // Kinetic friction
        // Friction force magnitude
        f32 friction_magnitude = kinetic_friction_coefficient * normal_force;

        // Normalize the velocity vector
        v3 velocity_normalized = v3_scale(velocity, 1.0f / speed);

        // Calculate the friction force
        v3 friction_force = v3_scale(velocity_normalized, -friction_magnitude);

        // Accumulate the friction force
        body->rigid_body->force = v3_add(body->rigid_body->force, friction_force);
    } else {
        // Static friction
        // Get the total external force (excluding friction)
        v3 total_external_force = body->rigid_body->force;

        // Compute the maximum static friction force
        f32 max_static_friction = static_friction_coefficient * normal_force;

        // Compute the magnitude of external force
        f32 external_force_magnitude = v3_magnitude(total_external_force);

        if (external_force_magnitude < max_static_friction) {
            // Cancel out external forces to prevent motion
            body->rigid_body->force = (v3){0.0f, 0.0f, 0.0f};
        } else {
            // Apply static friction force to oppose motion
            v3 external_force_normalized = v3_scale(total_external_force, 1.0f / external_force_magnitude);
            v3 friction_force = v3_scale(external_force_normalized, -max_static_friction);
            body->rigid_body->force = v3_add(body->rigid_body->force, friction_force);
        }
    }
}

void apply_surface_restitution(PhysicsBody* body, const PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL || body->collider == NULL) {
        return;
    }

    // Get the restitution coefficients
    f32 body_restitution = body->collider->restitution_coefficient; // Ensure this field exists
    f32 floor_restitution = floor_zone->restitution_coefficient;    // Ensure this field exists

    // Combine the restitution coefficients (assuming average)
    f32 restitution_coefficient = (body_restitution + floor_restitution) * 0.5f;

    // Determine the floor's Y position
    f32 floor_y = floor_zone->plane_y; // Ensure this field exists

    // Get the body's lowest point along the Y-axis
    f32 body_radius = body->collider->bounding_sphere_radius; // Ensure this field exists

    // Calculate the lowest point of the body
    f32 lowest_point_y = body->transform.position.y - body_radius;

    // Check if the body is penetrating the floor
    if (lowest_point_y < floor_y) {
        // Correct the position to prevent interpenetration
        f32 penetration_depth = floor_y - lowest_point_y;
        body->transform.position.y += penetration_depth;
        body->rigid_body->position.y = body->transform.position.y;

        // Reflect the Y component of velocity and apply restitution
        if (body->rigid_body->velocity.y < 0.0f) {
            body->rigid_body->velocity.y = -body->rigid_body->velocity.y * restitution_coefficient;
        }

        // Optionally, apply friction to the tangential velocities (X and Z)
        f32 friction_damping = 0.8f; // Adjust as needed
        body->rigid_body->velocity.x *= friction_damping;
        body->rigid_body->velocity.z *= friction_damping;
    }
}

void physics_body_adjust_for_wet_surface(PhysicsBody* body, PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL) {
        return;
    }

    // Increase slipperiness by reducing friction
    f32 wet_friction_factor = 0.5f; // Example value; adjust as needed
    floor_zone->friction_override *= wet_friction_factor;

    // Re-apply friction with the updated value
    apply_friction_for_surface(body, floor_zone);
}

void physics_body_apply_ice_friction(PhysicsBody* body, PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL) {
        return;
    }

    // Set friction to a very low value
    floor_zone->friction_override = 0.01f; // Near frictionless

    // Re-apply friction with the new value
    apply_friction_for_surface(body, floor_zone);
}

void physics_body_update_surface_interaction(PhysicsBody* body, PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL) {
        return;
    }

    // Apply friction based on the surface
    apply_friction_for_surface(body, floor_zone);

    // Apply restitution effects if any
    apply_surface_restitution(body, floor_zone);

    // Adjust for specific surface types
    if (floor_zone->is_wet) {
        physics_body_adjust_for_wet_surface(body, floor_zone);
    }

    if (floor_zone->is_icy) {
        physics_body_apply_ice_friction(body, floor_zone);
    }

    // Add more surface interactions here
}

void physics_body_adjust_for_slippery_surface(PhysicsBody* body, const PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL) {
        return;
    }

    // Use a lower friction coefficient to simulate slipperiness
    f32 slippery_friction_factor = 0.3f; // Adjust this value as needed

    // Adjust the friction coefficient
    f32 adjusted_friction_coefficient = floor_zone->friction_override * slippery_friction_factor;

    // Apply friction force
    f32 normal_force = body->rigid_body->mass * v3_magnitude(physics_state.settings.gravity);
    f32 friction_magnitude = adjusted_friction_coefficient * normal_force;
    v3 velocity_normalized = v3_normalize(body->rigid_body->velocity);
    v3 friction_force = v3_scale(velocity_normalized, -friction_magnitude);

    // Accumulate the friction force
    body->rigid_body->force = v3_add(body->rigid_body->force, friction_force);
}

void physics_body_interact_with_sticky_surface(PhysicsBody* body, const PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL) {
        return;
    }

    // Increase friction coefficient to simulate stickiness
    f32 sticky_friction_factor = 2.0f; // Adjust this value as needed

    // Adjust the friction coefficient
    f32 adjusted_friction_coefficient = floor_zone->friction_override * sticky_friction_factor;

    // Apply friction force
    f32 normal_force = body->rigid_body->mass * v3_magnitude(physics_state.settings.gravity);
    f32 friction_magnitude = adjusted_friction_coefficient * normal_force;
    v3 velocity_normalized = v3_normalize(body->rigid_body->velocity);
    v3 friction_force = v3_scale(velocity_normalized, -friction_magnitude);

    // Accumulate the friction force
    body->rigid_body->force = v3_add(body->rigid_body->force, friction_force);

    // Apply additional damping
    body->rigid_body->velocity = v3_scale(body->rigid_body->velocity, 0.5f); // Config or script to set
}

void physics_body_simulate_rough_surface(PhysicsBody* body, const PhysicsZone* floor_zone, f32 roughness_factor) {
    if (body == NULL || floor_zone == NULL || body->rigid_body == NULL) {
        return;
    }

    // Increase friction coefficient to simulate roughness
    f32 adjusted_friction_coefficient = floor_zone->friction_override * roughness_factor;

    // Apply friction force
    f32 normal_force = body->rigid_body->mass * v3_magnitude(physics_state.settings.gravity);
    f32 friction_magnitude = adjusted_friction_coefficient * normal_force;
    v3 velocity_normalized = v3_normalize(body->rigid_body->velocity);
    v3 friction_force = v3_scale(velocity_normalized, -friction_magnitude);

    // Accumulate the friction force
    body->rigid_body->force = v3_add(body->rigid_body->force, friction_force);
}

void physics_body_apply_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor) {
    if (body == NULL || force == NULL) {
        return;
    }

    // Clamp amplification_factor to max_force_amplification
    if (body->max_force_amplification > 0.0f) {
        amplification_factor = fminf(amplification_factor, body->max_force_amplification);
    }

    // Amplify the force
    v3 amplified_force = v3_scale(*force, amplification_factor);

    // Accumulate the amplified force
    body->force = v3_add(body->force, amplified_force);
}

void physics_body_apply_amplified_impulse(RigidBody* body, const v3* impulse, const v3* contact_point, f32 amplification_factor) {
    if (body == NULL || impulse == NULL || contact_point == NULL) {
        return;
    }

    // Clamp amplification_factor to max_impulse_amplification
    if (body->max_impulse_amplification > 0.0f) {
        amplification_factor = fminf(amplification_factor, body->max_impulse_amplification);
    }

    // Amplify the impulse
    v3 amplified_impulse = v3_scale(*impulse, amplification_factor);

    // Update linear velocity
    v3 delta_v = v3_scale(amplified_impulse, body->inverse_mass);
    body->velocity = v3_add(body->velocity, delta_v);

    // Calculate torque
    v3 r = v3_subtract(*contact_point, body->position);
    v3 torque = v3_cross_product(r, amplified_impulse);

    // Update angular velocity
    v3 delta_omega = m3_multiply_v3(body->inverse_inertia_tensor, torque);
    body->angular_velocity = v3_add(body->angular_velocity, delta_omega);
}

void physics_body_apply_combined_surface_effects(PhysicsBody* body, PhysicsZone* floor_zone) {
    if (body == NULL || floor_zone == NULL) {
        return;
    }

    // Apply friction based on the surface
    apply_friction_for_surface(body, floor_zone);

    // Apply restitution effects if any
    apply_surface_restitution(body, floor_zone);

    // Adjust for specific surface types
    if (floor_zone->is_wet) {
        physics_body_adjust_for_wet_surface(body, floor_zone);
    }

    if (floor_zone->is_icy) {
        physics_body_apply_ice_friction(body, floor_zone);
    }

    if (floor_zone->is_slippery) {
        physics_body_adjust_for_slippery_surface(body, floor_zone);
    }

    if (floor_zone->is_sticky) {
        physics_body_interact_with_sticky_surface(body, floor_zone);
    }

    if (floor_zone->is_rough) {
        f32 roughness_factor = 1.5f; // Config or script to drive this
        physics_body_simulate_rough_surface(body, floor_zone, roughness_factor);
    }

    // Add more surface interactions
}

void physics_body_check_surface_transition(PhysicsBody* body, const PhysicsZone* previous_zone, const PhysicsZone* current_zone) {
    if (body == NULL || previous_zone == NULL || current_zone == NULL) {
        return;
    }

    if (previous_zone != current_zone) {
        // Handle transitions between different surface types
        // Such as eset friction coefficients or adjust velocities

        // TODO: Transition from sticky to slippery surface
        if (previous_zone->is_sticky && current_zone->is_slippery) {
            // Reset velocity damping applied by sticky surface
            // Apply any necessary adjustments
        }

        // Implement other transition logic
    }
}

void physics_body_apply_directional_force(RigidBody* body, const v3* force, const v3* direction, f32 amplification_factor) {
    if (body == NULL || force == NULL || direction == NULL) {
        return;
    }

    // Normalize the direction vector
    v3 normalized_direction = v3_normalize(*direction);

    // Calculate the magnitude of the force
    f32 force_magnitude = v3_magnitude(*force);

    // Apply the force in the specified direction
    v3 directional_force = v3_scale(normalized_direction, force_magnitude * amplification_factor);

    // Accumulate the amplified force
    body->force = v3_add(body->force, directional_force);
}

void physics_body_apply_timed_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor, f32 duration) {
    if (body == NULL || force == NULL || duration <= 0.0f) {
        return;
    }

    // Clamp amplification_factor to max_force_amplification
    if (body->max_force_amplification > 0.0f) {
        amplification_factor = fminf(amplification_factor, body->max_force_amplification);
    }

    // Allocate memory for a new timed force
    const u32 MAX_TIMED_FORCES = 10;
    if (body->timed_force_count >= MAX_TIMED_FORCES) {
        // Cannot add more timed forces
        return;
    }

    if (body->timed_forces == NULL) {
        body->timed_forces = (TimedForce*)malloc(MAX_TIMED_FORCES * sizeof(TimedForce));
        if (body->timed_forces == NULL) {
            // Memory allocation failed
            return;
        }
        memset(body->timed_forces, 0, MAX_TIMED_FORCES * sizeof(TimedForce));
    }

    // Add the timed force
    TimedForce* new_force = &body->timed_forces[body->timed_force_count++];
    new_force->force = v3_scale(*force, amplification_factor);
    new_force->duration = duration;
    new_force->amplification_factor = amplification_factor;
}

void physics_body_apply_area_impulse(PhysicsBody* body, const v3* origin, f32 radius, f32 amplification_factor) {
    if (body == NULL || origin == NULL || body->rigid_body == NULL) {
        return;
    }

    // Calculate the distance from the body to the origin
    v3 distance_vector = v3_subtract(body->transform.position, *origin);
    f32 distance = v3_magnitude(distance_vector);

    if (distance <= radius) {
        // Calculate impulse direction (away from the origin)
        v3 impulse_direction = v3_normalize(distance_vector);

        // Calculate impulse magnitude
        f32 impulse_magnitude = (1.0f - (distance / radius)) * amplification_factor;

        // Create the impulse vector
        v3 impulse = v3_scale(impulse_direction, impulse_magnitude);

        // Apply the impulse to the rigid body
        body->rigid_body->velocity = v3_add(body->rigid_body->velocity, v3_scale(impulse, body->rigid_body->inverse_mass));
    }
}

void physics_body_set_amplification_limits(RigidBody* body, f32 max_force, f32 max_impulse) {
    if (body == NULL) {
        return;
    }

    body->max_force_amplification = max_force;
    body->max_impulse_amplification = max_impulse;
}

void physics_body_apply_decay_amplified_force(RigidBody* body, const v3* force, f32 amplification_factor, f32 decay_rate) {
    if (body == NULL || force == NULL || decay_rate <= 0.0f) {
        return;
    }

    // Clamp amplification_factor to max_force_amplification
    if (body->max_force_amplification > 0.0f) {
        amplification_factor = fminf(amplification_factor, body->max_force_amplification);
    }

    // Allocate memory for a new decay force
    const u32 MAX_DECAY_FORCES = 10;
    if (body->decay_force_count >= MAX_DECAY_FORCES) {
        // Cannot add more decay forces
        return;
    }

    if (body->decay_forces == NULL) {
        body->decay_forces = (DecayForce*)malloc(MAX_DECAY_FORCES * sizeof(DecayForce));
        if (body->decay_forces == NULL) {
            // Memory allocation failed
            return;
        }
        memset(body->decay_forces, 0, MAX_DECAY_FORCES * sizeof(DecayForce));
    }

    // Add the decay force
    DecayForce* new_force = &body->decay_forces[body->decay_force_count++];
    new_force->force = *force;
    new_force->amplification_factor = amplification_factor;
    new_force->decay_rate = decay_rate;
}

void physics_body_update_timed_forces(RigidBody* body, f32 delta_time) {
    if (body == NULL || body->timed_forces == NULL || body->timed_force_count == 0) {
        return;
    }

    for (u32 i = 0; i < body->timed_force_count; ) {
        TimedForce* timed_force = &body->timed_forces[i];

        // Apply the force
        body->force = v3_add(body->force, timed_force->force);

        // Decrease the duration
        timed_force->duration -= delta_time;

        if (timed_force->duration <= 0.0f) {
            // Remove the timed force by swapping with the last one
            body->timed_forces[i] = body->timed_forces[--body->timed_force_count];
        } else {
            i++;
        }
    }

    // If all timed forces are used up, free the memory
    if (body->timed_force_count == 0) {
        free(body->timed_forces);
        body->timed_forces = NULL;
    }
}

void physics_body_update_decay_forces(RigidBody* body, f32 delta_time) {
    if (body == NULL || body->decay_forces == NULL || body->decay_force_count == 0) {
        return;
    }

    for (u32 i = 0; i < body->decay_force_count; ) {
        DecayForce* decay_force = &body->decay_forces[i];

        // Apply the force
        v3 amplified_force = v3_scale(decay_force->force, decay_force->amplification_factor);
        body->force = v3_add(body->force, amplified_force);

        // Decay the amplification factor
        decay_force->amplification_factor -= decay_force->decay_rate * delta_time;

        if (decay_force->amplification_factor <= 0.0f) {
            // Remove the decay force by swapping with the last one
            body->decay_forces[i] = body->decay_forces[--body->decay_force_count];
        } else {
            i++;
        }
    }

    // If all decay forces are used up, free the memory
    if (body->decay_force_count == 0) {
        free(body->decay_forces);
        body->decay_forces = NULL;
    }
}
