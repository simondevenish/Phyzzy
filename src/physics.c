
/*
 * Phyzzy: Playful Physics for Imaginative Games
 *
 * Licensed under the GNU General Public License, Version 3.
 * For license details, visit: https://www.gnu.org/licenses/gpl-3.0.html
 *
 * Questions or contributions? Reach out to Simon Devenish:
 * simon.devenish@outlook.com
 */

#include "physics.h"
#include "vector_math.h"
#include <stdlib.h>
#include <string.h> // For memset

// Define global Phyzzy state variables
PhysicsSettings g_physics_settings;

void physics_initialize(PhysicsSettings* settings)
{
    if (physics_state.initialized) {
        // Phyzzy system is already initialized
        return;
    }

    // Copy the provided settings or use default settings if NULL
    if (settings != NULL) {
        physics_state.settings = *settings;
    } else {
        // Set default physics settings
        physics_state.settings.time_step = 0.016f;  // Approximate for 60 FPS
        physics_state.settings.max_iterations = 10;
        physics_state.settings.enable_gravity = true;
        physics_state.settings.gravity = (v3){0.0f, -9.81f, 0.0f};
        physics_state.settings.max_bodies = 1024;
        physics_state.settings.max_colliders = 2048;
        physics_state.settings.max_joints = 256;
        physics_state.settings.collision_tolerance = 0.01f;
        physics_state.settings.enable_sleeping = true;
        physics_state.settings.sleep_threshold = 0.05f;
        physics_state.settings.linear_damping = 0.1f;
        physics_state.settings.angular_damping = 0.1f;
        physics_state.settings.max_zones = 10;
        physics_state.settings.enable_debug_logs = false;
    }

    // Initialize internal data structures
    // For example, allocate memory for bodies, colliders, and joints
    // Since internal_data is a placeholder, we leave it as NULL for now
    physics_state.internal_data = NULL;

    // Mark the Phyzzy as initialized
    physics_state.initialized = true;
}

void physics_cleanup()
{
    if (!physics_state.initialized) {
        // Phyzzy is not initialized
        return;
    }

    // Clean up internal data structures
    // Since internal_data is a placeholder, we leave it as NULL for now
    if (physics_state.internal_data != NULL) {
        // Free or deallocate internal_data
        physics_state.internal_data = NULL;
    }

    // Reset physics settings to defaults if necessary
    // For now, we can zero out the settings
    memset(&physics_state.settings, 0, sizeof(PhysicsSettings));

    // Mark the physics system as uninitialized
    physics_state.initialized = false;
}
