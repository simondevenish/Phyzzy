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

#ifdef __cplusplus
extern "C" {
#endif

// Internal physics state
static struct {
    PhysicsSettings settings; // Stores configuration settings
    bool initialized;         // Tracks if the system is initialized
    void* internal_data;      // Placeholder for internal state (e.g., spatial grids, caches)
} physics_state;

// Declare global physics state variables
extern PhysicsSettings g_physics_settings;

// Physics system lifecycle
void physics_initialize(PhysicsSettings* settings);
void physics_cleanup(void);

#ifdef __cplusplus
}
#endif
