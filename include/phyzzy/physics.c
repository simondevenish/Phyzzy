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
#include <stdlib.h>
#include <string.h> // For memset

// Define global physics state variables
PhysicsSettings g_physics_settings;
bool g_is_initialized = false;