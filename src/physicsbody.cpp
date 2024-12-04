/*
 * Phyzzy: Playful Physics for Imaginative Games
 *
 * Licensed under the GNU General Public License, Version 3.
 * For license details, visit: https://www.gnu.org/licenses/gpl-3.0.html
 *
 * Questions or contributions? Reach out to Simon Devenish:
 * simon.devenish@outlook.com
 */

#include <utility>      // For std::move

#include "physicsbody.h"
#include "physicsbody_physics.h"    // For C functions
#include "vector_math.h" // For vector and quaternion operations

namespace phyzzy
{
    // Constructor
PhysicsBodyThing::PhysicsBodyThing(::PhysicsBody* body)
    : c_body(body, physics_body_destroy) // Set the custom deleter
{
}

// Move constructor
PhysicsBodyThing::PhysicsBodyThing(PhysicsBodyThing&& other) noexcept
    : c_body(std::move(other.c_body))
{
    // The moved-from object now has a nullptr in c_body
}

// Move assignment operator
PhysicsBodyThing& PhysicsBodyThing::operator=(PhysicsBodyThing&& other) noexcept
{
    if (this != &other)
    {
        c_body = std::move(other.c_body);
    }
    return *this;
}

// Private helper function to get the RigidBody*
::RigidBody* PhysicsBodyThing::GetRigidBody() const
{
    return c_body ? c_body->rigid_body : nullptr;
}

void PhysicsBodyThing::Destroy()
{
    c_body.reset(); // Calls the custom deleter if c_body is not nullptr
}

void PhysicsBodyThing::SetMass(const f32 mass)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_set_mass(rigidBody, mass);
    }
}

f32 PhysicsBodyThing::GetMass() const
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        f32 mass = 0.0f;
        physics_body_get_mass(rigidBody, &mass);
        return mass;
    }
    return 0.0f;
}

void PhysicsBodyThing::SetVelocity(const v3& velocity)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_set_velocity(rigidBody, &velocity);
    }
}

const v3 PhysicsBodyThing::GetVelocity() const
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        v3 velocity = {0.0f, 0.0f, 0.0f};
        physics_body_get_velocity(rigidBody, &velocity);
        return velocity;
    }
    return {0.0f, 0.0f, 0.0f};
}

void PhysicsBodyThing::SetAngularVelocity(const v3& angularVelocity)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_set_angular_velocity(rigidBody, &angularVelocity);
    }
}

const v3 PhysicsBodyThing::GetAngularVelocity() const
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        v3 angularVelocity = {0.0f, 0.0f, 0.0f};
        physics_body_get_angular_velocity(rigidBody, &angularVelocity);
        return angularVelocity;
    }
    return {0.0f, 0.0f, 0.0f};
}

void PhysicsBodyThing::ResetForces()
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_reset_forces(rigidBody);
    }
}

void PhysicsBodyThing::Translate(const v3& translation)
{
    if (c_body)
    {
        physics_body_translate(c_body.get(), &translation);
    }
}

void PhysicsBodyThing::Rotate(const quat& rotation)
{
    if (c_body)
    {
        physics_body_rotate(c_body.get(), &rotation);
    }
}

void PhysicsBodyThing::SetPosition(const v3& position)
{
    if (c_body)
    {
        physics_body_set_position(c_body.get(), &position);
    }
}

const v3 PhysicsBodyThing::GetPosition() const
{
    if (c_body)
    {
        v3 position = {0.0f, 0.0f, 0.0f};
        physics_body_get_position(c_body.get(), &position);
        return position;
    }
    return {0.0f, 0.0f, 0.0f};
}

void PhysicsBodyThing::SetStatic(const bool isStatic)
{
    if (c_body)
    {
        physics_body_set_static(c_body.get(), isStatic);
    }
}

bool PhysicsBodyThing::IsStatic() const
{
    if (c_body)
    {
        return physics_body_is_static(c_body.get());
    }
    return false;
}

void PhysicsBodyThing::ApplyDamping(const f32 linearDamping, const f32 angularDamping)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_damping(rigidBody, linearDamping, angularDamping);
    }
}

void PhysicsBodyThing::ApplyDrag(const f32 dragCoefficient)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_drag(rigidBody, dragCoefficient);
    }
}

void PhysicsBodyThing::Freeze()
{
    if (c_body)
    {
        physics_body_freeze(c_body.get());
    }
}

void PhysicsBodyThing::Unfreeze()
{
    if (c_body)
    {
        physics_body_unfreeze(c_body.get());
    }
}

void PhysicsBodyThing::ApplyDirectionalWobble(const v3& direction, f32 amplitude, f32 frequency)
{
    if (c_body)
    {
        physics_body_apply_directional_wobble(c_body.get(), &direction, amplitude, frequency);
    }
}

void PhysicsBodyThing::SyncWobbleWithAnimation(f32 animationPhase)
{
    if (c_body)
    {
        physics_body_sync_wobble_with_animation(c_body.get(), animationPhase);
    }
}

void PhysicsBodyThing::ApplyImpactWobble(const v3& impactForce, f32 decayRate)
{
    if (c_body)
    {
        physics_body_apply_impact_wobble(c_body.get(), &impactForce, decayRate);
    }
}

void PhysicsBodyThing::GetWobbleState(f32& amplitudeOut, f32& frequencyOut) const
{
    if (c_body)
    {
        physics_body_get_wobble_state(c_body.get(), &amplitudeOut, &frequencyOut);
    }
}

void PhysicsBodyThing::ToggleWobble(bool enable)
{
    if (c_body)
    {
        physics_body_toggle_wobble(c_body.get(), enable);
    }
}

void PhysicsBodyThing::UpdateSurfaceInteraction(PhysicsZone& floorZone)
{
    if (c_body)
    {
        physics_body_update_surface_interaction(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::ApplyFrictionForSurface(PhysicsZone& floorZone)
{
    if (c_body)
    {
        apply_friction_for_surface(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::ApplySurfaceRestitution(PhysicsZone& floorZone)
{
    if (c_body)
    {
        apply_surface_restitution(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::AdjustForWetSurface(PhysicsZone& floorZone)
{
    if (c_body)
    {
        physics_body_adjust_for_wet_surface(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::ApplyIceFriction(PhysicsZone& floorZone)
{
    if (c_body)
    {
        physics_body_apply_ice_friction(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::AdjustForSlipperySurface(const PhysicsZone& floorZone)
{
    if (c_body)
    {
        physics_body_adjust_for_slippery_surface(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::InteractWithStickySurface(const PhysicsZone& floorZone)
{
    if (c_body)
    {
        physics_body_interact_with_sticky_surface(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::SimulateRoughSurface(const PhysicsZone& floorZone, f32 roughnessFactor)
{
    if (c_body)
    {
        physics_body_simulate_rough_surface(c_body.get(), &floorZone, roughnessFactor);
    }
}

void PhysicsBodyThing::ApplyCombinedSurfaceEffects(PhysicsZone& floorZone)
{
    if (c_body)
    {
        physics_body_apply_combined_surface_effects(c_body.get(), &floorZone);
    }
}

void PhysicsBodyThing::CheckSurfaceTransition(const PhysicsZone& previousZone, const PhysicsZone& currentZone)
{
    if (c_body)
    {
        physics_body_check_surface_transition(c_body.get(), &previousZone, &currentZone);
    }
}

void PhysicsBodyThing::ApplyAmplifiedForce(const v3& force, f32 amplificationFactor)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_amplified_force(rigidBody, &force, amplificationFactor);
    }
}

void PhysicsBodyThing::ApplyAmplifiedImpulse(const v3& impulse, const v3& contactPoint, f32 amplificationFactor)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_amplified_impulse(rigidBody, &impulse, &contactPoint, amplificationFactor);
    }
}

void PhysicsBodyThing::ApplyDirectionalForce(const v3& force, const v3& direction, f32 amplificationFactor)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_directional_force(rigidBody, &force, &direction, amplificationFactor);
    }
}

void PhysicsBodyThing::ApplyTimedAmplifiedForce(const v3& force, f32 amplificationFactor, f32 duration)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_timed_amplified_force(rigidBody, &force, amplificationFactor, duration);
    }
}

void PhysicsBodyThing::ApplyAreaImpulse(const v3& origin, f32 radius, f32 amplificationFactor)
{
    if (c_body)
    {
        physics_body_apply_area_impulse(c_body.get(), &origin, radius, amplificationFactor);
    }
}

void PhysicsBodyThing::SetAmplificationLimits(f32 maxForce, f32 maxImpulse)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_set_amplification_limits(rigidBody, maxForce, maxImpulse);
    }
}

void PhysicsBodyThing::ApplyDecayAmplifiedForce(const v3& force, f32 amplificationFactor, f32 decayRate)
{
    ::RigidBody* rigidBody = GetRigidBody();
    if (rigidBody)
    {
        physics_body_apply_decay_amplified_force(rigidBody, &force, amplificationFactor, decayRate);
    }
}
} // namespace phyzzy