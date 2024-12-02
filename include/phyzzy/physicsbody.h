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

#include <memory>

#include "physics_types.h"
#include "vector_math.h"

namespace phyzzy
{
class PhysicsBodyThing
{
private:
    std::unique_ptr<::PhysicsBody, void(*)(::PhysicsBody*)> c_body;

       ::RigidBody* GetRigidBody() const
    {
        return c_body ? c_body->rigid_body : nullptr;
    }
public:
    // Constructor
    explicit PhysicsBodyThing(::PhysicsBody* body);

    // Destructor
    ~PhysicsBodyThing() = default;

    // Deleted copy constructor and assignment operator to enforce unique ownership
    PhysicsBodyThing(const PhysicsBodyThing&) = delete;
    PhysicsBodyThing& operator=(const PhysicsBodyThing&) = delete;

    // Move constructor and assignment operator for transferring ownership
    PhysicsBodyThing(PhysicsBodyThing&& other) noexcept;
    PhysicsBodyThing& operator=(PhysicsBodyThing&& other) noexcept;

    // Physics body management
    void Destroy();
    void SetMass(const f32 mass);
    const f32 GetMass() const;
    void SetVelocity(const v3& velocity);
    const v3 GetVelocity() const;
    void SetAngularVelocity(const v3& angularVelocity);
    const v3 GetAngularVelocity() const;
    void ResetForces();
    void Translate(const v3& translation);
    void Rotate(const quat& rotation);
    void SetPosition(const v3& position);
    const v3 GetPosition() const;
    void SetStatic(const bool isStatic);
    bool IsStatic() const;
    void ApplyDamping(const f32 linearDamping, const f32 angularDamping);
    void ApplyDrag(const f32 dragCoefficient);
    void Freeze();
    void Unfreeze();

     // Wobble and oscillation effects
    void ApplyDirectionalWobble(const v3& direction, f32 amplitude, f32 frequency);
    void SyncWobbleWithAnimation(f32 animationPhase);
    void ApplyImpactWobble(const v3& impactForce, f32 decayRate);
    void GetWobbleState(f32& amplitudeOut, f32& frequencyOut) const;
    void ToggleWobble(bool enable);

    // Dynamic floor effects
    void UpdateSurfaceInteraction(const PhysicsZone& floorZone);
    void ApplyFrictionForSurface(const PhysicsZone& floorZone);
    void ApplySurfaceRestitution(const PhysicsZone& floorZone);
    void AdjustForWetSurface(const PhysicsZone& floorZone);
    void ApplyIceFriction(const PhysicsZone& floorZone);

     // Advanced floor interaction
    void AdjustForSlipperySurface(const PhysicsZone& floorZone);
    void InteractWithStickySurface(const PhysicsZone& floorZone);
    void SimulateRoughSurface(const PhysicsZone& floorZone, f32 roughnessFactor);

    // Combined floor effects
    void ApplyCombinedSurfaceEffects(const PhysicsZone& floorZone);
    void CheckSurfaceTransition(const PhysicsZone& previousZone, const PhysicsZone& currentZone);

    // Force and impulse amplification
    void ApplyAmplifiedForce(const v3& force, f32 amplificationFactor);
    void ApplyAmplifiedImpulse(const v3& impulse, const v3& contactPoint, f32 amplificationFactor);

    // Advanced Force and Impulse Amplification Features
    void ApplyDirectionalForce(const v3& force, const v3& direction, f32 amplificationFactor);
    void ApplyTimedAmplifiedForce(const v3& force, f32 amplificationFactor, f32 duration);
    void ApplyAreaImpulse(const v3& origin, f32 radius, f32 amplificationFactor);
    void SetAmplificationLimits(f32 maxForce, f32 maxImpulse);
    void ApplyDecayAmplifiedForce(const v3& force, f32 amplificationFactor, f32 decayRate);
};
} // namespace phyzzy