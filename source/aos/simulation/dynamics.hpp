#pragma once

#include "../components/spacecraft.hpp"
#include "../core/state.hpp"
#include "../environment/environment.hpp"

namespace aos::simulation {

using components::Spacecraft;
using core::SystemState;
using environment::Environment;

class SpacecraftDynamics {
public:

    SpacecraftDynamics(const Spacecraft& spacecraft, const Environment& environment) : spacecraft_(spacecraft), environment_(environment) {}

    // The operator now only deals with physics, not environment calculation
    void operator()(const SystemState& x, SystemState& dxdt, double t) const {
        // 1. Get Environment State
        const auto b_inertial = environment_.inertial_magnetic_field_at(t);

        // 2. Kinematics (unchanged)
        // ...

        // 3. Dynamics
        const auto r_inertial_to_body = x.attitude.toRotationMatrix().transpose();
        const auto b_body             = r_inertial_to_body * b_inertial;

        core::Vector3d total_torque_body = core::Vector3d::Zero();

        // Accumulate torques
        for (const auto& magnet : spacecraft_.magnets()) {
            total_torque_body += magnet.magnetic_moment().cross(b_body);
        }
        for (size_t i = 0; i < spacecraft_.rods().size(); ++i) {
            total_torque_body += spacecraft_.rods()[i].magnetic_moment(x.rod_magnetizations(i)).cross(b_body);
        }
        total_torque_body -= x.angular_velocity.cross(spacecraft_.inertia_tensor() * x.angular_velocity);

        dxdt.angular_velocity = spacecraft_.inertia_tensor_inverse() * total_torque_body;

        // 4. Hysteresis
        for (size_t i = 0; i < spacecraft_.rods().size(); ++i) {
            dxdt.rod_magnetizations(i) = spacecraft_.rods()[i].magnetization_derivative(x.rod_magnetizations(i), b_body, x.angular_velocity);
        }
    }

private:

    const Spacecraft&  spacecraft_;
    const Environment& environment_;
};

}  // namespace aos::simulation
