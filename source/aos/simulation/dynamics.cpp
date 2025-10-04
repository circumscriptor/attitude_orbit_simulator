#include "dynamics.hpp"

#include "aos/core/state.hpp"
#include "aos/core/types.hpp"

#include <cstddef>

void aos::simulation::spacecraft_dynamics::operator()(const system_state& current_state, system_state& state_derivative, double t) const {
    // --- 1. Get Environment State ---
    // Calculate the Earth's magnetic field in the inertial frame at the current time.
    const vec3 b_inertial = _environment->inertial_magnetic_field_at(t);

    // --- 2. Attitude Kinematics ---
    // Calculate the rate of change of the attitude quaternion (dq/dt).
    // This is based on the equation: q_dot = 0.5 * q * omega_q
    const auto& q     = current_state.attitude;
    const vec3& omega = current_state.angular_velocity;

    // Represent angular velocity as a pure quaternion [0, ωx, ωy, ωz]
    quat omega_q(0, omega.x(), omega.y(), omega.z());

    // The result of the quaternion multiplication gives the derivative
    state_derivative.attitude.coeffs() = 0.5 * (q * omega_q).coeffs();

    // --- 3. Rotational Dynamics ---
    // Calculate the angular acceleration (dω/dt) by summing all torques.
    const auto r_inertial_to_body = q.toRotationMatrix().transpose();
    const auto b_body             = r_inertial_to_body * b_inertial;

    vec3 total_torque_body = vec3::Zero();

    // a) Torque from the permanent magnet
    total_torque_body += _spacecraft->magnet().magnetic_moment().cross(b_body);

    // b) Torque from all hysteresis rods
    for (std::size_t i = 0; i < _spacecraft->rods().size(); ++i) {
        total_torque_body += _spacecraft->rods()[i].magnetic_moment(current_state.rod_magnetizations(static_cast<int>(i))).cross(b_body);
    }

    // c) Gyroscopic torque (always present for a rotating body)
    total_torque_body -= omega.cross(_spacecraft->inertia_tensor() * omega);

    // Store angular velocity derivative
    state_derivative.angular_velocity = _spacecraft->inertia_tensor_inverse() * total_torque_body;

    // --- 4. Hysteresis Rod Dynamics ---
    // Calculate the rate of change of each rod's internal magnetization state (dM/dt).
    state_derivative.rod_magnetizations.resize(static_cast<int>(_spacecraft->rods().size()));
    for (std::size_t i = 0; i < _spacecraft->rods().size(); ++i) {
        const auto j{static_cast<int>(i)};
        state_derivative.rod_magnetizations(j) = _spacecraft->rods()[i].magnetization_derivative(current_state.rod_magnetizations(j), b_body, omega);
    }
}
