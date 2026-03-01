#include "dynamics.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace aos {

void spacecraft_dynamics::operator()(const system_state& current_state, system_state& state_derivative, double t_sec) const {
    const double t_global   = _global_time_offset + t_sec;
    const vec3&  r_eci      = current_state.position_m;
    const vec3&  v_eci      = current_state.velocity_m_s;
    const quat   q_att      = current_state.attitude.normalized();  // normalize to prevent drift
    const vec3&  omega_body = current_state.angular_velocity_m_s;
    const auto   env_data   = _environment->calculate(t_global, r_eci, v_eci);

    state_derivative.position_m   = v_eci;
    state_derivative.velocity_m_s = env_data.gravity_eci_m_s2;

    const mat3x3 R_eci_to_body    = q_att.toRotationMatrix().transpose();  // NOLINT(readability-identifier-naming)
    const vec3   b_body           = R_eci_to_body * env_data.magnetic_field_eci_T;
    const vec3   b_dot_orbital    = R_eci_to_body * env_data.magnetic_field_dot_eci_T_s;
    const vec3   b_dot_rotational = -omega_body.cross(b_body);
    const vec3   b_dot_body       = b_dot_orbital + b_dot_rotational;
    const vec3   rods_torque      = compute_rod_effects(current_state.rod_magnetizations, b_body, b_dot_body, state_derivative.rod_magnetizations);
    const auto   face_effects     = compute_face_effects(env_data.atmospheric_density_kg_m3, q_att, v_eci, omega_body);
    const vec3   faces_torque     = R_eci_to_body * face_effects.torque;
    const vec3   net_torque       = compute_other_torques(omega_body, b_body, r_eci, q_att) + rods_torque + faces_torque;

    // TODO: compute solar pressure

    state_derivative.velocity_m_s += face_effects.force / (_spacecraft->mass() * gram_to_kilogram);
    state_derivative.angular_velocity_m_s = _spacecraft->inertia_tensor_inverse() * net_torque;
    state_derivative.attitude.coeffs()    = compute_attitude_derivative(q_att, omega_body);
}

auto spacecraft_dynamics::compute_rod_effects(const vecX& rod_magnetizations, const vec3& b_body, const vec3& b_dot_body, vecX& dm_dt_out) const -> vec3 {
    const auto& rods     = _spacecraft->rods();
    const auto  num_rods = std::min(static_cast<std::ptrdiff_t>(rods.size()), rod_magnetizations.size());

    // ensure output vector is correctly sized
    // if (dm_dt_out.size() != num_rods) {
    //     dm_dt_out.resize(num_rods);
    // }

    vec3 total_torque = vec3::Zero();
    for (std::ptrdiff_t i = 0; i < num_rods; ++i) {
        const double m_irr = rod_magnetizations(static_cast<int>(i));

        // dM_irr/dt
        dm_dt_out(i) = rods[i].magnetization_derivative(m_irr, b_body, b_dot_body);

        // tau = m_dipole x B
        total_torque += rods[i].magnetic_moment(m_irr, b_body).cross(b_body);
    }
    return total_torque;
}

auto spacecraft_dynamics::compute_face_effects(double density, const quat& q_att, const vec3& v_eci, const vec3& omega_body) const -> face_effects {
    vec3 total_torque = vec3::Zero();
    vec3 total_force  = vec3::Zero();

    for (const auto& face : _spacecraft->faces()) {
        const auto force_drag = face.compute_force_drag(density, q_att, v_eci, omega_body);
        const auto face_eci   = (q_att * face.center_of_pressure_m);
        total_force += force_drag;
        total_torque += face_eci.cross(force_drag);
    }

    return {
        .torque = total_torque,
        .force  = total_force,
    };
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto spacecraft_dynamics::compute_other_torques(const vec3& omega, const vec3& b_body, const vec3& r_eci, const quat& q_att) const -> vec3 {
    vec3 torque = vec3::Zero();

    // permanent magnet
    torque += _spacecraft->magnet().magnetic_moment().cross(b_body);

    // gyroscopic (rigid body coupling): -omega x (I * omega)
    torque -= omega.cross(_spacecraft->inertia_tensor() * omega);

    // gravity gradient
    torque += compute_gravity_gradient_torque(r_eci, q_att);
    return torque;
}

auto spacecraft_dynamics::compute_gravity_gradient_torque(const vec3& r_eci, const quat& q_att) const -> vec3 {
    // position to body frame: r_body = R * r_eci
    const vec3 r_body = q_att.toRotationMatrix().transpose() * r_eci;

    const double r_sq  = r_body.squaredNorm();
    const double r_val = std::sqrt(r_sq);

    // tau_gg = (3*mu / r^5) * (r_body x (I * r_body))
    const double coef = (3.0 * _environment->earth_mu()) / (r_sq * r_sq * r_val);
    return coef * r_body.cross(_spacecraft->inertia_tensor() * r_body);
}

auto spacecraft_dynamics::compute_attitude_derivative(const quat& q_att, const vec3& omega) -> vecX {
    // dq/dt = 0.5 * q * omega_quat
    const quat omega_q(0, omega.x(), omega.y(), omega.z());
    return 0.5 * (q_att * omega_q).coeffs();
}

}  // namespace aos
