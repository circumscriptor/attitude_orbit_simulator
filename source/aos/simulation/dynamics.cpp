#include "dynamics.hpp"

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
    const auto   face_effects     = compute_face_effects(env_data, q_att, v_eci, r_eci, omega_body);
    const vec3   net_torque       = compute_other_torques(omega_body, b_body, r_eci, R_eci_to_body) + rods_torque + face_effects.torque_body;

    // if (std::isnan(net_torque.x()) || std::isnan(state_derivative.angular_velocity_m_s.x())) {
    //     std::println(stderr, "[DYNAMICS ERROR] NaN detected at t={:.4f}", t_global);
    //     std::println(stderr, "  Rods Torque:  [{:f}, {:f}, {:f}]", rods_torque.x(), rods_torque.y(), rods_torque.z());
    //     std::println(stderr, "  Faces Torque: [{:f}, {:f}, {:f}]", faces_torque.x(), faces_torque.y(), faces_torque.z());
    //     std::println(stderr, "  Net Torque:   [{:f}, {:f}, {:f}]", net_torque.x(), net_torque.y(), net_torque.z());
    //     std::println(stderr, "  Omega Body:   [{:f}, {:f}, {:f}]", omega_body.x(), omega_body.y(), omega_body.z());
    //     for (int i = 0; i < current_state.rod_magnetizations.size(); ++i) {
    //         std::println(stderr, "  Rod[{}] M_dot: {:f}", i, state_derivative.rod_magnetizations[i]);
    //     }
    // }

    state_derivative.velocity_m_s += face_effects.force_eci / (_spacecraft->mass_g() * gram_to_kilogram);
    state_derivative.angular_velocity_m_s = _spacecraft->inertia_tensor_kg_m2_inverse() * net_torque;
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

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
auto spacecraft_dynamics::compute_face_effects(const environment_data& data,
                                               const quat&             q_att,
                                               const vec3&             v_eci,
                                               const vec3&             r_eci,
                                               const vec3&             omega_body) const -> face_effects {
    vec3 total_torque_body = vec3::Zero();
    vec3 total_force_body  = vec3::Zero();

    const quat q_inv  = q_att.conjugate();
    const vec3 s_body = (q_inv * data.r_sun_eci).normalized();
    const vec3 v_body = q_inv * environment_model::earth_relative_v(v_eci, r_eci);

    for (const auto& face : _spacecraft->faces()) {
        const vec3 v_rel_body   = face.compute_v_rel_body(v_body, omega_body);
        const vec3 f_drag_body  = face.compute_force_drag_body(data.atmospheric_density_kg_m3, v_rel_body);
        const vec3 f_srp_body   = face.compute_force_srp_body(data.solar_pressure_Pa, s_body, data.shadow_factor);
        const vec3 f_total_body = f_drag_body + f_srp_body;

        total_force_body += f_total_body;
        total_torque_body += face.center_of_pressure_m.cross(f_total_body);
    }

    return {
        .torque_body = total_torque_body,
        .force_eci   = q_att * total_force_body,
    };
}
// NOLINTEND(bugprone-easily-swappable-parameters)

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto spacecraft_dynamics::compute_other_torques(const vec3& omega, const vec3& b_body, const vec3& r_eci, const mat3x3& eci_to_body) const -> vec3 {
    vec3 torque = vec3::Zero();

    // permanent magnet
    torque += _spacecraft->magnet().magnetic_moment().cross(b_body);

    // gyroscopic (rigid body coupling): -omega x (I * omega)
    torque -= omega.cross(_spacecraft->inertia_tensor_kg_m2() * omega);

    // gravity gradient
    torque += compute_gravity_gradient_torque(r_eci, eci_to_body);
    return torque;
}

auto spacecraft_dynamics::compute_gravity_gradient_torque(const vec3& r_eci, const mat3x3& eci_to_body) const -> vec3 {
    // position to body frame: r_body = R * r_eci
    const vec3 r_body = eci_to_body * r_eci;

    const double r_sq  = r_body.squaredNorm();
    const double r_val = std::sqrt(r_sq);

    // tau_gg = (3*mu / r^5) * (r_body x (I * r_body))
    const double coef = (3.0 * _environment->earth_mu()) / (r_sq * r_sq * r_val);
    return coef * r_body.cross(_spacecraft->inertia_tensor_kg_m2() * r_body);
}

auto spacecraft_dynamics::compute_attitude_derivative(const quat& q_att, const vec3& omega) -> vecX {
    // dq/dt = 0.5 * q * omega_quat
    const quat omega_q(0, omega.x(), omega.y(), omega.z());
    return 0.5 * (q_att * omega_q).coeffs();
}

}  // namespace aos
