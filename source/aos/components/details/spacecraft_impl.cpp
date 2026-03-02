#include "spacecraft_impl.hpp"

#include "aos/components/hysteresis_rods.hpp"
#include "aos/components/inertia_tensor.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft.hpp"
#include "aos/components/spacecraft_faces.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

namespace aos {

spacecraft_impl::spacecraft_impl(const spacecraft_properties& properties)
    : _mass_kg(properties.mass_kg),
      _inertia(properties.mass_kg, properties.shape),
      _faces(properties.shape),
      _magnet(properties.magnet),
      _hystresis(properties.rods, properties.hysteresis) {}

spacecraft_impl::~spacecraft_impl() = default;

auto spacecraft_impl::mass_kg() const -> double {
    return _mass_kg;
}

auto spacecraft_impl::inertia() const -> const inertia_tensor& {
    return _inertia;
}

auto spacecraft_impl::faces() const -> const spacecraft_faces& {
    return _faces;
}

auto spacecraft_impl::magnet() const -> const permanent_magnet& {
    return _magnet;
}

auto spacecraft_impl::hystresis() const -> const hysteresis_rods& {
    return _hystresis;
}

void spacecraft_impl::derivative(const environment_effects& env, const system_state& current_state, system_state& state_derivative) const {
    const vec3& r_eci            = current_state.position_m;
    const vec3& v_eci            = current_state.velocity_m_s;
    const quat  q_att            = current_state.attitude.normalized();  // normalize to prevent drift
    const vec3& omega_body       = current_state.angular_velocity_m_s;
    const quat  q_inv            = q_att.conjugate();
    const vec3  r_body           = q_inv * r_eci;
    const vec3  b_body           = q_inv * env.magnetic_field_eci_T;
    const vec3  b_dot_orbital    = q_inv * env.magnetic_field_dot_eci_T_s;
    const vec3  b_dot_rotational = -omega_body.cross(b_body);
    const vec3  b_dot_body       = b_dot_orbital + b_dot_rotational;
    const vec3  rods_torque      = _hystresis.compute_rod_torques(current_state.rod_magnetizations, b_body);
    const auto  face_effects     = _faces.compute_face_effects(env, q_att, q_inv, omega_body);
    const vec3  net_torque       = compute_torques(omega_body, b_body, r_body, env.earth_mu) + rods_torque + face_effects.torque_body;

    state_derivative.position_m   = v_eci;
    state_derivative.velocity_m_s = env.gravity_eci_m_s2;
    state_derivative.velocity_m_s += face_effects.force_eci / _mass_kg;
    state_derivative.angular_velocity_m_s = _inertia.inverse() * net_torque;
    state_derivative.attitude.coeffs()    = system_state::compute_attitude_derivative(q_att, omega_body);
    _hystresis.compute_rod_derivatives(current_state.rod_magnetizations, b_body, b_dot_body, state_derivative.rod_magnetizations);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto spacecraft_impl::compute_torques(const vec3& omega, const vec3& b_body, const vec3& r_body, double earth_mu) const -> vec3 {
    vec3 torque = vec3::Zero();
    torque += _magnet.compute_torque(b_body);
    torque += _inertia.compute_gyroscopic_torque(omega);
    torque += _inertia.compute_gravity_gradient_torque(r_body, earth_mu);
    return torque;
}

}  // namespace aos
