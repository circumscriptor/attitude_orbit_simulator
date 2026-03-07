#include "state.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

namespace aos {

auto system_state::operator+=(const system_state& other) -> system_state& {
    position_m += other.position_m;
    velocity_m_s += other.velocity_m_s;
    attitude.coeffs() += other.attitude.coeffs();
    angular_velocity_m_s += other.angular_velocity_m_s;
    rod_magnetizations += other.rod_magnetizations;
    return *this;
}

auto system_state::operator+=(real scalar) -> system_state& {
    position_m.array() += scalar;
    velocity_m_s.array() += scalar;
    attitude.coeffs().array() += scalar;
    angular_velocity_m_s.array() += scalar;
    rod_magnetizations.array() += scalar;
    return *this;
}

auto system_state::operator-=(const system_state& other) -> system_state& {
    position_m -= other.position_m;
    velocity_m_s -= other.velocity_m_s;
    attitude.coeffs() -= other.attitude.coeffs();
    angular_velocity_m_s -= other.angular_velocity_m_s;
    rod_magnetizations -= other.rod_magnetizations;
    return *this;
}

auto system_state::operator*=(real scalar) -> system_state& {
    position_m *= scalar;
    velocity_m_s *= scalar;
    attitude.coeffs() *= scalar;
    angular_velocity_m_s *= scalar;
    rod_magnetizations *= scalar;
    return *this;
}

auto system_state::operator/=(const system_state& other) -> system_state& {
    position_m.array() /= other.position_m.array();
    velocity_m_s.array() /= other.velocity_m_s.array();
    attitude.coeffs().array() /= other.attitude.coeffs().array();
    angular_velocity_m_s.array() /= other.angular_velocity_m_s.array();
    rod_magnetizations.array() /= other.rod_magnetizations.array();
    return *this;
}

auto system_state::altitude_m() const -> real {
    return position_m.norm() - earth_radius_m;
}

auto system_state::has_nan() const -> bool {
    return position_m.hasNaN() || velocity_m_s.hasNaN() || attitude.coeffs().hasNaN() || angular_velocity_m_s.hasNaN() || rod_magnetizations.hasNaN();
}

auto system_state::compute_attitude_derivative(const quat& q_att, const vec3& omega) -> vec4 {
    // dq/dt = 0.5 * q * omega_quat
    const quat omega_q(0, omega.x(), omega.y(), omega.z());
    return 0.5 * (q_att * omega_q).coeffs();
}

}  // namespace aos
