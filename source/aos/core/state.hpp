#pragma once

#include "types.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resize.hpp>
#include <boost/numeric/odeint/util/same_size.hpp>
#include <boost/type_traits/integral_constant.hpp>

#include <algorithm>
#include <limits>
// #include <type_traits>

namespace aos {

struct system_state {
    vec3 position_m;            //!< [m] Spacecraft position in ECI
    vec3 velocity_m_s;          //!< [m/s] Spacecraft velocity in ECI
    quat attitude;              //!< [-] Spacecraft orientation (body-to-ECI)
    vec3 angular_velocity_m_s;  //!< [m/s] Spacecraft angular velocity in body frame
    vecX rod_magnetizations;    //!< [?] Hysteresis rod magnetization state (array)

    auto operator+=(const system_state& other) -> system_state&;
    auto operator+=(real_t scalar) -> system_state&;
    auto operator-=(const system_state& other) -> system_state&;
    auto operator*=(real_t scalar) -> system_state&;
    auto operator/=(const system_state& other) -> system_state&;

    [[nodiscard]] auto altitude_m() const -> real_t;
    [[nodiscard]] auto has_nan() const -> bool;

    // quaternion derivative: 0.5 * q * omega
    [[nodiscard]] static auto compute_attitude_derivative(const quat& q_att, const vec3& omega) -> vec4;
};

inline auto operator+(real_t scalar, const system_state& state) -> system_state {
    return system_state{state} += scalar;
}

inline auto operator+(const system_state& lhs, const system_state& rhs) -> system_state {
    return system_state{lhs} += rhs;
}

inline auto operator*(const system_state& state, real_t scalar) -> system_state {
    return system_state{state} *= scalar;
}

inline auto operator*(real_t scalar, const system_state& state) -> system_state {
    return state * scalar;
}

inline auto operator/(const system_state& lhs, const system_state& rhs) -> system_state {
    return system_state{lhs} /= rhs;
}

inline auto abs(const aos::system_state& state) -> aos::system_state {
    aos::system_state result;
    result.rod_magnetizations.resize(state.rod_magnetizations.size());
    result.position_m           = state.position_m.cwiseAbs();
    result.velocity_m_s         = state.velocity_m_s.cwiseAbs();
    result.attitude.coeffs()    = state.attitude.coeffs().cwiseAbs();
    result.angular_velocity_m_s = state.angular_velocity_m_s.cwiseAbs();
    result.rod_magnetizations   = state.rod_magnetizations.cwiseAbs();
    return result;
}

}  // namespace aos

namespace boost::numeric::odeint {

// Tell odeint that our struct is a valid state type
template <>
struct is_resizeable<aos::system_state> : boost::true_type {};  // when switching to boost 1.90, use std::true_type

// Define how to check if two states have the same size
template <>
struct same_size_impl<aos::system_state, aos::system_state> {
    static auto same_size(const aos::system_state& s1, const aos::system_state& s2) -> bool {
        return s1.rod_magnetizations.size() == s2.rod_magnetizations.size();
    }
};

// Define how to resize one state to match another
template <>
struct resize_impl<aos::system_state, aos::system_state> {
    static void resize(aos::system_state& s1, const aos::system_state& s2) { s1.rod_magnetizations.resize(s2.rod_magnetizations.size()); }
};

// Define the "norm" used by adaptive steppers for error control
template <>
struct vector_space_norm_inf<aos::system_state> {
    using result_type = aos::real_t;
    auto operator()(const aos::system_state& s) const -> aos::real_t {
        return std::max({
            s.position_m.cwiseAbs().maxCoeff(),
            s.velocity_m_s.cwiseAbs().maxCoeff(),
            s.attitude.coeffs().cwiseAbs().maxCoeff(),
            s.angular_velocity_m_s.cwiseAbs().maxCoeff(),
            s.rod_magnetizations.size() > 0 ? s.rod_magnetizations.cwiseAbs().maxCoeff() : std::numeric_limits<aos::real_t>::min(),
        });
    }
};

}  // namespace boost::numeric::odeint
