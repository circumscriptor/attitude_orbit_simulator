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

namespace aos {

struct system_state {
    vec3 position;
    vec3 velocity;
    quat attitude;
    vec3 angular_velocity;
    vecX rod_magnetizations;

    auto operator+=(const system_state& other) -> system_state& {
        position += other.position;
        velocity += other.velocity;
        attitude.coeffs() += other.attitude.coeffs();
        angular_velocity += other.angular_velocity;
        rod_magnetizations += other.rod_magnetizations;
        return *this;
    }

    auto operator+=(double scalar) -> system_state& {
        position.array() += scalar;
        velocity.array() += scalar;
        attitude.coeffs().array() += scalar;
        angular_velocity.array() += scalar;
        rod_magnetizations.array() += scalar;
        return *this;
    }

    auto operator-=(const system_state& other) -> system_state& {
        position -= other.position;
        velocity -= other.velocity;
        attitude.coeffs() -= other.attitude.coeffs();
        angular_velocity -= other.angular_velocity;
        rod_magnetizations -= other.rod_magnetizations;
        return *this;
    }

    auto operator*=(double scalar) -> system_state& {
        position *= scalar;
        velocity *= scalar;
        attitude.coeffs() *= scalar;
        angular_velocity *= scalar;
        rod_magnetizations *= scalar;
        return *this;
    }

    auto operator/=(const system_state& other) -> system_state& {
        position.cwiseQuotient(other.position);
        velocity.cwiseQuotient(other.velocity);
        attitude.coeffs().cwiseQuotient(other.attitude.coeffs());
        angular_velocity.cwiseQuotient(other.angular_velocity);
        rod_magnetizations.cwiseQuotient(other.rod_magnetizations);
        return *this;
    }
};

inline auto operator+(double scalar, const system_state& state) -> system_state {
    return system_state{state} += scalar;
}

inline auto operator+(const system_state& lhs, const system_state& rhs) -> system_state {
    return system_state{lhs} += rhs;
}

inline auto operator*(const system_state& state, double scalar) -> system_state {
    return system_state{state} *= scalar;
}

inline auto operator*(double scalar, const system_state& state) -> system_state {
    return state * scalar;
}

inline auto operator/(const system_state& lhs, const system_state& rhs) -> system_state {
    return system_state{lhs} /= rhs;
}

inline auto abs(const aos::system_state& state) -> aos::system_state {
    aos::system_state result;
    result.rod_magnetizations.resize(state.rod_magnetizations.size());
    result.position           = state.position.cwiseAbs();
    result.velocity           = state.velocity.cwiseAbs();
    result.attitude.coeffs()  = state.attitude.coeffs().cwiseAbs();
    result.angular_velocity   = state.angular_velocity.cwiseAbs();
    result.rod_magnetizations = state.rod_magnetizations.cwiseAbs();
    return result;
}

}  // namespace aos

namespace boost::numeric::odeint {

// Tell odeint that our struct is a valid state type
template <>
struct is_resizeable<aos::system_state> : boost::false_type {};

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
    using result_type = double;
    auto operator()(const aos::system_state& s) const -> double {
        return std::max({
            s.position.cwiseAbs().maxCoeff(),
            s.velocity.cwiseAbs().maxCoeff(),
            s.attitude.coeffs().cwiseAbs().maxCoeff(),
            s.angular_velocity.cwiseAbs().maxCoeff(),
            s.rod_magnetizations.size() > 0 ? s.rod_magnetizations.cwiseAbs().maxCoeff() : std::numeric_limits<double>::min(),
        });
    }
};

}  // namespace boost::numeric::odeint
