#pragma once

#include "types.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resize.hpp>
#include <boost/numeric/odeint/util/same_size.hpp>
#include <boost/type_traits/integral_constant.hpp>

#include <algorithm>

namespace aos {

struct system_state {
    quat attitude;
    vec3 angular_velocity;
    vecX rod_magnetizations;

    system_state& operator+=(const system_state& other) {
        attitude.coeffs() += other.attitude.coeffs();
        angular_velocity += other.angular_velocity;
        rod_magnetizations += other.rod_magnetizations;
        return *this;
    }

    system_state& operator+=(double scalar) {
        attitude.coeffs().array() += scalar;
        angular_velocity.array() += scalar;
        rod_magnetizations.array() += scalar;
        return *this;
    }

    system_state& operator-=(const system_state& other) {
        attitude.coeffs() -= other.attitude.coeffs();
        angular_velocity -= other.angular_velocity;
        rod_magnetizations -= other.rod_magnetizations;
        return *this;
    }

    system_state& operator*=(double scalar) {
        attitude.coeffs() *= scalar;
        angular_velocity *= scalar;
        rod_magnetizations *= scalar;
        return *this;
    }

    system_state& operator/=(const system_state& other) {
        attitude.coeffs().cwiseQuotient(other.attitude.coeffs());
        angular_velocity.cwiseQuotient(other.angular_velocity);
        rod_magnetizations.cwiseQuotient(other.rod_magnetizations);
        return *this;
    }
};

inline system_state operator+(double scalar, const system_state& state) {
    return system_state{state} += scalar;
}

inline system_state operator+(const system_state& lhs, const system_state& rhs) {
    return system_state{lhs} += rhs;
}

inline system_state operator*(const system_state& state, double scalar) {
    return system_state{state} *= scalar;
}

inline system_state operator*(double scalar, const system_state& state) {
    return state * scalar;
}

inline system_state operator/(const system_state& lhs, const system_state& rhs) {
    return system_state{lhs} /= rhs;
}

inline aos::system_state abs(const aos::system_state& state) {
    aos::system_state result;
    result.rod_magnetizations.resize(state.rod_magnetizations.size());
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
    static bool same_size(const aos::system_state& s1, const aos::system_state& s2) { return s1.rod_magnetizations.size() == s2.rod_magnetizations.size(); }
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
    double operator()(const aos::system_state& s) const {
        return std::max({s.attitude.coeffs().cwiseAbs().maxCoeff(), s.angular_velocity.cwiseAbs().maxCoeff(), s.rod_magnetizations.cwiseAbs().maxCoeff()});
    }
};

}  // namespace boost::numeric::odeint
