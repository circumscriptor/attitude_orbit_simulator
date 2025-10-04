#pragma once

#include "types.hpp"

namespace aos::core {

struct system_state {
    quat attitude;
    vec3 angular_velocity;
    vecX rod_magnetizations;

    system_state operator+(const system_state& rhs) const {
        system_state result;
        result.attitude.coeffs()  = attitude.coeffs() + rhs.attitude.coeffs();
        result.angular_velocity   = angular_velocity + rhs.angular_velocity;
        result.rod_magnetizations = rod_magnetizations + rhs.rod_magnetizations;
        return result;
    }

    system_state operator-(const system_state& rhs) const {
        system_state result;
        result.attitude.coeffs()  = attitude.coeffs() - rhs.attitude.coeffs();
        result.angular_velocity   = angular_velocity - rhs.angular_velocity;
        result.rod_magnetizations = rod_magnetizations - rhs.rod_magnetizations;
        return result;
    }

    system_state operator*(double scalar) const {
        system_state result;
        result.attitude.coeffs()  = attitude.coeffs() * scalar;
        result.angular_velocity   = angular_velocity * scalar;
        result.rod_magnetizations = rod_magnetizations * scalar;
        return result;
    }

    system_state& operator+=(const system_state& rhs) {
        attitude.coeffs() += rhs.attitude.coeffs();
        angular_velocity += rhs.angular_velocity;
        rod_magnetizations += rhs.rod_magnetizations;
        return *this;
    }

    system_state& operator-=(const system_state& rhs) {
        attitude.coeffs() -= rhs.attitude.coeffs();
        angular_velocity -= rhs.angular_velocity;
        rod_magnetizations -= rhs.rod_magnetizations;
        return *this;
    }

    system_state& operator*=(double scalar) {
        attitude.coeffs() *= scalar;
        angular_velocity *= scalar;
        rod_magnetizations *= scalar;
        return *this;
    }
};

}  // namespace aos::core
