#pragma once

#include "aos/core/types.hpp"

#include <cmath>
#include <cstdlib>

namespace aos::components {

using core::vec3;

class hysteresis_rod {
public:

    static constexpr double vacuum_permeability = 1.25663706212e-6;

    struct ja_parameters {
        double ms;     // Saturation Magnetization [A/m]
        double a;      // Anhysteretic shape parameter [A/m]
        double k;      // Pinning energy density [A/m]
        double c;      // Reversibility coefficient [0-1]
        double alpha;  // Inter-domain coupling
    };

    hysteresis_rod(double volume, const vec3& orientation, const ja_parameters& ja_params)
        : m_volume(volume), m_orientation_body(orientation.normalized()), m_params(ja_params) {}

    [[nodiscard]]
    vec3 magnetic_moment(double m_scalar_am) const {
        return m_scalar_am * m_volume * m_orientation_body;
    }

    // Implementation of the J-A model
    [[nodiscard]]
    double magnetization_derivative(double m_scalar_am, const vec3& b_body_t, const vec3& omega_rad_s) const {
        static constexpr double absolute_error      = 1e-6;
        static constexpr double denominator_epsilon = 1e-9;

        // Convert B [Tesla] to H [A/m]
        double h_along_rod = b_body_t.dot(m_orientation_body) / vacuum_permeability;

        // 1. Calculate Anhysteretic Magnetization (Man)
        double heff = h_along_rod + (m_params.alpha * m_scalar_am);
        double man  = 0.0;
        if (std::abs(heff / m_params.a) > absolute_error) {  // Avoid division by zero
            man = m_params.ms * (1.0 / std::tanh(heff / m_params.a) - m_params.a / heff);
        }

        // 2. Calculate dM/dH
        double dman_dheff = m_params.ms / m_params.a * (1.0 - 1.0 / std::pow(std::tanh(heff / m_params.a), 2) + 1.0 / std::pow(heff / m_params.a, 2));
        double dmrev_dh   = m_params.c * dman_dheff;

        // This term prevents unphysical jumps when Man - M changes sign
        double sign_term   = (h_along_rod > 0) ? 1.0 : -1.0;
        double denominator = (m_params.k * sign_term) - (m_params.alpha * (man - m_scalar_am));
        if (std::abs(denominator) < denominator_epsilon) {
            denominator = denominator_epsilon;
        }

        double dmirr_dh = (man - m_scalar_am) / denominator;
        double dm_dh    = dmrev_dh + dmirr_dh;

        // 3. Calculate dH/dt
        double dh_dt = (-omega_rad_s.cross(b_body_t)).dot(m_orientation_body) / vacuum_permeability;

        // 4. Return dM/dt using the chain rule
        return dm_dh * dh_dt;
    }

private:

    double        m_volume;
    vec3          m_orientation_body;
    ja_parameters m_params;
};

}  // namespace aos::components
