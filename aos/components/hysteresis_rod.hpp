#pragma once

#include "aos/core/types.hpp"

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
    double magnetization_derivative(double m_scalar_am, const vec3& b_body_t, const vec3& omega_rad_s) const;

    [[nodiscard]]
    double magnetization_derivative_from_h(double m_scalar_am, double h_along_rod, double dh_dt) const;

private:

    double        m_volume;
    vec3          m_orientation_body;
    ja_parameters m_params;
};

}  // namespace aos::components
