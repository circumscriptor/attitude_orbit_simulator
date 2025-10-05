#pragma once

#include "aos/core/types.hpp"

namespace aos::components {

class hysteresis_rod {
public:

    struct ja_parameters {
        double ms;     // Saturation Magnetization [A/m]
        double a;      // Anhysteretic shape parameter [A/m]
        double k;      // Pinning energy density [A/m]
        double c;      // Reversibility coefficient [0-1]
        double alpha;  // Inter-domain coupling

        void debug_print() const;

        static ja_parameters hymu80();
    };

    hysteresis_rod(double volume, const vec3& orientation, const ja_parameters& ja_params)
        : _volume(volume), _orientation_body(orientation.normalized()), _params(ja_params) {}

    [[nodiscard]]
    vec3 magnetic_moment(double m_scalar_am) const {
        return m_scalar_am * _volume * _orientation_body;
    }

    // Implementation of the J-A model
    [[nodiscard]]
    double magnetization_derivative(double m_scalar_am, const vec3& b_body_t, const vec3& omega_rad_s) const;

    [[nodiscard]]
    double magnetization_derivative_from_h(double m_scalar_am, double h_along_rod, double dh_dt) const;

private:

    double        _volume;
    vec3          _orientation_body;
    ja_parameters _params;
};

}  // namespace aos::components
