#pragma once

#include "aos/core/types.hpp"

namespace aos {

class hysteresis_rod {
public:

    static constexpr double absolute_error      = 1e-6;
    static constexpr double denominator_epsilon = 1e-9;
    static constexpr double dh_dt_threshold     = 1e-12;
    static constexpr double vector_norm_epsilon = 1e-12;

    struct ja_parameters {
        double ms;     // Saturation Magnetization [A/m]
        double a;      // Anhysteretic shape parameter [A/m]
        double k;      // Pinning energy density [A/m]
        double c;      // Reversibility coefficient [0-1]
        double alpha;  // Inter-domain coupling

        void debug_print() const;

        static ja_parameters hymu80();
    };

    struct anhysteretic {
        double man;         // Anhysteretic Magnetization
        double dman_dheff;  // Derivative dMan/dHeff
    };

    hysteresis_rod(double volume, const vec3& orientation, const ja_parameters& ja_params);

    [[nodiscard]]
    vec3 magnetic_moment(double m_scalar_am) const;

    [[nodiscard]]
    double effective_field(double h_along_rod, double m_clamped) const;

    [[nodiscard]]
    anhysteretic anhysteretic_magnetization(double h_eff) const;

    [[nodiscard]]
    double irreversible_susceptibility(double man, double m_clamped, double dh_dt) const;

    [[nodiscard]]
    double total_susceptibility(double dmirr_dh, double dman_dheff) const;

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

}  // namespace aos
