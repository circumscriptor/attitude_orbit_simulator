#include "hysteresis_rod.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace aos::components {

void hysteresis_rod::ja_parameters::debug_print() const {
    std::cout << "-- hysteresis properties --"                           //
              << "\nhysteresis saturation magnetization:     " << ms     //
              << "\nhysteresis anhysteretic shape parameter: " << a      //
              << "\nhysteresis pinning energy density:       " << k      //
              << "\nhysteresis reversibility coefficient:    " << c      //
              << "\nhysteresis inter-domain coupling:        " << alpha  //
              << '\n';
}

hysteresis_rod::ja_parameters hysteresis_rod::ja_parameters::hymu80() {
    // NOLINTBEGIN(readability-magic-numbers)
    return {
        .ms    = 6.0e5,   // Saturation Magnetization [A/m]
        .a     = 6.5,     // Anhysteretic shape parameter [A/m]
        .k     = 4.0,     // Pinning energy density (coercivity) [A/m]
        .c     = 0.05,    // Reversibility coefficient
        .alpha = 1.0e-5,  // Inter-domain coupling
    };
    // NOLINTEND(readability-magic-numbers)
}

double hysteresis_rod::magnetization_derivative(double m_scalar_am, const vec3& b_body_t, const vec3& omega_rad_s) const {
    static constexpr double absolute_error      = 1e-6;
    static constexpr double denominator_epsilon = 1e-9;

    // Convert B [Tesla] to H [A/m]
    const double h_along_rod = b_body_t.dot(_orientation_body) / vacuum_permeability;

    // Calculate dH/dt first, as its sign is needed
    const double dh_dt = (-omega_rad_s.cross(b_body_t)).dot(_orientation_body) / vacuum_permeability;

    // --- Jiles-Atherton Model Calculations ---

    // 1. Anhysteretic Magnetization (Man) and its derivative
    const double heff       = h_along_rod + (_params.alpha * m_scalar_am);
    double       man        = 0.0;
    double       dman_dheff = 0.0;

    if (std::abs(heff / _params.a) > absolute_error) {
        const double x     = heff / _params.a;
        const double tanhx = std::tanh(x);
        man                = _params.ms * (1.0 / tanhx - 1.0 / x);
        dman_dheff         = _params.ms / _params.a * (1.0 - 1.0 / (tanhx * tanhx) + 1.0 / (x * x));
    }

    // 2. Reversible and Irreversible components of dM/dH

    // CORRECTION 1: The sign term MUST be based on dH/dt, not H.
    const double sign_term = (dh_dt >= 0.0) ? 1.0 : -1.0;

    double denominator = (_params.k * sign_term) - (_params.alpha * (man - m_scalar_am));

    // CORRECTION 3: Use copysign for the stability check to preserve the sign.
    if (std::abs(denominator) < denominator_epsilon) {
        denominator = std::copysign(denominator_epsilon, denominator);
    }

    const double dmirr_dh = (man - m_scalar_am) / denominator;
    const double dmrev_dh = _params.c * dman_dheff;

    // CORRECTION 2: The total dM/dH is a weighted sum of reversible and irreversible parts.
    const double dm_dh = ((1.0 - _params.c) * dmirr_dh) + dmrev_dh;

    // 3. Return dM/dt using the chain rule
    return dm_dh * dh_dt;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double hysteresis_rod::magnetization_derivative_from_h(double m_scalar_am, double h_along_rod, double dh_dt) const {
    static constexpr double absolute_error      = 1e-6;
    static constexpr double denominator_epsilon = 1e-9;

    // --- Jiles-Atherton Model Calculations ---

    // 1. Anhysteretic Magnetization (Man) and its derivative
    const double heff       = h_along_rod + (_params.alpha * m_scalar_am);
    double       man        = 0.0;
    double       dman_dheff = 0.0;

    // CORRECTION 3: dman_dheff calculation is now protected.
    if (std::abs(heff / _params.a) > absolute_error) {
        const double x     = heff / _params.a;
        const double tanhx = std::tanh(x);
        man                = _params.ms * (1.0 / tanhx - 1.0 / x);
        dman_dheff         = _params.ms / _params.a * (1.0 - 1.0 / (tanhx * tanhx) + 1.0 / (x * x));
    }

    // 2. Reversible and Irreversible components of dM/dH

    // This was already correct: sign is based on dH/dt
    const double sign_term = (dh_dt >= 0.0) ? 1.0 : -1.0;

    double denominator = (_params.k * sign_term) - (_params.alpha * (man - m_scalar_am));

    // CORRECTION 2: Use copysign for the stability check to preserve the sign.
    if (std::abs(denominator) < denominator_epsilon) {
        denominator = std::copysign(denominator_epsilon, denominator);
    }

    const double dmirr_dh = (man - m_scalar_am) / denominator;
    const double dmrev_dh = _params.c * dman_dheff;

    // CORRECTION 1: The total dM/dH is a weighted sum.
    const double dm_dh = ((1.0 - _params.c) * dmirr_dh) + dmrev_dh;

    // 3. Return dM/dt using the chain rule
    return dm_dh * dh_dt;
}

}  // namespace aos::components
