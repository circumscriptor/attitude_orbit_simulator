#include "hysteresis_rod.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

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

hysteresis_rod::hysteresis_rod(double volume, const vec3& orientation, const ja_parameters& ja_params)
    : _volume(volume), _orientation_body(orientation), _params(ja_params) {
    if (orientation.norm() < vector_norm_epsilon) {
        throw std::runtime_error("Orientation vector must be non-zero.");
    }
    _orientation_body.normalize();

    if (_params.ms <= 0.0) {
        throw std::runtime_error("Saturation magnetization (ms) must be positive.");
    }
    if (_params.a <= 0.0) {
        throw std::runtime_error("Anhysteretic parameter (a) must be positive.");
    }
    if (_params.k <= 0.0) {
        throw std::runtime_error("Pinning energy (k) must be positive.");
    }
    if (_params.c < 0.0 || _params.c > 1.0) {
        throw std::runtime_error("Reversible fraction (c) must be in [0, 1].");
    }
    if (_params.alpha < 0.0) {
        throw std::runtime_error("Inter-domain coupling (alpha) must be non-negative.");
    }
    if (_volume <= 0.0) {
        throw std::runtime_error("Volume must be positive.");
    }
}

vec3 hysteresis_rod::magnetic_moment(double m_scalar_am) const {
    return m_scalar_am * _volume * _orientation_body;
}

double hysteresis_rod::effective_field(double h_along_rod, double m_clamped) const {
    return h_along_rod + (_params.alpha * m_clamped);
}

hysteresis_rod::anhysteretic hysteresis_rod::anhysteretic_magnetization(double h_eff) const {
    const double x = h_eff / _params.a;

    // Use Taylor series for small x to maintain numerical stability
    if (std::abs(x) < absolute_error) {
        // Man ≈ Ms * (x/3 - x³/45)
        // dMan/dH_eff ≈ Ms/a * (1/3 - x²/15)
        return {
            .man        = _params.ms * (x / 3.0 - (x * x * x) / 45.0),            // NOLINT(readability-magic-numbers)
            .dman_dheff = _params.ms / _params.a * (1.0 / 3.0 - (x * x) / 15.0),  // NOLINT(readability-magic-numbers)
        };
    }

    // Alternative (more numerically stable for large |x|):
    // const double cothx = (x > 0) ? (1.0 + 2.0 / (std::exp(2.0 * x) - 1.0)) : (1.0 - 2.0 / (std::exp(-2.0 * x) - 1.0));

    // Brillouin/Langevin function: Man = Ms * [coth(x) - 1/x]
    const double cothx = 1.0 / std::tanh(x);
    const double man   = _params.ms * (cothx - 1.0 / x);

    // Exact derivative: d/dx[coth(x) - 1/x] = -csch²(x) + 1/x²
    const double sinhx      = std::sinh(x);
    const double cschx_sq   = 1.0 / (sinhx * sinhx);
    const double dman_dheff = _params.ms / _params.a * (-cschx_sq + 1.0 / (x * x));

    return {
        .man        = man,
        .dman_dheff = dman_dheff,
    };
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double hysteresis_rod::irreversible_susceptibility(double man, double m_clamped, double dh_dt) const {
    const double m_diff = man - m_clamped;

    // Sign term based on dH/dt direction
    const double sign_term = (dh_dt >= 0.0) ? 1.0 : -1.0;

    // Calculate denominator with stability check
    const double denominator = (_params.k * sign_term) - (_params.alpha * m_diff);

    // Preserve sign while ensuring numerical stability
    if (std::abs(denominator) < denominator_epsilon) {
        return m_diff / std::copysign(denominator_epsilon, denominator);
    }
    return m_diff / denominator;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double hysteresis_rod::total_susceptibility(double dmirr_dh, double dman_dheff) const {
    const double dmrev_dh = _params.c * dman_dheff;

    // Total dM/dH is weighted sum of reversible and irreversible parts
    return ((1.0 - _params.c) * dmirr_dh) + dmrev_dh;
}

double hysteresis_rod::magnetization_derivative(double m_scalar_am, const vec3& b_body_t, const vec3& omega_rad_s) const {
    // Convert B [Tesla] to H [A/m]
    const double h_along_rod = b_body_t.dot(_orientation_body) / vacuum_permeability;

    // Calculate dH/dt first, as its sign is needed
    const double dh_dt = (-omega_rad_s.cross(b_body_t)).dot(_orientation_body) / vacuum_permeability;

    return magnetization_derivative_from_h(m_scalar_am, h_along_rod, dh_dt);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double hysteresis_rod::magnetization_derivative_from_h(double m_scalar_am, double h_along_rod, double dh_dt) const {
    if (std::abs(dh_dt) < dh_dt_threshold) {
        return 0.0;
    }

    const double m_clamped       = std::clamp(m_scalar_am, -_params.ms, _params.ms);
    const double heff            = effective_field(h_along_rod, m_clamped);
    const auto [man, dman_dheff] = anhysteretic_magnetization(heff);
    const double dmirr_dh        = irreversible_susceptibility(man, m_clamped, dh_dt);
    const double dm_dh           = total_susceptibility(dmirr_dh, dman_dheff);
    return dm_dh * dh_dt;  // Return dM/dt using the chain rule
}

}  // namespace aos::components
