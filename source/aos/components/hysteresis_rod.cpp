#include "hysteresis_rod.hpp"

#include "aos/core/constants.hpp"

#include <cmath>
#include <cstdlib>

double aos::components::hysteresis_rod::magnetization_derivative(double m_scalar_am, const vec3& b_body_t, const vec3& omega_rad_s) const {
    using aos::core::vacuum_permeability;

    static constexpr double absolute_error      = 1e-6;
    static constexpr double denominator_epsilon = 1e-9;

    // Convert B [Tesla] to H [A/m]
    const double h_along_rod = b_body_t.dot(m_orientation_body) / vacuum_permeability;

    // Calculate dH/dt first, as its sign is needed
    const double dh_dt = (-omega_rad_s.cross(b_body_t)).dot(m_orientation_body) / vacuum_permeability;

    // --- Jiles-Atherton Model Calculations ---

    // 1. Anhysteretic Magnetization (Man) and its derivative
    const double heff       = h_along_rod + (m_params.alpha * m_scalar_am);
    double       man        = 0.0;
    double       dman_dheff = 0.0;

    if (std::abs(heff / m_params.a) > absolute_error) {
        const double x     = heff / m_params.a;
        const double tanhx = std::tanh(x);
        man                = m_params.ms * (1.0 / tanhx - 1.0 / x);
        dman_dheff         = m_params.ms / m_params.a * (1.0 - 1.0 / (tanhx * tanhx) + 1.0 / (x * x));
    }

    // 2. Reversible and Irreversible components of dM/dH

    // CORRECTION 1: The sign term MUST be based on dH/dt, not H.
    const double sign_term = (dh_dt >= 0.0) ? 1.0 : -1.0;

    double denominator = (m_params.k * sign_term) - (m_params.alpha * (man - m_scalar_am));

    // CORRECTION 3: Use copysign for the stability check to preserve the sign.
    if (std::abs(denominator) < denominator_epsilon) {
        denominator = std::copysign(denominator_epsilon, denominator);
    }

    const double dmirr_dh = (man - m_scalar_am) / denominator;
    const double dmrev_dh = m_params.c * dman_dheff;

    // CORRECTION 2: The total dM/dH is a weighted sum of reversible and irreversible parts.
    const double dm_dh = ((1.0 - m_params.c) * dmirr_dh) + dmrev_dh;

    // 3. Return dM/dt using the chain rule
    return dm_dh * dh_dt;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double aos::components::hysteresis_rod::magnetization_derivative_from_h(double m_scalar_am, double h_along_rod, double dh_dt) const {
    static constexpr double absolute_error      = 1e-6;
    static constexpr double denominator_epsilon = 1e-9;

    // --- Jiles-Atherton Model Calculations ---

    // 1. Anhysteretic Magnetization (Man) and its derivative
    const double heff       = h_along_rod + (m_params.alpha * m_scalar_am);
    double       man        = 0.0;
    double       dman_dheff = 0.0;

    // CORRECTION 3: dman_dheff calculation is now protected.
    if (std::abs(heff / m_params.a) > absolute_error) {
        const double x     = heff / m_params.a;
        const double tanhx = std::tanh(x);
        man                = m_params.ms * (1.0 / tanhx - 1.0 / x);
        dman_dheff         = m_params.ms / m_params.a * (1.0 - 1.0 / (tanhx * tanhx) + 1.0 / (x * x));
    }

    // 2. Reversible and Irreversible components of dM/dH

    // This was already correct: sign is based on dH/dt
    const double sign_term = (dh_dt >= 0.0) ? 1.0 : -1.0;

    double denominator = (m_params.k * sign_term) - (m_params.alpha * (man - m_scalar_am));

    // CORRECTION 2: Use copysign for the stability check to preserve the sign.
    if (std::abs(denominator) < denominator_epsilon) {
        denominator = std::copysign(denominator_epsilon, denominator);
    }

    const double dmirr_dh = (man - m_scalar_am) / denominator;
    const double dmrev_dh = m_params.c * dman_dheff;

    // CORRECTION 1: The total dM/dH is a weighted sum.
    const double dm_dh = ((1.0 - m_params.c) * dmirr_dh) + dmrev_dh;

    // 3. Return dM/dt using the chain rule
    return dm_dh * dh_dt;
}
