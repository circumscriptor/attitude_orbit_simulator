#include "hysteresis_rod.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/toml.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace aos {

void hysteresis_parameters::from_toml(const toml_table& table) {
    ms    = table["ms"].value_or(0.0);
    a     = table["a"].value_or(0.0);
    k     = table["k"].value_or(0.0);
    c     = table["c"].value_or(0.0);
    alpha = table["alpha"].value_or(0.0);
}

void hysteresis_parameters::debug_print() const {
    std::cout << "-- hysteresis properties --"     //
              << "\n  Ms (Saturation): " << ms     //
              << "\n  a (Shape):       " << a      //
              << "\n  k (Coercivity):  " << k      //
              << "\n  c (Reversible):  " << c      //
              << "\n  alpha (Coupling):" << alpha  //
              << '\n';
}

void hysteresis_rod_properties::from_toml(const toml_table& table) {
    volume_m3 = table["volume"].value_or(0.0);

    if (const auto* vec = table["orientation"].as_array()) {
        orientation <<                   //
            vec->get(0)->value_or(0.0),  //
            vec->get(1)->value_or(0.0),  //
            vec->get(2)->value_or(0.0);
    }

    if (const auto* hyst = table["hysteresis"].as_table()) {
        hysteresis.emplace().from_toml(*hyst);
    }
}

void hysteresis_rod_properties::debug_print() const {
    std::cout << "-- hysteresis rod properties --"                                                           //
              << "\n  volume:      " << volume_m3                                                            //
              << "\n  orientation: " << orientation.x() << ' ' << orientation.y() << ' ' << orientation.z()  //
              << '\n';

    if (hysteresis) {
        hysteresis->debug_print();
    }
}

auto hysteresis_parameters::hymu80() -> hysteresis_parameters {
    // NOLINTBEGIN(readability-magic-numbers)
    return {
        .ms    = 6.0e5,  // ~0.75 Tesla / mu0
        .a     = 6.5,
        .k     = 4.0,
        .c     = 0.05,
        .alpha = 1.0e-5,
    };
    // NOLINTEND(readability-magic-numbers)
}

hysteresis_rod::hysteresis_rod(const hysteresis_rod_properties& properties, const hysteresis_parameters& params)
    : _volume(properties.volume_m3),
      _orientation_body(properties.orientation.normalized()),
      _hysteresis(properties.hysteresis ? *properties.hysteresis : params) {
    if (_volume <= 0.0) {
        throw std::runtime_error("Volume must be positive");
    }
    if (_hysteresis.ms <= 0.0) {
        throw std::runtime_error("Ms must be positive");
    }
    if (_hysteresis.a <= 0.0) {
        throw std::runtime_error("Parameter 'a' must be positive");
    }
    if (_hysteresis.k <= 0.0) {
        throw std::runtime_error("Parameter 'k' must be positive");
    }
    if (_hysteresis.c < 0.0 || _hysteresis.c > 1.0) {
        throw std::runtime_error("Parameter 'c' must be [0, 1]");
    }
}

hysteresis_rod::hysteresis_rod(const hysteresis_rod_properties& properties) : hysteresis_rod(properties, *properties.hysteresis) {}

auto hysteresis_rod::hysteresis() const noexcept -> const hysteresis_parameters& {
    return _hysteresis;
}

auto hysteresis_rod::calculate_h_eff(real h_along_rod, real m_val) const -> real {
    // H_eff = H + alpha * M
    return h_along_rod + (_hysteresis.alpha * m_val);
}

auto hysteresis_rod::calculate_anhysteretic(real h_eff_am) const -> real {
    const real ratio = h_eff_am / _hysteresis.a;

    // numerical stability for langevin function near zero
    if (std::abs(ratio) < epsilon_langevin) {
        // taylor expansion: L(x) approx x/3 - x^3/45
        return _hysteresis.ms * (ratio / 3.0);  // NOLINT(readability-magic-numbers)
        // - (ratio*ratio*ratio)/45.0;
    }

    // langevin: L(x) = coth(x) - 1/x
    // M_an = Ms * L(x)
    return _hysteresis.ms * ((1.0 / std::tanh(ratio)) - (1.0 / ratio));
}

auto hysteresis_rod::magnetic_moment(real m_irr_am, const vec3& b_body_t) const -> vec3 {
    // get H field along the rod
    const real h_applied = b_body_t.dot(_orientation_body) / vacuum_permeability;

    // clamp M_irr to physical limits to prevent divergence in H_eff
    const real m_irr_clamped = std::clamp(m_irr_am, -_hysteresis.ms, _hysteresis.ms);

    // calculate H_eff based on current M_irr
    const real h_eff = calculate_h_eff(h_applied, m_irr_clamped);
    const real m_an  = calculate_anhysteretic(h_eff);

    // M_tot = M_irr + M_rev
    // M_rev = c * (M_an - M_irr)
    // => M_tot = (1 - c) * M_irr + c * M_an
    const real m_total = ((1.0 - _hysteresis.c) * m_irr_clamped) + (_hysteresis.c * m_an);

    // moment: m = M_tot * volume * direction
    return m_total * _volume * _orientation_body;
}

auto hysteresis_rod::magnetization_derivative(real m_irr_am, const vec3& b_body_t, const vec3& b_dot_body_t) const -> real {
    // calculate H and dH/dt along the rod
    const real h_applied = b_body_t.dot(_orientation_body) / vacuum_permeability;
    const real dh_dt     = b_dot_body_t.dot(_orientation_body) / vacuum_permeability;
    return magnetization_derivative(m_irr_am, h_applied, dh_dt);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto hysteresis_rod::magnetization_derivative(real m_irr_am, real h_along_rod, real dh_dt) const -> real {
    // If saturated and driving further into saturation, no change possible.
    if (m_irr_am >= _hysteresis.ms && dh_dt > 0.0) {
        return 0.0;
    }
    if (m_irr_am <= -_hysteresis.ms && dh_dt < 0.0) {
        return 0.0;
    }

    // Skip calculation if field change is negligible (Static field)
    if (std::abs(dh_dt) < epsilon_dh_dt) {
        return 0.0;
    }

    const real m_irr_clamped = std::clamp(m_irr_am, -_hysteresis.ms, _hysteresis.ms);
    const real h_eff         = calculate_h_eff(h_along_rod, m_irr_clamped);
    const real m_an          = calculate_anhysteretic(h_eff);
    const real delta         = (dh_dt > 0.0) ? 1.0 : -1.0;
    const real numerator     = m_an - m_irr_clamped;
    const real denominator   = (_hysteresis.k * delta) - (_hysteresis.alpha * numerator);
    const real max_chi       = _hysteresis.ms / std::max(_hysteresis.k, min_k_value);

    real dmirr_dh = 0.0;
    if (std::abs(denominator) < epsilon_denominator) {
        // Singularity Handling (0/0 check)
        if (std::abs(numerator) < epsilon_denominator) {
            dmirr_dh = 0.0;
        } else {
            // Cap to physical max limit, preserving direction
            dmirr_dh = std::copysign(max_chi, numerator);
        }
    } else {
        // Standard J-A Equation
        dmirr_dh = numerator / denominator;

        // Numerical Stability Clamp (Sign-Preserving)
        if (std::abs(dmirr_dh) > max_chi) {
            dmirr_dh = std::copysign(max_chi, dmirr_dh);
        }
    }

    const real dm_irr_dt = dmirr_dh * dh_dt;

    // Enforce that magnetization changes in the direction driven by the field.
    // If H increases (dh_dt > 0), M_irr should not decrease.
    // If H decreases (dh_dt < 0), M_irr should not increase.
    // This prevents unphysical "active" behavior (energy generation) which
    // can occur numerically in certain parameter regimes.

    if (dh_dt > 0.0 && dm_irr_dt < -tolerance_causality) {
        return 0.0;  // Field Increasing, but Math says M decreases -> Clamp
    }

    if (dh_dt < 0.0 && dm_irr_dt > tolerance_causality) {
        return 0.0;  // Field Decreasing, but Math says M increases -> Clamp
    }

    return dm_irr_dt;
}

}  // namespace aos
