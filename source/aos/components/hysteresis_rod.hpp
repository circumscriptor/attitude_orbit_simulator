#pragma once

#include "aos/core/types.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
#include <optional>
// clang-format on

namespace aos {

struct hysteresis_parameters {
    double ms;     // [A/m] Saturation Magnetization
    double a;      // [A/m] Anhysteretic shape parameter
    double k;      // [A/m] Pinning energy density (coercivity)
    double c;      // [-] Reversibility coefficient (0..1)
    double alpha;  // [-] Inter-domain coupling coefficient

    void from_toml(const toml::table& table);

    void debug_print() const;

    static auto hymu80() -> hysteresis_parameters;
};

struct hysteresis_rod_properties {
    double volume_m3;
    vec3   orientation;
    // [optional] custom hysteresis for this rod
    std::optional<hysteresis_parameters> hysteresis;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

class hysteresis_rod {
public:

    // Stability thresholds
    static constexpr double epsilon_langevin = 1e-6;
    static constexpr double epsilon_vector   = 1e-12;
    // Threshold below which dh/dt is treated as static
    static constexpr double epsilon_dh_dt = 1e-9;
    // Threshold for singularity (denominator -> 0)
    static constexpr double epsilon_denominator = 1e-9;
    // Tolerance for causality checks (preventing noise triggers)
    static constexpr double tolerance_causality = 1e-12;
    // Physical floor for 'k' to prevent division by zero in max_chi calculation
    static constexpr double min_k_value = 1e-3;

    hysteresis_rod(const hysteresis_rod_properties& properties, const hysteresis_parameters& params);
    explicit hysteresis_rod(const hysteresis_rod_properties& properties);

    [[nodiscard]] auto hysteresis() const noexcept -> const hysteresis_parameters&;

    /**
     * @brief Calculates the TOTAL magnetic dipole moment (Irreversible + Reversible).
     *
     * M_total = (1-c)*M_irr + c*M_an
     * Torque = M_total * Volume * Orientation X B_body
     *
     * @param m_irr_am The current state scalar (Irreversible Magnetization) [A/m].
     * @param b_body_t The magnetic field vector in the body frame [T].
     * @return Total magnetic dipole moment vector [A*m^2].
     */
    [[nodiscard]] auto magnetic_moment(double m_irr_am, const vec3& b_body_t) const -> vec3;

    /**
     * @brief Calculates the time derivative of the irreversible magnetization (dM_irr/dt).
     *
     * This solves the Jiles-Atherton differential equation:
     * dM_irr/dt = (M_an - M_irr) / (k*delta - alpha*(M_an - M_irr)) * dH/dt
     *
     * @param m_irr_am Current scalar irreversible magnetization [A/m].
     * @param b_body_t Current magnetic field in the body frame [T].
     * @param b_dot_body_t Rate of change of the magnetic field in the body frame [T/s].
     * @return The rate of change dM_irr/dt [A/m/s].
     */
    [[nodiscard]] auto magnetization_derivative(double m_irr_am, const vec3& b_body_t, const vec3& b_dot_body_t) const -> double;

    /**
     * @brief Calculates derivative based on scalar H-Field.
     *
     * @param m_irr_am Current irreversible magnetization [A/m].
     * @param h_along_rod H-Field intensity along the rod [A/m].
     * @param dh_dt Rate of change of H-field [A/m/s].
     */
    [[nodiscard]] auto magnetization_derivative(double m_irr_am, double h_along_rod, double dh_dt) const -> double;

protected:

    /**
     * @brief Computes the Anhysteretic Magnetization M_an(H_eff).
     * Uses the Langevin function: M_an = Ms * (coth(Heff/a) - a/Heff)
     */
    [[nodiscard]] auto calculate_anhysteretic(double h_eff_am) const -> double;

    /**
     * @brief Computes Effective Field H_eff = H + alpha * M.
     */
    [[nodiscard]] auto calculate_h_eff(double h_along_rod, double m_val) const -> double;

private:

    double                _volume;
    vec3                  _orientation_body;
    hysteresis_parameters _hysteresis;
};

}  // namespace aos
