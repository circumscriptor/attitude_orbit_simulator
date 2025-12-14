#pragma once

#include "aos/core/types.hpp"

#include <numbers>

namespace aos {

class hysteresis_rod {
public:

    static constexpr double vacuum_permeability = 4.0 * std::numbers::pi * 1e-7;  // [H/m] or [T*m/A]

    // Stability thresholds
    static constexpr double epsilon_denominator = 1e-9;
    static constexpr double epsilon_langevin    = 1e-6;
    static constexpr double epsilon_dh_dt       = 1e-12;
    static constexpr double epsilon_vector      = 1e-12;

    struct ja_parameters {
        double ms;     // Saturation Magnetization [A/m]
        double a;      // Anhysteretic shape parameter [A/m]
        double k;      // Pinning energy density (coercivity) [A/m]
        double c;      // Reversibility coefficient [0-1]
        double alpha;  // Inter-domain coupling coefficient

        void debug_print() const;

        static ja_parameters hymu80();
    };

    /**
     * @brief Constructor for a hysteresis rod.
     * @param volume Volume of the rod [m^3].
     * @param orientation Unit vector defining the rod's axis in the body frame.
     * @param ja_params J-A model parameters.
     */
    hysteresis_rod(double volume, const vec3& orientation, const ja_parameters& ja_params);

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
    [[nodiscard]] vec3 magnetic_moment(double m_irr_am, const vec3& b_body_t) const;

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
    [[nodiscard]] double magnetization_derivative(double m_irr_am, const vec3& b_body_t, const vec3& b_dot_body_t) const;

    /**
     * @brief Calculates derivative based on scalar H-Field.
     *
     * @param m_irr_am Current irreversible magnetization [A/m].
     * @param h_along_rod H-Field intensity along the rod [A/m].
     * @param dh_dt Rate of change of H-field [A/m/s].
     */
    [[nodiscard]] double magnetization_derivative(double m_irr_am, double h_along_rod, double dh_dt) const;

protected:

    /**
     * @brief Computes the Anhysteretic Magnetization M_an(H_eff).
     * Uses the Langevin function: M_an = Ms * (coth(Heff/a) - a/Heff)
     */
    [[nodiscard]] double calculate_anhysteretic(double h_eff_am) const;

    /**
     * @brief Computes Effective Field H_eff = H + alpha * M.
     */
    [[nodiscard]] double calculate_h_eff(double h_along_rod, double m_val) const;

private:

    double        _volume;
    vec3          _orientation_body;
    ja_parameters _params;
};

}  // namespace aos
