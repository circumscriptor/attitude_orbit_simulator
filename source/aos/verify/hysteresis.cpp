#include "hysteresis.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/constants.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>

namespace aos::verify {

void hysteresis_loop_dynamics::operator()(const hysteresis_state_type& m, hysteresis_state_type& dm_dt, double t) const {
    // Applied H-field and its time derivative
    const double h     = h_max * sin(2.0 * std::numbers::pi * frequency * t);
    const double dh_dt = h_max * 2.0 * std::numbers::pi * frequency * cos(2.0 * std::numbers::pi * frequency * t);
    dm_dt              = _rod->magnetization_derivative_from_h(m, h, dh_dt);
}

void bh_observer::operator()(const hysteresis_state_type& m, double t) const {
    const double h = hysteresis_loop_dynamics::h_max * sin(2.0 * std::numbers::pi * hysteresis_loop_dynamics::frequency * t);
    // Calculate B = μ₀ * (H + M)
    const double b = vacuum_permeability * (h + m);
    (*output) << t << "," << h << "," << m << "," << b << "\n";
}

int verify_hysteresis(const std::string& output_filename, const hysteresis_rod::ja_parameters& params) {
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;

    static constexpr double epsilon = 1e-6;

    std::cout << "Starting HyMu-80 B-H Curve Verification..." << '\n';

    // --- Test Configuration ---
    static constexpr const double t_end = 2.0;    // [s] Simulate for 2 full cycles.
    static constexpr const double dt    = 0.001;  // [s] Initial time step.

    // 1. Instantiate the HyMu-80 Rod
    hysteresis_rod hymu80_rod(1.0, {1, 0, 0}, params);  // Volume and orientation don't matter here

    // 2. Setup the test
    hysteresis_loop_dynamics dynamics(hymu80_rod);

    std::ofstream data_file(output_filename);
    data_file << "time,H_Am,M_Am,B_T\n";

    bh_observer observer(data_file);

    // 3. Set initial conditions (demagnetized)
    hysteresis_state_type m_initial = {0.0};

    // 4. Run the integrator
    using stepper_type = runge_kutta_dopri5<hysteresis_state_type>;
    auto stepper       = make_controlled<stepper_type>(epsilon, epsilon, stepper_type{});
    integrate_adaptive(stepper, dynamics, m_initial, 0.0, t_end, dt, observer);

    std::cout << "Verification finished. Data saved to " << output_filename << '\n';
    return 0;
}

}  // namespace aos::verify
