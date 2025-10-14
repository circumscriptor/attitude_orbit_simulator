#include "hysteresis.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/constants.hpp"
#include "aos/verify/bh_observer.hpp"
#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include <string>

namespace aos {

void verify_hysteresis(const std::string& output_filename, const hysteresis_rod::ja_parameters& params) {
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;

    static constexpr const double t_start = 0.0;
    static constexpr const double t_end   = 2.0;    // [s] Simulate for 2 full cycles.
    static constexpr const double dt      = 0.001;  // [s] Initial time step.

    const hysteresis_rod           hymu80_rod(1.0, {1, 0, 0}, params);  // Volume and orientation don't matter here
    const hysteresis_loop_dynamics dynamics(hymu80_rod);
    const bh_observer              observer(output_filename);

    hysteresis_state_type m_initial{0.0};

    using stepper_type = runge_kutta_dopri5<hysteresis_state_type>;
    auto stepper       = make_controlled<stepper_type>(epsilon, epsilon, stepper_type{});
    integrate_adaptive(stepper, dynamics, m_initial, t_start, t_end, dt, observer);
}

}  // namespace aos
