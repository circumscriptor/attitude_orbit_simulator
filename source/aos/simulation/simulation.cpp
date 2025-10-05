#include "simulation.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/dynamics.hpp"
#include "aos/simulation/observer.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include <memory>
#include <string>

namespace aos::simulation {

void run_simulation(const std::string& output_filename, const simulation_parameters& params) {
    using aos::abs;
    using aos::components::spacecraft;
    using aos::environment::wmm2020_environment;
    using aos::simulation::csv_state_observer;
    using aos::simulation::simulation_parameters;
    using aos::simulation::spacecraft_dynamics;
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::vector_space_algebra;

    auto satellite   = std::make_shared<spacecraft>(params.spacecraft);
    auto environment = std::make_shared<wmm2020_environment>(params.environment);

    spacecraft_dynamics dynamics{satellite, environment};
    csv_state_observer  observer(output_filename, satellite->rods().size());
    system_state        initial;
    initial.attitude         = aos::quat::Identity();
    initial.angular_velocity = params.initial_angular_velocity;
    initial.rod_magnetizations.resize(4);
    initial.rod_magnetizations.setZero();

    using stepper_type = runge_kutta_dopri5<system_state, double, system_state, double, vector_space_algebra>;
    auto stepper       = make_controlled<stepper_type>(epsilon, epsilon);
    integrate_adaptive(stepper, dynamics, initial, params.t_start, params.t_end, params.dt_initial, observer);
}

}  // namespace aos::simulation
