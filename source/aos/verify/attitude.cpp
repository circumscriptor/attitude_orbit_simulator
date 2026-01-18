#include "attitude.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/environment/orbital_mechanics.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/dynamics.hpp"
#include "aos/verify/observers.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include <memory>
#include <string>

namespace aos {

void verify_attitude(const std::string& output_filename, const simulation_parameters& params) {
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::vector_space_algebra;

    auto satellite = std::make_shared<spacecraft>(params.satellite);
    auto env       = std::make_shared<environment_model>(params.simulation_year, params.gravity_model_degree);

    spacecraft_dynamics dynamics(satellite, env);
    attitude_observer   observer(output_filename);

    const auto [position, velocity] = orbital_converter::to_cartesian(params.orbit);

    system_state state;
    state.position         = position;
    state.velocity         = velocity;
    state.attitude         = quat::Identity();
    state.angular_velocity = vec3::Zero();

    using stepper_type = runge_kutta_dopri5<system_state, double, system_state, double, vector_space_algebra>;
    auto stepper       = make_controlled<stepper_type>(default_absolute_error, default_relative_error);
    integrate_adaptive(stepper, dynamics, state, params.t_start, params.t_end, params.dt_initial, observer);
}

}  // namespace aos
