#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/dynamics.hpp"
#include "aos/simulation/observer.hpp"

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/numeric/odeint.hpp>
#include <iostream>

// ... (init_logging function) ...

int main(int argc, char* argv[]) {
    init_logging();

    try {
        // Use the namespace for all project types
        using namespace aos;

        simulation::SimulationParameters params;
        if (!simulation::load_parameters(argc, argv, params)) {
            return 1;
        }

        BOOST_LOG_TRIVIAL(info) << "Configuration loaded successfully.";

        // Initialize models
        environment::Environment env(params);

        components::Spacecraft my_satellite(params.inertia_tensor);
        // ... (initialize satellite components with HyMu-80 as before) ...

        simulation::SpacecraftDynamics dynamics(my_satellite, env);

        // ... (set initial conditions as before) ...

        simulation::CsvStateObserver observer(params.output_filename);

        BOOST_LOG_TRIVIAL(info) << "Starting simulation...";

        // --- odeint integration ---
        // The stepper type definition is now more explicit
        using stepper_type =
            boost::numeric::odeint::runge_kutta_dopri5<core::SystemState, double, core::SystemState, double, boost::numeric::odeint::vector_space_algebra>;

        auto stepper = boost::numeric::odeint::make_controlled<stepper_type>(1.0e-6, 1.0e-6);

        boost::numeric::odeint::integrate_adaptive(stepper, dynamics, initial_state, params.t_start, params.t_end, params.dt_initial, observer);

        BOOST_LOG_TRIVIAL(info) << "Simulation finished successfully.";

    } catch (const std::exception& e) {
        BOOST_LOG_TRIVIAL(fatal) << "An unrecoverable error occurred: " << e.what();
        return 1;
    }

    return 0;
}
