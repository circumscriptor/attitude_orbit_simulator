// #include "aos/components/hysteresis_rod.hpp"
// #include "aos/components/permanent_magnet.hpp"
// #include "aos/components/spacecraft.hpp"
// #include "aos/core/state.hpp"
// #include "aos/environment/environment.hpp"
// #include "aos/simulation/config.hpp"
// #include "aos/simulation/dynamics.hpp"
// #include "aos/simulation/observer.hpp"

// #include <boost/log/core.hpp>
// #include <boost/log/expressions.hpp>
// #include <boost/log/trivial.hpp>
// #include <boost/log/utility/setup/common_attributes.hpp>
// #include <boost/log/utility/setup/file.hpp>
// #include <boost/numeric/odeint.hpp>
// #include <iostream>

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/dynamics.hpp"
#include "aos/simulation/observer.hpp"
// #include "aos/verify/hysteresis.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>

#include <memory>
#include <numbers>

int main(int /*argc*/, char** /*argv*/) {
    // init_logging();

    [[maybe_unused]] const aos::components::hysteresis_rod::ja_parameters hymu80_params{
        .ms    = 6.0e5,   // Saturation Magnetization [A/m]
        .a     = 6.5,     // Anhysteretic shape parameter [A/m]
        .k     = 4.0,     // Pinning energy density (coercivity) [A/m]
        .c     = 0.05,    // Reversibility coefficient
        .alpha = 1.0e-5,  // Inter-domain coupling
    };

    // try {
    //     // Use the namespace for all project types
    //     using namespace aos;

    //     simulation::SimulationParameters params;
    //     if (!simulation::load_parameters(argc, argv, params)) {
    //         return 1;
    //     }

    //     BOOST_LOG_TRIVIAL(info) << "Configuration loaded successfully.";

    //     // Initialize models
    //     environment::Environment env(params);

    //     components::Spacecraft my_satellite(params.inertia_tensor);
    //     // ... (initialize satellite components with HyMu-80 as before) ...

    //     simulation::SpacecraftDynamics dynamics(my_satellite, env);

    //     // ... (set initial conditions as before) ...

    //     simulation::CsvStateObserver observer(params.output_filename);

    //     BOOST_LOG_TRIVIAL(info) << "Starting simulation...";

    //     // --- odeint integration ---
    //     // The stepper type definition is now more explicit
    //     using stepper_type =
    //         boost::numeric::odeint::runge_kutta_dopri5<core::SystemState, double, core::SystemState, double, boost::numeric::odeint::vector_space_algebra>;

    //     auto stepper = boost::numeric::odeint::make_controlled<stepper_type>(1.0e-6, 1.0e-6);

    //     boost::numeric::odeint::integrate_adaptive(stepper, dynamics, initial_state, params.t_start, params.t_end, params.dt_initial, observer);

    //     BOOST_LOG_TRIVIAL(info) << "Simulation finished successfully.";

    // } catch (const std::exception& e) {
    //     BOOST_LOG_TRIVIAL(fatal) << "An unrecoverable error occurred: " << e.what();
    //     return 1;
    // }

    using aos::system_state;
    using aos::components::spacecraft;
    using aos::environment::wmm2020_environment;
    using aos::simulation::csv_state_observer;
    using aos::simulation::simulation_parameters;
    using aos::simulation::spacecraft_dynamics;
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::vector_space_algebra;

    const double epsilon       = 1e-6;
    const double mass_kg       = 1.3;
    const double mass_g        = mass_kg * 1000.;
    const double size_cm       = 10.;
    const double size_m        = size_cm / 100.;
    const double time          = 2. * 7. * 24. * 60. * 60.;
    const double rod_radius_cm = 0.5;
    const double rod_radius_m  = rod_radius_cm / 100.;
    const double rod_length_cm = 10.;
    const double rod_length_m  = rod_length_cm / 100.;
    const double rod_volume_m3 = rod_radius_m * rod_radius_m * std::numbers::pi * rod_length_m;

    // {{"N35", 1.21},  // Using nominal Br in Tesla
    //  {"N42", 1.32},
    //  {"N52", 1.45},
    //  {"N35SH", 1.19}};

    simulation_parameters params;
    params.dt_initial = 0.1;
    // params.initial_angular_velocity = {1.0, -0.5, 0.8};  // A significant initial tumble [rad/s]
    params.initial_angular_velocity = {0.1, -0.05, 0.08};
    params.orbit_altitude_km        = 6371. + 650.;
    params.orbit_inclination_deg    = 51.;
    params.t_start                  = 0.;
    params.t_end                    = time;
    params.output_filename          = "dynamics.csv";

    spacecraft::properties properties;

    properties.hysteresis_rod_volume       = rod_volume_m3;
    properties.hysteresis_params           = hymu80_params;
    properties.hysteresis_rod_orientations = {{0.0, 1.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0}};

    // Permanent Magnet: A Grade N35 NdFeB magnet
    properties.magnet_remanence   = 1.21;             // [T] for Grade N35
    properties.magnet_length      = 0.05;             // 5 cm long
    properties.magnet_diameter    = 0.01;             // 1 cm diameter
    properties.magnet_orientation = {0.0, 0.0, 1.0};  // Aligned with the body Z-axis

    auto satellite   = std::make_shared<spacecraft>(mass_g, size_m, properties);
    auto environment = std::make_shared<wmm2020_environment>(params);

    spacecraft_dynamics dynamics{satellite, environment};
    csv_state_observer  observer(params.output_filename, satellite->rods().size());
    system_state        initial;
    initial.attitude         = aos::quat::Identity();
    initial.angular_velocity = params.initial_angular_velocity;
    initial.rod_magnetizations.resize(4);
    initial.rod_magnetizations.setZero();

    using aos::abs;

    using stepper_type = runge_kutta_dopri5<system_state, double, system_state, double, vector_space_algebra>;
    auto stepper       = make_controlled<stepper_type>(epsilon, epsilon);
    integrate_adaptive(stepper, dynamics, initial, params.t_start, params.t_end, params.dt_initial, observer);

    // return aos::verify::verify_hysteresis("hymu80_hysteresis.csv", hymu80_params);
    return 0;
}
