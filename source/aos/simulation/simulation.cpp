#include "simulation.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/environment/orbital_mechanics.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/dynamics.hpp"
#include "aos/simulation/observer.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <print>
#include <string>

namespace aos {

void run_simulation(const std::string& output_filename, const simulation_parameters& params) {
    using aos::abs;
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_cash_karp54;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::runge_kutta_fehlberg78;
    using boost::numeric::odeint::vector_space_algebra;
    using stepper_type_f78 = runge_kutta_fehlberg78<system_state, double, system_state, double, vector_space_algebra>;
    using stepper_type_dp5 = runge_kutta_dopri5<system_state, double, system_state, double, vector_space_algebra>;
    using stepper_type_k54 = runge_kutta_cash_karp54<system_state, double, system_state, double, vector_space_algebra>;

    auto satellite   = std::make_shared<spacecraft>(params.satellite);
    auto environment = std::make_shared<environment_model>(params.environment);

    spacecraft_dynamics dynamics{satellite, environment};
    state_observer      observer(output_filename, satellite->rods().size(), params.observer);

    const auto [position, velocity] = orbital_converter::to_cartesian(params.orbit);
    const double total_duration     = params.t_end - params.t_start;
    const bool   use_checkpoints    = (params.checkpoint_interval >= 1.0);

    system_state current_state;
    current_state.position_m           = position;
    current_state.velocity_m_s         = velocity;
    current_state.attitude             = aos::quat::Identity();
    current_state.angular_velocity_m_s = params.angular_velocity;
    current_state.rod_magnetizations.resize(static_cast<std::ptrdiff_t>(params.satellite.rods.size()));
    current_state.rod_magnetizations.setZero();

    observer(current_state, params.t_start);

    auto run_integration_loop = [&](auto& stepper) {
        if (not use_checkpoints) {
            std::println("Starting simulation");
            dynamics.set_global_time_offset(0.0);
            integrate_adaptive(stepper, dynamics, current_state, params.t_start, params.t_end, params.dt_initial, observer);
        } else {
            std::println("Starting simulation with checkpoints");
            double global_time_accum = params.t_start;
            double remaining_time    = total_duration;
            double current_dt        = params.dt_initial;

            while (remaining_time > 1.0) {
                const double section_period = std::min(params.checkpoint_interval, remaining_time);

                dynamics.set_global_time_offset(global_time_accum);
                integrate_adaptive(stepper, dynamics, current_state, 0.0, section_period, current_dt);
                current_state.attitude.normalize();  // fix drift

                // in case of integrator overshot
                const auto rods = satellite->rods();
                for (std::ptrdiff_t i = 0; i < current_state.rod_magnetizations.size(); ++i) {
                    const auto& hysteresis = rods[i].hysteresis();

                    current_state.rod_magnetizations(i) = std::clamp(  //
                        current_state.rod_magnetizations(i),           //
                        -hysteresis.ms,                                //
                        hysteresis.ms                                  //
                    );
                }

                global_time_accum += section_period;
                remaining_time -= section_period;
                observer(current_state, global_time_accum);
                // observer.flush();  // comment when not needed
                std::print("Checkpoint: {} s / {} s\r", global_time_accum, params.t_end);
            }
        }
    };

    switch (params.stepper_function) {
        case 2: {
            const auto stepper = make_controlled<stepper_type_f78>(params.absolute_error, params.relative_error);
            run_integration_loop(stepper);
        } break;
        case 1: {
            const auto stepper = make_controlled<stepper_type_dp5>(params.absolute_error, params.relative_error);
            run_integration_loop(stepper);
        } break;
        case 0: {
            const auto stepper = make_controlled<stepper_type_k54>(params.absolute_error, params.relative_error);
            run_integration_loop(stepper);
        } break;
        default: {
            std::println("unknown stepper function: {}", params.stepper_function);
            break;
        }
    }

    std::print("\n");
}

}  // namespace aos
