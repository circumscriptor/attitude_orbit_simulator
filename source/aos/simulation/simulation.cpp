#include "simulation.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/constants.hpp"
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
#include <stdexcept>
#include <string>
#include <utility>

namespace aos {

simulation::simulation(const std::string& output_filename, const simulation_properties& properties)
    : simulation(output_filename, properties, spacecraft::create(properties.satellite), environment::create(properties.environment)) {}

simulation::simulation(const std::string&                  output_filename,
                       const simulation_properties&        properties,
                       const std::shared_ptr<spacecraft>&  satellite,
                       const std::shared_ptr<environment>& environment)
    : simulation(properties,
                 satellite,
                 environment,
                 dynamics::create(satellite, environment),
                 observer::create(output_filename, properties.satellite.rods.size(), properties.observer)) {}

simulation::simulation(const simulation_properties& properties,
                       std::shared_ptr<spacecraft>  satellite,
                       std::shared_ptr<environment> environment,
                       std::shared_ptr<dynamics>    dynamics,
                       std::shared_ptr<observer>    observer)
    : _satellite(std::move(satellite)),
      _environment(std::move(environment)),
      _dynamics(std::move(dynamics)),
      _observer(std::move(observer)),
      _t_start(properties.t_start),
      _t_end(properties.t_end),
      _t_now(_t_start),
      _dt_initial(properties.dt_initial),
      _checkpoint_interval(properties.checkpoint_interval),
      _absolute_error(properties.absolute_error),
      _relative_error(properties.relative_error),
      _stepper_function(properties.stepper_function) {
    const auto [position, velocity]     = orbital_converter::to_cartesian(properties.orbit);
    _current_state.position_m           = position;
    _current_state.velocity_m_s         = velocity;
    _current_state.attitude             = aos::quat::Identity();
    _current_state.angular_velocity_m_s = properties.angular_velocity;
    _current_state.rod_magnetizations.resize(static_cast<std::ptrdiff_t>(properties.satellite.rods.size()));
    _current_state.rod_magnetizations.setZero();
}

void simulation::run() {
    using aos::abs;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_cash_karp54;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::runge_kutta_fehlberg78;
    using boost::numeric::odeint::vector_space_algebra;
    using stepper_type_f78 = runge_kutta_fehlberg78<system_state, real, system_state, real, vector_space_algebra>;
    using stepper_type_dp5 = runge_kutta_dopri5<system_state, real, system_state, real, vector_space_algebra>;
    using stepper_type_k54 = runge_kutta_cash_karp54<system_state, real, system_state, real, vector_space_algebra>;

    _observer->write_header() << '\n';

    auto system = [this](const system_state& current_state, system_state& state_derivative, real t_sec) {
        _dynamics->step(current_state, state_derivative, t_sec);
    };

    auto observe = [this](const system_state& state, real time) {
        _observer->write(state, time) << '\n';

        if (state.altitude_m() <= reentry_altitude_m) {
            throw std::runtime_error("Deorbited");
        }

        if (state.has_nan()) {
            throw std::runtime_error("Numerical Instability");
        }
    };

    auto run_integration_loop = [&](auto& stepper) {
        using boost::numeric::odeint::integrate_adaptive;

        if (_checkpoint_interval < 1.0) {
            std::println("Starting simulation");
            _dynamics->set_time_offset(0.0);

            try {
                integrate_adaptive(stepper, system, _current_state, _t_start, _t_end, _dt_initial, observe);
            } catch (const std::runtime_error& e) {
                std::println("\n[Terminated] Simulation stopped early: {}", e.what());
            }
        } else {
            std::println("Starting simulation with checkpoints");

            _observer->write(_current_state, _t_start) << '\n';
            while (_t_now < _t_end) {
                const auto section_period = std::min(_checkpoint_interval, _t_end - _t_now);

                _dynamics->set_time_offset(_t_now);
                integrate_adaptive(stepper, system, _current_state, 0.0, section_period, _dt_initial);

                fix_integration_errors();

                _t_now += section_period;
                _observer->write(_current_state, _t_now) << '\n';
                // observer.flush();  // comment when not needed
                std::print("Checkpoint: {} s / {} s\r", _t_now, _t_end);

                if (const auto altitude_m = _current_state.altitude_m(); altitude_m <= reentry_altitude_m) {
                    const auto altitude_km = altitude_m * meter_to_kilometer;
                    std::println("\n[Terminated] Satellite deorbited at t = {:.1f} s. Altitude: {:.2f} km", _t_now + section_period, altitude_km);
                    break;
                }

                if (_current_state.has_nan()) {
                    std::println("\n[Terminated] Numerical instability (NaN detected) at t = {:.1f} s.", _t_now + section_period);
                    break;
                }
            }
        }
    };

    switch (_stepper_function) {
        case 2: {
            const auto stepper = make_controlled<stepper_type_f78>(_absolute_error, _relative_error);
            run_integration_loop(stepper);
        } break;
        case 1: {
            const auto stepper = make_controlled<stepper_type_dp5>(_absolute_error, _relative_error);
            run_integration_loop(stepper);
        } break;
        case 0: {
            const auto stepper = make_controlled<stepper_type_k54>(_absolute_error, _relative_error);
            run_integration_loop(stepper);
        } break;
        default: {
            std::println("Error: Unknown stepper function: {}", _stepper_function);
            break;
        }
    }

    std::print("\n");
}

void simulation::fix_integration_errors() {
    _current_state.attitude.normalize();  // fix drift

    // in case of integrator overshot
    const auto rods = _satellite->hystresis().rods();
    for (std::ptrdiff_t i = 0; i < _current_state.rod_magnetizations.size(); ++i) {
        const auto& hysteresis = rods[i].hysteresis();

        _current_state.rod_magnetizations(i) = std::clamp(  //
            _current_state.rod_magnetizations(i),           //
            -hysteresis.ms,                                 //
            hysteresis.ms                                   //
        );
    }
}

}  // namespace aos
