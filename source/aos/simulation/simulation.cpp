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
#include <utility>

namespace aos {

simulation::simulation(const std::string& output_filename, const simulation_properties& properties)
    : simulation(properties,
                 spacecraft::create(properties.satellite),
                 environment::create(properties.environment),
                 dynamics::create(_satellite, _environment),
                 observer::create(output_filename, properties.satellite.rods.size(), properties.observer)) {}

simulation::simulation(const simulation_properties& properties,
                       std::shared_ptr<spacecraft>  sat,
                       std::shared_ptr<environment> env,
                       std::shared_ptr<dynamics>    dyn,
                       std::shared_ptr<observer>    obs)
    : _satellite(std::move(sat)),
      _environment(std::move(env)),
      _dynamics(std::move(dyn)),
      _observer(std::move(obs)),
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
    using stepper_type_f78 = runge_kutta_fehlberg78<system_state, double, system_state, double, vector_space_algebra>;
    using stepper_type_dp5 = runge_kutta_dopri5<system_state, double, system_state, double, vector_space_algebra>;
    using stepper_type_k54 = runge_kutta_cash_karp54<system_state, double, system_state, double, vector_space_algebra>;

    _observer->write_header() << '\n';

    auto system = [this](const system_state& current_state, system_state& state_derivative, double t_sec) {
        _dynamics->step(current_state, state_derivative, t_sec);
    };

    auto observe = [this](const system_state& state, double time) { _observer->write(state, time) << '\n'; };

    auto run_integration_loop = [&](auto& stepper) {
        using boost::numeric::odeint::integrate_adaptive;

        if (_checkpoint_interval < 1.0) {
            std::println("Starting simulation");
            _dynamics->set_time_offset(0.0);
            integrate_adaptive(stepper, system, _current_state, _t_start, _t_end, _dt_initial, observe);
        } else {
            std::println("Starting simulation with checkpoints");

            _observer->write(_current_state, _t_start) << '\n';
            while (_t_now < _t_end) {
                const double section_period = std::min(_checkpoint_interval, _t_end - _t_now);

                _dynamics->set_time_offset(_t_now);
                integrate_adaptive(stepper, system, _current_state, 0.0, section_period, _dt_initial);

                fix_integration_errors();

                _t_now += section_period;
                _observer->write(_current_state, _t_now) << '\n';
                // observer.flush();  // comment when not needed
                std::print("Checkpoint: {} s / {} s\r", _t_now, _t_end);
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
