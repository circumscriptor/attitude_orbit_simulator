#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"

#include <memory>
#include <utility>

namespace aos {

/**
 * @class spacecraft_dynamics
 * @brief A functor that calculates the time derivative of the spacecraft's state.
 *
 * This class implements the core differential equations for the spacecraft's
 * attitude and rotational motion, designed to be used by a numerical integrator
 * like boost::odeint.
 */
class spacecraft_dynamics {
public:

    /**
     * @brief Constructs the dynamics model.
     * @param spacecraft_model A const shared pointer to the spacecraft's physical properties.
     * @param environment_model A const shared pointer to the environmental models (e.g., magnetic field).
     */
    spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const magnetic_model> environment_model)
        : _spacecraft(std::move(spacecraft_model)), _environment(std::move(environment_model)) {}

    /**
     * @brief The main operator called by the ODE solver.
     * @param current_state The current state vector (x) of the system.
     * @param state_derivative The output state derivative vector (dxdt) to be calculated.
     * @param t The current simulation time in seconds.
     */
    void operator()(const system_state& current_state, system_state& state_derivative, double t) const;

private:

    std::shared_ptr<const spacecraft>     _spacecraft;
    std::shared_ptr<const magnetic_model> _environment;
};

}  // namespace aos
