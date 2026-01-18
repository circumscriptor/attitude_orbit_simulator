#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <Eigen/src/Core/Matrix.h>

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
    spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment_model> environment_model)
        : _spacecraft(std::move(spacecraft_model)), _environment(std::move(environment_model)) {}

    /**
     * @brief The main operator called by the ODE solver.
     * @param current_state The current state vector (x) of the system.
     * @param state_derivative The output state derivative vector (dxdt) to be calculated.
     * @param t The current simulation time in seconds.
     */
    void operator()(const system_state& current_state, system_state& state_derivative, double t) const;

    void set_global_time_offset(double offset_s) { _global_time_offset = offset_s; }

protected:

    // compute total rod torque and dM/dt for each rod, returns the total torque exerted by all rods, writes dM/dt values into the dM_dt_out
    [[nodiscard]] auto compute_rod_effects(const system_state& state, const vec3& b_body, const vec3& b_dot_body, vecX& dm_dt_out) const -> vec3;

    // sums permanent magnet, hysteresis, gyroscopic, and gravity gradient torques
    [[nodiscard]] auto compute_net_torque(const vec3& omega, const vec3& b_body, const vec3& rod_torque, const vec3& r_eci, const quat& q_att) const -> vec3;

    // gravity gradient torque
    [[nodiscard]] auto compute_gravity_gradient_torque(const vec3& r_eci, const quat& q_att) const -> vec3;

    // quaternion derivative: 0.5 * q * omega
    [[nodiscard]] static auto compute_attitude_derivative(const quat& q, const vec3& omega) -> vecX;

private:

    std::shared_ptr<const spacecraft>        _spacecraft;
    std::shared_ptr<const environment_model> _environment;
    double                                   _global_time_offset{};
};

}  // namespace aos
