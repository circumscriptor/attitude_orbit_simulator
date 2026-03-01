#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <Eigen/src/Core/Matrix.h>

#include <memory>
#include <utility>

namespace aos {

/** @brief Functor that calculates the time derivative of the spacecraft's state */
class spacecraft_dynamics {
public:

    struct face_effects {
        vec3 torque_body;
        vec3 force_eci;
    };

    spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment_model> environment_model)
        : _spacecraft(std::move(spacecraft_model)), _environment(std::move(environment_model)) {}

    /** @brief The main operator called by the ODE solver */
    void operator()(const system_state& current_state, system_state& state_derivative, double t_sec) const;

    void set_global_time_offset(double offset_s) { _global_time_offset = offset_s; }

protected:

    // compute total rod torque and dM/dt for each rod, return the total torque exerted by all rods, write dM/dt values into the dM_dt_out
    [[nodiscard]] auto compute_rod_effects(const vecX& rod_magnetizations, const vec3& b_body, const vec3& b_dot_body, vecX& dm_dt_out) const -> vec3;

    // compute total face (drag) torque, return the total torque exerted by all faces
    [[nodiscard]] auto compute_face_effects(const environment_data& data, const quat& q_att, const vec3& v_eci, const vec3& r_eci, const vec3& omega_body) const
        -> face_effects;

    // sums permanent magnet, gyroscopic, and gravity gradient torques
    [[nodiscard]] auto compute_other_torques(const vec3& omega, const vec3& b_body, const vec3& r_eci, const mat3x3& eci_to_body) const -> vec3;

    // gravity gradient torque
    [[nodiscard]] auto compute_gravity_gradient_torque(const vec3& r_eci, const mat3x3& eci_to_body) const -> vec3;

    // quaternion derivative: 0.5 * q * omega
    [[nodiscard]] static auto compute_attitude_derivative(const quat& q_att, const vec3& omega) -> vecX;

private:

    std::shared_ptr<const spacecraft>        _spacecraft;
    std::shared_ptr<const environment_model> _environment;
    double                                   _global_time_offset{};
};

}  // namespace aos
