#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"

#include <Eigen/src/Core/Matrix.h>

#include <memory>

namespace aos {

/** @brief Functor that calculates the time derivative of the spacecraft's state */
class spacecraft_dynamics {
public:

    spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment_model> environment_model);

    /** @brief The main operator called by the ODE solver */
    void operator()(const system_state& current_state, system_state& state_derivative, double t_sec) const;

    void set_global_time_offset(double offset_s);

private:

    std::shared_ptr<const spacecraft>        _spacecraft;
    std::shared_ptr<const environment_model> _environment;
    double                                   _global_time_offset{};
};

}  // namespace aos
