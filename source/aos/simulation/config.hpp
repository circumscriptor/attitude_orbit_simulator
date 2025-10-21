#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/observer.hpp"

namespace aos {

struct simulation_parameters {
    spacecraft::properties             spacecraft;
    wmm2025_magnetic_model::properties environment;
    csv_state_observer::properties     observer;

    vec3   angular_velocity;
    double t_start{};
    double t_end{};
    double dt_initial{};
    double absolute_error{};
    double relative_error{};
    bool   higher_order{};  // Use higher order solver

    void debug_print() const;

    static simulation_parameters get_default();
};

}  // namespace aos
