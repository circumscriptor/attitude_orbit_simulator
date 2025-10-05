#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

namespace aos::simulation {

using components::spacecraft;
using environment::wmm2020_environment;

struct simulation_parameters {
    spacecraft::properties          spacecraft;
    wmm2020_environment::properties environment;

    vec3 initial_angular_velocity;
    // vec3   permanent_magnet_moment;
    double t_start{};
    double t_end{};
    double dt_initial{};

    static simulation_parameters get_default();
};

}  // namespace aos::simulation
