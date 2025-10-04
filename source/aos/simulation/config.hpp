#pragma once

#include "aos/core/types.hpp"

#include <string>

namespace aos::simulation {

using core::mat3x3;
using core::vec3;

struct simulation_parameters {
    mat3x3 inertia_tensor;
    vec3   initial_angular_velocity;
    vec3   permanent_magnet_moment;

    double t_start{};
    double t_end{};
    double dt_initial{};

    double orbit_altitude_km{};
    double orbit_inclination_deg{};

    std::string output_filename;
};

}  // namespace aos::simulation
