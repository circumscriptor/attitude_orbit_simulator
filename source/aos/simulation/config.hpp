#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/environment/orbital_mechanics.hpp"
#include "aos/simulation/observer.hpp"

namespace aos {

struct simulation_properties {
    spacecraft_properties  satellite;
    keplerian_elements     orbit;
    observer_properties    observer;
    environment_properties environment;

    vec3   angular_velocity;
    double t_start{};
    double t_end{};
    double dt_initial{};
    double absolute_error{};
    double relative_error{};
    int    stepper_function{};
    double checkpoint_interval{};

    void from_toml(const toml_table& table);
    void debug_print() const;
};

}  // namespace aos
