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
    real_t t_start{};
    real_t t_end{};
    real_t dt_initial{};
    real_t absolute_error{};
    real_t relative_error{};
    int    stepper_function{};
    real_t checkpoint_interval{};

    void from_toml(const toml_table& table);
    void debug_print() const;
};

}  // namespace aos
