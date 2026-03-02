#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/environment/orbital_mechanics.hpp"
#include "aos/simulation/observer.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

namespace aos {

struct simulation_properties {
    spacecraft_properties        satellite;
    keplerian_elements           orbit;
    state_observer_properties    observer;
    environment_model_properties environment;

    vec3   angular_velocity;
    double t_start{};
    double t_end{};
    double dt_initial{};
    double absolute_error{};
    double relative_error{};
    int    stepper_function{};
    double checkpoint_interval{};

    void from_toml(const toml::table& table);

    void debug_print() const;
};

}  // namespace aos
