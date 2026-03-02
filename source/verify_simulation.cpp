#include "aos/cli.hpp"
#include "aos/components/spacecraft.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/simulation.hpp"
#include "aos/simulation/spacecraft_dynamics.hpp"
#include "aos/verify/verification_observer.hpp"

#include <exception>
#include <memory>
#include <print>
#include <string>

auto main(int argc, char** argv) -> int {
    aos::simulation_properties properties;
    std::string                output_path;
    if (not aos::parse_cli(argc, argv, properties, output_path)) {
        return 1;
    }

    try {
        auto satellite   = std::make_shared<aos::spacecraft>(properties.satellite);
        auto environment = aos::environment::create(properties.environment);
        auto dynamics    = std::make_shared<aos::spacecraft_dynamics>(satellite, environment);
        auto observer    = std::make_shared<aos::verification_observer>(output_path, satellite, environment, properties.observer);
        std::make_unique<aos::simulation>(properties, satellite, environment, dynamics, observer)->run();
        return 0;
    } catch (const std::exception& ex) {
        std::println(stderr, "Error: {}", ex.what());
        return 1;
    }
}
