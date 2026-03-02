#include "aos/cli.hpp"
#include "aos/simulation/config.hpp"
#include "aos/verify/hysteresis.hpp"

#include <exception>
#include <print>
#include <string>

auto main(int argc, char** argv) -> int {
    aos::simulation_properties properties;
    std::string                output_path;
    if (not aos::parse_cli(argc, argv, properties, output_path)) {
        return 1;
    }

    try {
        aos::verify_hysteresis(output_path, properties.satellite.hysteresis);
        return 0;
    } catch (const std::exception& ex) {
        std::println(stderr, "Error: {}", ex.what());
        return 1;
    }
}
