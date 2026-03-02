#include "aos/cli.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/simulation.hpp"

#include <exception>
#include <print>
#include <string>

auto main(int argc, char** argv) -> int {
    aos::simulation_parameters params;
    std::string                output_path;
    if (not aos::parse_cli(argc, argv, params, output_path)) {
        return 1;
    }

    try {
        aos::run_simulation(output_path, params);
        return 0;
    } catch (const std::exception& ex) {
        std::println(stderr, "Error: {}", ex.what());
        return 1;
    }
}
