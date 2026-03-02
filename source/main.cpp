#include "aos/cli.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/simulation.hpp"

#include <exception>
#include <memory>
#include <print>
#include <string>

auto main(int argc, char** argv) -> int {
    aos::simulation_parameters params;
    std::string                output_path;
    if (not aos::parse_cli(argc, argv, params, output_path)) {
        return 1;
    }

    try {
        std::make_unique<aos::simulation>(output_path, params)->run();
        return 0;
    } catch (const std::exception& ex) {
        std::println(stderr, "Error: {}", ex.what());
        return 1;
    }
}
