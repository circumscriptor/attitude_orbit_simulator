#include "aos/cli.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/simulation.hpp"

#include <exception>
#include <iostream>
#include <string>

auto main(int argc, char** argv) -> int {
    aos::simulation_parameters params;
    std::string                output_path;
    if (not aos::parse_cli(argc, argv, params, output_path)) {
        return 1;
    }

    try {
        aos::run_simulation(output_path, params);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
