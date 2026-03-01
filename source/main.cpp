#include "aos/simulation/config.hpp"
#include "aos/simulation/simulation.hpp"
#include "aos/verify/attitude.hpp"
#include "aos/verify/hysteresis.hpp"
#include "aos/verify/orbit.hpp"

#include <toml++/impl/parse_error.hpp>
#include <toml++/impl/parser.hpp>

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

auto main(int argc, char** argv) -> int {
    std::string config_path       = "config.toml";
    std::string output_path       = "output.csv";
    bool        verify_attitude   = false;
    bool        verify_hysteresis = false;
    bool        verify_orbit      = false;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];  // NOLINT

        if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: simulator [config.toml] [options]\n"
                      << "Options:\n"
                      << "  -o, --output <file>      Output file (default: output.csv)\n"
                      << "  --verify-attitude        Run attitude verification\n"
                      << "  --verify-hysteresis      Run hysteresis verification\n"
                      << "  --verify-orbit           Run orbit verification\n";
            return 0;
        }

        if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                output_path = argv[++i];  // NOLINT
            }
        } else if (arg == "--verify-attitude") {
            verify_attitude = true;
        } else if (arg == "--verify-hysteresis") {
            verify_hysteresis = true;
        } else if (arg == "--verify-orbit") {
            verify_orbit = true;
        } else if (arg[0] != '-') {
            config_path = arg;
        }
    }

    aos::simulation_parameters params;
    if (not std::filesystem::exists(config_path)) {
        std::cerr << "Error: File '" << config_path << "' not found.\n";
        return 1;
    }

    try {
        auto table = toml::parse_file(config_path);
        params.from_toml(table);
    } catch (const toml::parse_error& err) {
        std::cerr << "TOML Error: " << err << "\n";
        return 1;
    }

    try {
        params.debug_print();

        if (verify_hysteresis) {
            if (params.satellite.rods.empty()) {
                return 1;
            }
            aos::verify_hysteresis(output_path, params.satellite.hysteresis);
        } else if (verify_attitude) {
            aos::verify_attitude(output_path, params);
        } else if (verify_orbit) {
            aos::verify_orbit(output_path, params);
        } else {
            aos::run_simulation(output_path, params);
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
