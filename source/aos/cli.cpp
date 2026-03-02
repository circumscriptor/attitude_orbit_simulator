#include "cli.hpp"

#include "aos/simulation/config.hpp"

#include <toml++/impl/parse_error.hpp>
#include <toml++/impl/parser.hpp>

#include <cstddef>
#include <filesystem>
#include <print>
#include <ranges>
#include <span>
#include <string>
#include <string_view>

namespace aos {

auto parse_cli(int argc, char** argv, simulation_properties& properties, std::string& output_path) -> bool {
    std::string config_path = "config.toml";
    output_path             = "output.csv";

    auto args = std::span(argv, argc)                                                     //
                | std::views::transform([](char* arg) { return std::string_view(arg); })  //
                | std::views::drop(1);

    size_t skip_count = 0;
    for (auto [i, arg] : args | std::views::enumerate) {
        if (skip_count > 0) {
            --skip_count;
            continue;
        }

        if (arg == "-h" || arg == "--help") {
            std::println(
                "Usage: simulator [config.toml] [options]\n"
                "Options:\n"
                "  -o, --output <file>      Output file (default: output.csv)");
            return false;
        }

        if (arg == "-o" || arg == "--output") {
            if (i + 1 < std::ranges::ssize(args)) {
                output_path = args[i + 1];
                skip_count  = 1;
            }
        } else if (not arg.starts_with('-')) {
            config_path = arg;
        } else {
            std::println(stderr, "Error: Unknown option '{}'", arg);
            return false;
        }
    }

    if (not std::filesystem::exists(config_path)) {
        std::println(stderr, "Error: File '{}' not found.", config_path);
        return false;
    }

    try {
        auto table = toml::parse_file(config_path);
        properties.from_toml(table);
        properties.debug_print();
        return true;
    } catch (const toml::parse_error& err) {
        std::println(stderr, "TOML Error: {}", err.description());
        return false;
    }
}

}  // namespace aos
