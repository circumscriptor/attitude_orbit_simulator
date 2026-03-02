#pragma once

#include "aos/simulation/config.hpp"

#include <string>

namespace aos {

auto parse_cli(int argc, char** argv, simulation_properties& properties, std::string& output_path) -> bool;

}  // namespace aos
