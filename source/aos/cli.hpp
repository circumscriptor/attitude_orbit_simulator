#pragma once

#include "aos/simulation/config.hpp"

#include <string>

namespace aos {

auto parse_cli(int argc, char** argv, simulation_parameters& params, std::string& output_path) -> bool;

}  // namespace aos
