#pragma once

#include "aos/simulation/config.hpp"

#include <string>

namespace aos::simulation {

void run_simulation(const std::string& output_filename, const simulation_parameters& params);

}  // namespace aos::simulation
