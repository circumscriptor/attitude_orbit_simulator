#pragma once

#include "aos/simulation/config.hpp"

#include <string>

namespace aos {

// Tests if gravity gradient pulls the satellite back to Nadir
void verify_attitude(const std::string& output_filename, const simulation_parameters& params);

}  // namespace aos
