#pragma once

#include "aos/simulation/config.hpp"

#include <string>

namespace aos {

// Integrates the orbit for several periods to visualize J2 effects
void verify_orbit(const std::string& output_filename, const simulation_parameters& params);

}  // namespace aos
