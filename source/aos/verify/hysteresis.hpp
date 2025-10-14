#pragma once

#include "aos/components/hysteresis_rod.hpp"

#include <string>

namespace aos {

void verify_hysteresis(const std::string& output_filename, const hysteresis_rod::ja_parameters& params);

}  // namespace aos
