#include "hysteresis_observer.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"
#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
hysteresis_observer::hysteresis_observer(const std::string& filename, real h_max, real frequency) : _h_max(h_max), _frequency(frequency) {
    std::filesystem::path file_path(filename);
    if (file_path.has_parent_path()) {
        std::filesystem::create_directories(file_path.parent_path());
    }

    _file = std::make_shared<std::ofstream>(filename);
    if (not _file->is_open()) {
        throw std::runtime_error("Observer could not open output file: " + filename);
    }

    *_file << std::fixed << std::setprecision(3);
    *_file << "time,H_Am,M_Am,B_T\n";
}

void hysteresis_observer::operator()(const hysteresis_state_type& m, real t) const {
    const real h = _h_max * std::sin(2.0 * pi * _frequency * t);
    const real b = vacuum_permeability * (h + m);
    (*_file) << t << "," << h << "," << m << "," << b << "\n";
}

}  // namespace aos
