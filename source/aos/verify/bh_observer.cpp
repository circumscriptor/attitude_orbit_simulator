#include "bh_observer.hpp"

#include "aos/core/constants.hpp"
#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>

namespace aos {

bh_observer::bh_observer(const std::string& filename) {
    boost::filesystem::path file_path(filename);
    if (file_path.has_parent_path()) {
        boost::filesystem::create_directories(file_path.parent_path());
    }

    _file = std::make_shared<std::ofstream>(filename);
    if (not _file->is_open()) {
        throw std::runtime_error("Observer could not open output file: " + filename);
    }

    *_file << std::fixed << std::setprecision(3);

    // --- Write the CSV Header ---
    *_file << "time,H_Am,M_Am,B_T\n";
}

void bh_observer::operator()(const hysteresis_state_type& m, double t) const {
    const double h = hysteresis_loop_dynamics::h_max * sin(2.0 * std::numbers::pi * hysteresis_loop_dynamics::frequency * t);
    // Calculate B = μ₀ * (H + M)
    const double b = vacuum_permeability * (h + m);
    (*_file) << t << "," << h << "," << m << "," << b << "\n";
}

}  // namespace aos
