#include "observer.hpp"

#include "aos/core/state.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>

aos::simulation::csv_state_observer::csv_state_observer(const std::string& filename, std::size_t num_rods) : _num_rods(num_rods) {
    // Use Boost.Filesystem to ensure the output directory exists
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
    *_file << "time,q_w,q_x,q_y,q_z,w_x,w_y,w_z";
    for (std::size_t i = 0; i < _num_rods; ++i) {
        *_file << ",M_" << (i + 1);
    }
    *_file << "\n";
}

void aos::simulation::csv_state_observer::operator()(const system_state& state, double time) const {
    *_file << time << ','                         //
           << state.attitude.coeffs().w() << ','  //
           << state.attitude.coeffs().x() << ','  //
           << state.attitude.coeffs().y() << ','  //
           << state.attitude.coeffs().z() << ','  //
           << state.angular_velocity.x() << ','   //
           << state.angular_velocity.y() << ','   //
           << state.angular_velocity.z();

    for (std::size_t i = 0; i < _num_rods; ++i) {
        *_file << ',' << state.rod_magnetizations(static_cast<int>(i));
    }
    *_file << '\n';
}
