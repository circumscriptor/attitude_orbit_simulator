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

namespace aos {

csv_state_observer::csv_state_observer(const std::string& filename, std::size_t num_rods, const properties& props)
    : _num_rods(num_rods), _include_elements(not props.exclude_elements), _include_magnitudes(not props.exclude_magnitudes) {
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
    *_file << "time";
    if (_include_magnitudes) {
        *_file
            // << ",q"  // always norm() == 1
            << ",r"
            << ",v"
            << ",w";
    }

    if (_include_elements) {
        *_file << ",r_x,r_y,r_z"
               << ",v_x,v_y,v_z"
               << ",q_w,q_x,q_y,q_z"
               << ",w_x,w_y,w_z";
    }

    for (std::size_t i = 0; i < _num_rods; ++i) {
        *_file << ",M_" << (i + 1);
    }
    *_file << "\n";
}

void csv_state_observer::operator()(const system_state& state, double time) const {
    *_file << time;

    if (_include_magnitudes) {
        *_file  //
                // << ',' << state.attitude.norm() // ignore
            << ',' << state.position.norm()           //
            << ',' << state.velocity.norm()           //
            << ',' << state.angular_velocity.norm();  //
    }

    if (_include_elements) {
        *_file                                     //
            << ',' << state.position.x()           //
            << ',' << state.position.y()           //
            << ',' << state.position.z()           //
            << ',' << state.velocity.x()           //
            << ',' << state.velocity.y()           //
            << ',' << state.velocity.z()           //
            << ',' << state.attitude.coeffs().w()  //
            << ',' << state.attitude.coeffs().x()  //
            << ',' << state.attitude.coeffs().y()  //
            << ',' << state.attitude.coeffs().z()  //
            << ',' << state.angular_velocity.x()   //
            << ',' << state.angular_velocity.y()   //
            << ',' << state.angular_velocity.z();
    }

    for (std::size_t i = 0; i < _num_rods; ++i) {
        *_file << ',' << state.rod_magnetizations(static_cast<int>(i));
    }
    *_file << '\n';
}

}  // namespace aos
