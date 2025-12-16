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
    *_file << "time";

    if (_include_magnitudes) {
        *_file << ",r"
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
        *_file << ',' << state.position.norm()           // r
               << ',' << state.velocity.norm()           // v
               << ',' << state.angular_velocity.norm();  // w
    }

    if (_include_elements) {
        *_file << ',' << state.position.x()           // r_x
               << ',' << state.position.y()           // r_y
               << ',' << state.position.z()           // r_z
               << ',' << state.velocity.x()           // v_x
               << ',' << state.velocity.y()           // v_y
               << ',' << state.velocity.z()           // v_z
               << ',' << state.attitude.coeffs().w()  // q_w
               << ',' << state.attitude.coeffs().x()  // q_x
               << ',' << state.attitude.coeffs().y()  // q_y
               << ',' << state.attitude.coeffs().z()  // q_z
               << ',' << state.angular_velocity.x()   // w_x
               << ',' << state.angular_velocity.y()   // w_y
               << ',' << state.angular_velocity.z();  // w_z
    }

    for (std::size_t i = 0; i < _num_rods; ++i) {
        *_file << ',' << state.rod_magnetizations(static_cast<int>(i));
    }
    *_file << '\n';
}

}  // namespace aos
