#include "observer_impl.hpp"

#include "aos/core/state.hpp"
#include "aos/simulation/observer.hpp"

#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <ios>
#include <ostream>
#include <stdexcept>
#include <string>

namespace aos {

observer_impl::observer_impl(const std::string& filename, std::size_t num_rods, const observer_properties& properties)
    : _num_rods(num_rods), _include_elements(not properties.exclude_elements), _include_magnitudes(not properties.exclude_magnitudes) {
    std::filesystem::path file_path(filename);
    if (file_path.has_parent_path()) {
        std::filesystem::create_directories(file_path.parent_path());
    }

    _file.open(filename);
    if (not _file.is_open()) {
        throw std::runtime_error("Observer could not open output file: " + filename);
    }

    _file << std::fixed << std::setprecision(properties.precission);
}

observer_impl::~observer_impl() = default;

auto observer_impl::write_header() -> std::ostream& {
    _file << "time";

    if (_include_magnitudes) {
        _file << ",r"
              << ",v"
              << ",w";
    }

    if (_include_elements) {
        _file << ",r_x,r_y,r_z"
              << ",v_x,v_y,v_z"
              << ",q_w,q_x,q_y,q_z"
              << ",w_x,w_y,w_z";
    }

    for (std::size_t i = 0; i < _num_rods; ++i) {
        _file << ",M_" << (i + 1);
    }

    return _file;
}

auto observer_impl::write(const system_state& state, double time) -> std::ostream& {
    _file << time;

    if (_include_magnitudes) {
        _file << ',' << state.position_m.norm()             // r
              << ',' << state.velocity_m_s.norm()           // v
              << ',' << state.angular_velocity_m_s.norm();  // w
    }

    if (_include_elements) {
        _file << ',' << state.position_m.x()             // r_x
              << ',' << state.position_m.y()             // r_y
              << ',' << state.position_m.z()             // r_z
              << ',' << state.velocity_m_s.x()           // v_x
              << ',' << state.velocity_m_s.y()           // v_y
              << ',' << state.velocity_m_s.z()           // v_z
              << ',' << state.attitude.coeffs().w()      // q_w
              << ',' << state.attitude.coeffs().x()      // q_x
              << ',' << state.attitude.coeffs().y()      // q_y
              << ',' << state.attitude.coeffs().z()      // q_z
              << ',' << state.angular_velocity_m_s.x()   // w_x
              << ',' << state.angular_velocity_m_s.y()   // w_y
              << ',' << state.angular_velocity_m_s.z();  // w_z
    }

    for (std::size_t i = 0; i < _num_rods; ++i) {
        _file << ',' << state.rod_magnetizations(static_cast<int>(i));
    }

    return _file;
}

}  // namespace aos
