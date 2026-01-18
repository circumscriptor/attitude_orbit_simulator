#include "observers.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <math.h>

#include <algorithm>
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
    *_file << "time,H_Am,M_Am,B_T\n";
}

void bh_observer::operator()(const hysteresis_state_type& m, double t) const {
    const double h = hysteresis_loop_dynamics::h_max * sin(2.0 * std::numbers::pi * hysteresis_loop_dynamics::frequency * t);
    const double b = vacuum_permeability * (h + m);
    (*_file) << t << "," << h << "," << m << "," << b << "\n";
}

orbit_observer::orbit_observer(const std::string& filename) {
    _file = std::make_shared<std::ofstream>(filename);
    *_file << std::fixed << std::setprecision(3);
    *_file << "time,r_x,r_y,r_z,r_mag,v_mag\n";
}

void orbit_observer::operator()(const system_state& state, double t) const {
    const double r_mag = state.position.norm();
    const double v_mag = state.velocity.norm();

    (*_file) << t << "," << state.position.x() << "," << state.position.y() << "," << state.position.z() << "," << r_mag << "," << v_mag << "\n";
}

attitude_observer::attitude_observer(const std::string& filename) {
    _file = std::make_shared<std::ofstream>(filename);
    *_file << std::fixed << std::setprecision(3);
    *_file << "time,q_w,q_x,q_y,q_z,roll_deg,pitch_deg,yaw_deg,omega_x,omega_y,omega_z,nadir_error_deg\n";
}

void attitude_observer::operator()(const system_state& state, double t) const {
    // Extract Euler angles (Z-Y-X sequence)
    // Eigen's eulerAngles returns [alpha, beta, gamma]
    // 2, 1, 0 corresponds to Yaw, Pitch, Roll
    auto         euler       = state.attitude.toRotationMatrix().eulerAngles(2, 1, 0);
    const double rad_to_deg  = 180.0 / std::numbers::pi;
    const double nadir_error = calculate_nadir_error(state);

    (*_file) << t << "," << state.attitude.w() << "," << state.attitude.x() << "," << state.attitude.y() << "," << state.attitude.z() << ","
             << euler[2] * rad_to_deg << ","       // Roll
             << euler[1] * rad_to_deg << ","       // Pitch
             << euler[0] * rad_to_deg << ","       // Yaw
             << state.angular_velocity.x() << ","  //
             << state.angular_velocity.y() << ","  //
             << state.angular_velocity.z() << ","  //
             << nadir_error << "\n";
}

auto attitude_observer::calculate_nadir_error(const system_state& state) -> double {
    const vec3   nadir_eci     = -state.position.normalized();
    const mat3x3 R_eci_to_body = state.attitude.toRotationMatrix().transpose();  // NOLINT(readability-identifier-naming)
    const vec3   nadir_body    = R_eci_to_body * nadir_eci;
    const double cos_theta     = std::clamp(nadir_body.z(), -1.0, 1.0);
    return std::acos(cos_theta) * (180.0 / std::numbers::pi);  // NOLINT(readability-magic-numbers)
}

}  // namespace aos
