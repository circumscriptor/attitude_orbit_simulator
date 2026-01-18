#include "observers.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/state.hpp"
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
    *_file << "time,q_w,q_x,q_y,q_z,roll_deg,pitch_deg,yaw_deg,omega_x,omega_y,omega_z\n";
}

void attitude_observer::operator()(const system_state& state, double t) const {
    // Extract Euler angles (Z-Y-X sequence)
    // Eigen's eulerAngles returns [alpha, beta, gamma]
    // 2, 1, 0 corresponds to Yaw, Pitch, Roll
    auto         euler      = state.attitude.toRotationMatrix().eulerAngles(2, 1, 0);
    const double rad_to_deg = 180.0 / std::numbers::pi;

    (*_file) << t << "," << state.attitude.w() << "," << state.attitude.x() << "," << state.attitude.y() << "," << state.attitude.z() << ","
             << euler[2] * rad_to_deg << ","  // Roll
             << euler[1] * rad_to_deg << ","  // Pitch
             << euler[0] * rad_to_deg << ","  // Yaw
             << state.angular_velocity.x() << "," << state.angular_velocity.y() << "," << state.angular_velocity.z() << "\n";
}

}  // namespace aos
