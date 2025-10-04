#include "permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <cmath>
#include <numbers>

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
aos::components::permanent_magnet::permanent_magnet(double remanence, double length_m, double diameter_m, const vec3& orientation) {
    using core::vacuum_permeability;

    double br_t = remanence;

    // Calculate physical properties
    double m_am      = br_t / vacuum_permeability;  // Magnetization [A/m]
    double volume_m3 = std::numbers::pi * std::pow(diameter_m / 2.0, 2) * length_m;

    // Calculate the magnetic moment vector
    double moment_magnitude      = m_am * volume_m3;
    this->m_magnetic_moment_body = moment_magnitude * orientation.normalized();
}
