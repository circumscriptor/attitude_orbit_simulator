#include "permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <cmath>
#include <numbers>

namespace aos::components {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
permanent_magnet::permanent_magnet(double remanence, double length_m, double diameter_m, const vec3& orientation)
    : permanent_magnet(remanence, std::numbers::pi * std::pow(diameter_m / 2.0, 2) * length_m, orientation) {}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
permanent_magnet::permanent_magnet(double br_t, double volume_m3, const vec3& orientation) {
    const double m_am             = br_t / vacuum_permeability;  // Magnetization [A/m]
    const double moment_magnitude = m_am * volume_m3;
    _magnetic_moment_body         = moment_magnitude * orientation.normalized();
}

}  // namespace aos::components
