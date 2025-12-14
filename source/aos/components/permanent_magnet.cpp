#include "permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <limits>
#include <numbers>
#include <stdexcept>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
permanent_magnet::permanent_magnet(double remanence_t, double volume_m3, const vec3& orientation) : _remanence(remanence_t), _volume(volume_m3) {
    if (volume_m3 <= 0.0) {
        throw std::invalid_argument("Magnet volume must be positive.");
    }

    const double norm = orientation.norm();
    if (norm <= std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument("Magnet orientation vector cannot be zero.");
    }

    _orientation_body = orientation.normalized();

    // M = Br / mu0  [A/m]
    // m = M * V     [A*m^2]
    const double magnetization = _remanence / vacuum_permeability;
    const double moment_mag    = magnetization * _volume;
    _magnetic_moment_body      = moment_mag * _orientation_body;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
permanent_magnet permanent_magnet::cylindrical(double remanence_t, double length_m, double diameter_m, const vec3& orientation) {
    double radius = diameter_m / 2.0;
    double volume = std::numbers::pi * (radius * radius) * length_m;
    return {
        remanence_t,
        volume,
        orientation,
    };
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
permanent_magnet permanent_magnet::rectangular(double remanence_t, double width_m, double height_m, double length_m, const vec3& orientation) {
    double volume = width_m * height_m * length_m;
    return {
        remanence_t,
        volume,
        orientation,
    };
}

void permanent_magnet::update_temperature(double temp_celsius, double temp_coeff, double ref_temp) {
    const double temp_factor        = 1.0 + (temp_coeff * (temp_celsius - ref_temp));
    const double adjusted_remanence = _remanence * temp_factor;
    const double magnetization      = adjusted_remanence / vacuum_permeability;
    _magnetic_moment_body           = (magnetization * _volume) * _orientation_body;
}

}  // namespace aos
