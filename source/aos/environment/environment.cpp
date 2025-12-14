#include "environment.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>

#include <cmath>
#include <print>
#include <utility>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
environment_model::environment_model(double start_year_decimal, int degree)
    : _start_year_decimal(start_year_decimal),
      _earth(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f()),
      _gravity_model("egm2008", "", degree),
      _magnetic_model("wmm2025", "") {
    const auto rotation_matrix_size = 3 * 3;
    _cache.rotation_matrix_buffer.resize(rotation_matrix_size, 0.0);

    std::println("Magnetic model: {} (from {}, degree {}, order {})", _magnetic_model.MagneticModelName(), _magnetic_model.MagneticModelDirectory(),
                 _magnetic_model.Degree(), _magnetic_model.Order());
    std::println("Gravity model: {} (from {}, degree {}, order {})", _gravity_model.GravityModelName(), _gravity_model.GravityModelDirectory(),
                 _gravity_model.Degree(), _gravity_model.Order());
}

environment_model::~environment_model() = default;

environment_data environment_model::calculate(double t_sec, const vec3& r_eci_m) const {
    // update earth rotation and ecef position
    step_eci_to_ecef(t_sec, r_eci_m);

    // update geodetic coords and local basis matrix
    step_geodetic_conversion();

    // get raw physics data in local tangent plane (ENU)
    const auto [mag_enu, grav_enu] = step_compute_physics_enu(t_sec);

    // transform to inertial frame
    return step_rotate_to_eci(mag_enu, grav_enu);
}

void environment_model::step_eci_to_ecef(double t_sec, const vec3& r_eci_m) const {
    const double theta = earth_rotation_rate_rad_s * t_sec;
    const double c     = std::cos(theta);
    const double s     = std::sin(theta);

    // ECI <-> ECEF (z-axis)
    _cache.R_ecef_to_eci << c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0;

    // r_ecef = R_ecef_to_eci.transpose() * r_eci
    _cache.r_ecef_m.x() = c * r_eci_m.x() + s * r_eci_m.y();
    _cache.r_ecef_m.y() = -s * r_eci_m.x() + c * r_eci_m.y();
    _cache.r_ecef_m.z() = r_eci_m.z();
}

void environment_model::step_geodetic_conversion() const {
    _earth.Reverse(_cache.r_ecef_m.x(), _cache.r_ecef_m.y(), _cache.r_ecef_m.z(), _cache.lat_deg, _cache.lon_deg, _cache.h_m, _cache.rotation_matrix_buffer);
    const auto& buf = _cache.rotation_matrix_buffer;
    _cache.R_enu_to_ecef << buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7], buf[8];  // NOLINT(readability-magic-numbers)
}

std::pair<vec3, vec3> environment_model::step_compute_physics_enu(double t_sec) const {
    const double current_year = _start_year_decimal + (t_sec / seconds_per_year);

    // Gravity
    double gx{};
    double gy{};
    double gz{};
    _gravity_model.Disturbance(_cache.lat_deg, _cache.lon_deg, _cache.h_m, gx, gy, gz);
    const vec3 grav_enu(gx, gy, gz);

    // Magnetic field
    double bx{};
    double by{};
    double bz{};
    _magnetic_model(current_year, _cache.lat_deg, _cache.lon_deg, _cache.h_m, bx, by, bz);
    const vec3 mag_enu(bx * nanotesla_to_tesla, by * nanotesla_to_tesla, bz * nanotesla_to_tesla);

    return {
        mag_enu,
        grav_enu,
    };
}

environment_data environment_model::step_rotate_to_eci(const vec3& mag_enu, const vec3& grav_enu) const {
    // ENU -> ECEF -> ECI
    _cache.R_enu_to_eci = _cache.R_ecef_to_eci * _cache.R_enu_to_ecef;
    return {
        .magnetic_field_eci_t         = _cache.R_enu_to_eci * mag_enu,
        .gravity_disturbance_eci_m_s2 = _cache.R_enu_to_eci * grav_enu,
    };
}

}  // namespace aos
