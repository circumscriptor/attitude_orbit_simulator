#include "environment.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>

#include <cmath>
#include <print>
#include <vector>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
environment_model::environment_model(double start_year_decimal, int degree)
    : _start_year_decimal(start_year_decimal),
      _earth(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f()),
      _gravity_model("egm2008", "", degree),
      _magnetic_model("wmm2025", "") {
    const auto rotation_matrix_size = 3 * 3;
    _cache.rotation_matrix_buffer.resize(rotation_matrix_size, 0.0);

    if (start_year_decimal < 1900.0 || start_year_decimal > 2100.0) {  // NOLINT(readability-magic-numbers)
        std::println(stderr, "Warning: Magnetic model year {} may be outside valid range", start_year_decimal);
    }

    std::println("Magnetic model: {} (from {}, degree {}, order {})",  //
                 _magnetic_model.MagneticModelName(),                  //
                 _magnetic_model.MagneticModelDirectory(),             //
                 _magnetic_model.Degree(),                             //
                 _magnetic_model.Order());
    std::println("Gravity model: {} (from {}, degree {}, order {})",  //
                 _gravity_model.GravityModelName(),                   //
                 _gravity_model.GravityModelDirectory(),              //
                 _gravity_model.Degree(),                             //
                 _gravity_model.Order());
}

environment_model::~environment_model() = default;

auto environment_model::calculate(double t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_data {
    // compute fields at current position
    const auto [b_current, g_current] = compute_fields_at(t_sec, r_eci_m);

    // compute fields at future position for gradient calculation
    const double t_next = t_sec + dt_gradient;
    const vec3   r_next = r_eci_m + v_eci_m_s * dt_gradient;

    const auto [b_next, g_next] = compute_fields_at(t_next, r_next);

    // compute material derivative
    const vec3 db_dt = (b_next - b_current) / dt_gradient;

    return {
        .magnetic_field_eci_t       = b_current,
        .magnetic_field_dot_eci_t_s = db_dt,
        .gravity_eci_m_s2           = g_current,
    };
}

auto environment_model::compute_fields_at(double t_sec, const vec3& r_eci_m) const -> field_at_point {
    update_ecef_transform(t_sec, r_eci_m);
    update_geodetic_conversion();

    // ENU -> ECEF -> ECI
    _cache.R_enu_to_eci = _cache.R_ecef_to_eci * _cache.R_enu_to_ecef;

    const double current_year = _start_year_decimal + (t_sec / seconds_per_year);

    double bx{};
    double by{};
    double bz{};
    double gx{};
    double gy{};
    double gz{};

    _magnetic_model(current_year, _cache.lat_deg, _cache.lon_deg, _cache.h_m, bx, by, bz);
    _gravity_model.Gravity(_cache.lat_deg, _cache.lon_deg, _cache.h_m, gx, gy, gz);

    const vec3 b_enu(bx * nanotesla_to_tesla, by * nanotesla_to_tesla, bz * nanotesla_to_tesla);
    const vec3 g_enu(gx, gy, gz);

    return {
        .b_eci = _cache.R_enu_to_eci * b_enu,
        .g_eci = _cache.R_enu_to_eci * g_enu,
    };
}

void environment_model::update_ecef_transform(double t_sec, const vec3& r_eci_m) const {
    const double theta     = earth_rotation_rate_rad_s * t_sec;
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    // ECEF -> ECI
    _cache.R_ecef_to_eci << cos_theta, -sin_theta, 0.0, sin_theta, cos_theta, 0.0, 0.0, 0.0, 1.0;

    // r_ecef = R_ecef_to_eci^T * r_eci
    _cache.r_ecef_m.x() = cos_theta * r_eci_m.x() + sin_theta * r_eci_m.y();
    _cache.r_ecef_m.y() = -sin_theta * r_eci_m.x() + cos_theta * r_eci_m.y();
    _cache.r_ecef_m.z() = r_eci_m.z();
}

void environment_model::update_geodetic_conversion() const {
    _earth.Reverse(_cache.r_ecef_m.x(), _cache.r_ecef_m.y(), _cache.r_ecef_m.z(), _cache.lat_deg, _cache.lon_deg, _cache.h_m, _cache.rotation_matrix_buffer);
    const auto& b = _cache.rotation_matrix_buffer;
    _cache.R_enu_to_ecef << b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8];  // NOLINT(readability-magic-numbers)
}

}  // namespace aos
