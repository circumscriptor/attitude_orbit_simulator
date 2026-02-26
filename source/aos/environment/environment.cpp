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
      _magnetic_model("wmm2025", ""),
      _atmospheric_model("../third-party/data/SW-All.csv") {
    if (start_year_decimal < 1900.0 || start_year_decimal > 2100.0) {  // NOLINT(readability-magic-numbers)
        std::println(stderr, "Warning: Magnetic model year {} may be outside valid range", start_year_decimal);
    }

    const auto rotation_matrix_size = 3 * 3;
    _cache.rotation_matrix_buffer.resize(rotation_matrix_size, 0.0);

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
    cache_transform(t_sec, r_eci_m);

    const auto b = magnetic_field();
    const auto g = gravitational_field();

    // compute fields at future position for gradient calculation
    const double t_next = t_sec + dt_gradient;
    const vec3   r_next = r_eci_m + (v_eci_m_s * dt_gradient);
    cache_transform(t_next, r_next);

    const auto b_next = magnetic_field();

    // compute derivative
    const vec3 db_dt = (b_next - b) / dt_gradient;

    return {
        .magnetic_field_eci_T       = b,
        .magnetic_field_dot_eci_T_s = db_dt,
        .gravity_eci_m_s2           = g,
    };
}

auto environment_model::atmospheric_density() const -> double {
    return _atmospheric_model.density_at(_cache.current_year, _cache.lat_deg, _cache.lon_deg, _cache.alt_m);
}

auto environment_model::magnetic_field() const -> vec3 {
    double bx{};
    double by{};
    double bz{};

    _magnetic_model(_cache.current_year,                           // time
                    _cache.lat_deg, _cache.lon_deg, _cache.alt_m,  // geodetic
                    bx, by, bz);                                   // field (ENU)
    return _cache.R_enu_to_eci * vec3(bx * nanotesla_to_tesla, by * nanotesla_to_tesla, bz * nanotesla_to_tesla);
}

auto environment_model::gravitational_field() const -> vec3 {
    double gx_ecef{};
    double gy_ecef{};
    double gz_ecef{};

    _gravity_model.V(_cache.r_ecef_m.x(), _cache.r_ecef_m.y(), _cache.r_ecef_m.z(),  // geocentric
                     gx_ecef, gy_ecef, gz_ecef);                                     // acceleration
    return _cache.R_ecef_to_eci * vec3(gx_ecef, gy_ecef, gz_ecef);
}

void environment_model::cache_transform(double t_sec, const vec3& r_eci_m) const {
    _cache.current_year = _start_year_decimal + (t_sec / seconds_per_year);

    const double theta     = earth_rotation_rate_rad_s * t_sec;  // TODO: calculate gmst
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    // clang-format off
    _cache.R_ecef_to_eci << cos_theta, -sin_theta, 0.0,
                            sin_theta,  cos_theta, 0.0,
                                  0.0,        0.0, 1.0;
    // clang-format on

    // r_ecef = R_ecef_to_eci^T * r_eci
    _cache.r_ecef_m.x() = (cos_theta * r_eci_m.x()) + (sin_theta * r_eci_m.y());
    _cache.r_ecef_m.y() = (-sin_theta * r_eci_m.x()) + (cos_theta * r_eci_m.y());
    _cache.r_ecef_m.z() = r_eci_m.z();

    _earth.Reverse(_cache.r_ecef_m.x(), _cache.r_ecef_m.y(), _cache.r_ecef_m.z(),  // geocentric
                   _cache.lat_deg, _cache.lon_deg, _cache.alt_m,                   // geodetic
                   _cache.rotation_matrix_buffer);                                 //

    const auto& m = _cache.rotation_matrix_buffer;
    _cache.R_enu_to_ecef << m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8];  // NOLINT(readability-magic-numbers)
    _cache.R_enu_to_eci = _cache.R_ecef_to_eci * _cache.R_enu_to_ecef;
}

}  // namespace aos
