#include "environment.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <toml++/impl/table.hpp>

#include <cmath>
#include <iostream>
#include <print>
#include <string>
#include <vector>

namespace aos {

void environment_model_properties::from_toml(const toml::table& table) {
    // NOLINTBEGIN(readability-magic-numbers)

    start_year_decimal    = table["start_year_decimal"].value_or(2026.0);
    gravity_model_name    = table["gravity_model_name"].value_or<std::string>("egm2008");
    gravity_model_path    = table["gravity_model_path"].value_or<std::string>("");
    magnetic_model_name   = table["magnetic_model_name"].value_or<std::string>("wmm2025");
    magnetic_model_path   = table["magnetic_model_path"].value_or<std::string>("");
    weather_data_path     = table["weather_data_path"].value_or<std::string>(WEATHER_DATA_PATH);
    gravity_model_degree  = table["gravity_model_degree"].value_or(-1);
    gravity_model_order   = table["gravity_model_order"].value_or(-1);
    magnetic_model_degree = table["magnetic_model_degree"].value_or(-1);
    magnetic_model_order  = table["magnetic_model_order"].value_or(-1);

    // NOLINTEND(readability-magic-numbers)
}

void environment_model_properties::debug_print() const {
    std::cout << "--  spacecraft (custom) shape  --"                     //
              << "\n  start year:            " << start_year_decimal     //
              << "\n  gravity model name:    " << gravity_model_name     //
              << "\n  gravity model path:    " << gravity_model_path     //
              << "\n  magnetic model name:   " << magnetic_model_name    //
              << "\n  magnetic model path:   " << magnetic_model_path    //
              << "\n  weather data path:     " << weather_data_path      //
              << "\n  gravity model degree:  " << gravity_model_degree   //
              << "\n  gravity model degree:  " << gravity_model_order    //
              << "\n  magnetic model degree: " << magnetic_model_degree  //
              << "\n  magnetic model degree: " << magnetic_model_order   //
              << '\n';
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
environment_model::environment_model(const environment_model_properties& properties)
    : _start_year_decimal(properties.start_year_decimal),
      _earth(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f()),
      _gravity_model(properties.gravity_model_name, properties.gravity_model_path, properties.gravity_model_degree, properties.gravity_model_order),
      _magnetic_model(properties.magnetic_model_name,
                      properties.magnetic_model_path,
                      GeographicLib::Geocentric::WGS84(),
                      properties.magnetic_model_degree,
                      properties.magnetic_model_order),
      _atmospheric_model(properties.weather_data_path) {
    if (_start_year_decimal < 1900.0 || _start_year_decimal > 2100.0) {  // NOLINT(readability-magic-numbers)
        std::println(stderr, "Warning: Magnetic model year {} may be outside valid range", _start_year_decimal);
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
    const auto d = atmospheric_density();
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
        .atmospheric_density_kg_m3  = d,
    };
}

auto environment_model::earth_mu() const -> double {
    return _gravity_model.MassConstant();
}

auto environment_model::earth_relative_v(const vec3& v_eci_m_s) -> vec3 {
    return v_eci_m_s - vec3(0., 0., earth_rotation_rate_rad_s);
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
