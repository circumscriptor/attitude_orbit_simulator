#include "environment_impl.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/environment/nrlmsise.hpp"

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include <algorithm>
#include <cmath>
#include <print>
#include <vector>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
environment_impl::environment_impl(const environment_properties& properties)
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

environment_impl::~environment_impl() = default;

auto environment_impl::compute_effects(real t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_effects {
    // compute fields at current position
    cache_transform(t_sec, r_eci_m);

    const auto b       = magnetic_field();
    const auto d       = atmospheric_density();
    const auto g_earth = gravitational_field();
    const auto g_sun   = solar_perturbation(r_eci_m, _cache.r_sun_eci);
    const auto v_rel   = earth_relative_v(v_eci_m_s, r_eci_m);
    const auto g_total = g_earth + g_sun;

    const vec3& r_sun    = _cache.r_sun_eci;
    const real  d_sun_sq = r_sun.squaredNorm();
    const real  pressure = solar_pressure_1au * (au_to_m_2 / d_sun_sq);
    const real  shadow   = earth_shadow_factor(r_eci_m, r_sun);
    const real  earth_mu = _gravity_model.MassConstant();

    // compute fields at future position for gradient calculation
    const real t_next = t_sec + dt_gradient;
    const vec3 r_next = r_eci_m + (v_eci_m_s * dt_gradient);
    cache_transform(t_next, r_next);

    const auto b_next = magnetic_field();

    // compute derivative
    const vec3 db_dt = (b_next - b) / dt_gradient;

    return {
        .magnetic_field_eci_T       = b,
        .magnetic_field_dot_eci_T_s = db_dt,
        .gravity_eci_m_s2           = g_total,
        .atmospheric_density_kg_m3  = d,
        .r_sun_eci                  = r_sun,
        .v_earth_rel                = v_rel,
        .shadow_factor              = shadow,
        .solar_pressure_Pa          = pressure,
        .earth_mu                   = earth_mu,
    };
}

void environment_impl::cache_transform(real t_sec, const vec3& r_eci_m) const {
    _cache.current_year = _start_year_decimal + (t_sec / seconds_per_year);

    const real days_j2000 = (_cache.current_year - 2000.0) * 365.25;
    _cache.r_sun_eci      = sun_position_eci(days_j2000);

    const real theta     = earth_rotation_rate_rad_s * t_sec;  // TODO: calculate gmst
    const real cos_theta = std::cos(theta);
    const real sin_theta = std::sin(theta);

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

auto environment_impl::atmospheric_density() const -> real {
    return _atmospheric_model.density_at(_cache.current_year, _cache.lat_deg, _cache.lon_deg, _cache.alt_m);
}

auto environment_impl::magnetic_field() const -> vec3 {
    double bx{};
    double by{};
    double bz{};

    _magnetic_model(_cache.current_year,                           // time
                    _cache.lat_deg, _cache.lon_deg, _cache.alt_m,  // geodetic
                    bx, by, bz);                                   // field (ENU)
    return _cache.R_enu_to_eci * vec3(bx * nanotesla_to_tesla, by * nanotesla_to_tesla, bz * nanotesla_to_tesla);
}

auto environment_impl::gravitational_field() const -> vec3 {
    double gx_ecef{};
    double gy_ecef{};
    double gz_ecef{};

    _gravity_model.V(_cache.r_ecef_m.x(), _cache.r_ecef_m.y(), _cache.r_ecef_m.z(),  // geocentric
                     gx_ecef, gy_ecef, gz_ecef);                                     // acceleration
    return _cache.R_ecef_to_eci * vec3(gx_ecef, gy_ecef, gz_ecef);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto environment_impl::earth_relative_v(const vec3& v_eci_m_s, const vec3& r_eci_m) -> vec3 {
    const vec3 omega_earth_eci(0.0, 0.0, earth_rotation_rate_rad_s);  // Earth rotation vector in ECI
    const vec3 v_atm_eci = omega_earth_eci.cross(r_eci_m);            // Linear velocity of the atmosphere at position r_eci
    return v_eci_m_s - v_atm_eci;                                     // Velocity of satellite relative to the air
}

auto environment_impl::sun_position_eci(real days_since_j2000) -> vec3 {
    const real g       = (357.528 + 0.9856003 * days_since_j2000) * deg_to_rad;
    const real l       = (280.460 + 0.9856474 * days_since_j2000) * deg_to_rad;
    const real sin_g   = std::sin(g);
    const real cos_g   = std::cos(g);
    const real lambda  = l + ((1.915 * sin_g + 0.020 * (2.0 * sin_g * cos_g)) * deg_to_rad);
    const real epsilon = (23.439 - 0.0000004 * days_since_j2000) * deg_to_rad;
    const real r_au    = 1.00014 - (0.01671 * cos_g) - (0.00014 * (2.0 * cos_g * cos_g - 1.0));
    const real r_m     = r_au * au_to_m;
    const real sin_l   = std::sin(lambda);
    const real cos_l   = std::cos(lambda);
    const real sin_e   = std::sin(epsilon);
    const real cos_e   = std::cos(epsilon);
    return {
        r_m * cos_l,
        r_m * cos_e * sin_l,
        r_m * sin_e * sin_l,
    };
}

auto environment_impl::solar_perturbation(const vec3& r_sat_eci, const vec3& r_sun_eci) -> vec3 {
    const vec3 r_rel    = r_sun_eci - r_sat_eci;
    const real d_rel_sq = r_rel.squaredNorm();
    const real d_sun_sq = r_sun_eci.squaredNorm();
    // (d^2)^(1.5) = d^3
    const real d_rel_cubed = d_rel_sq * std::sqrt(d_rel_sq);
    const real d_sun_cubed = d_sun_sq * std::sqrt(d_sun_sq);
    // a = mu * (r_rel/d_rel^3 - r_sun/d_sun^3)
    return sun_mu_m3_s2 * ((r_rel / d_rel_cubed) - (r_sun_eci / d_sun_cubed));
}

auto environment_impl::earth_shadow_factor(const vec3& r_sat, const vec3& r_sun) -> real {
    const real d_sat     = r_sat.norm();
    const vec3 unit_sat  = r_sat / d_sat;
    const vec3 unit_sun  = r_sun.normalized();
    const real cos_theta = unit_sat.dot(unit_sun);
    if (cos_theta > 0) {
        return 1.0;  // Sun is on the same side as satellite
    }

    const real d_sat_sun = (r_sun - r_sat).norm();

    // Apparent angular radii (half-angles)
    const real alpha = std::asin(sun_radius_m / d_sat_sun);
    const real beta  = std::asin(earth_radius_m / d_sat);
    const real gamma = std::acos(std::clamp(-cos_theta, -1.0, 1.0));

    if (gamma > alpha + beta) {
        return 1.0;  // Full sunlight
    }
    if (gamma < beta - alpha) {
        return 0.0;  // Umbra (total eclipse)
    }
    // Penumbra: Linear transition based on angular overlap
    return std::clamp((gamma + alpha - beta) / (2.0 * alpha), 0.0, 1.0);
}

}  // namespace aos
