#pragma once

#include "aos/core/types.hpp"
#include "aos/environment/nrlmsise.hpp"

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <string>
#include <vector>

namespace aos {

struct environment_effects {
    // NOLINTBEGIN(readability-identifier-naming)
    vec3   magnetic_field_eci_T;        // [Tesla] B Magnetic field
    vec3   magnetic_field_dot_eci_T_s;  // [Tesla/s] dB/dt Magnetic field time derivative
    vec3   gravity_eci_m_s2;            // [m/s^2] Total acceleration
    double atmospheric_density_kg_m3;   // [kg/m^3] Density of gasses including Anomalous Oxygen
    vec3   r_sun_eci;                   // [m] Position of Sun relative to Earth
    vec3   v_earth_rel;                 // [m/s] Earth-Spacecraft relative speed
    double shadow_factor;               // [-] 1.0 = Sun, 0.0 = Umbra, (0,1) = Penumbra
    double solar_pressure_Pa;           // [N/m^2] Radiation pressure at current distance
    // NOLINTEND(readability-identifier-naming)
};

struct environment_model_properties {
    double      start_year_decimal;
    std::string gravity_model_name;  // "egm2008"
    std::string gravity_model_path;
    std::string magnetic_model_name;  // "wmm2025"
    std::string magnetic_model_path;
    std::string weather_data_path;
    int         gravity_model_degree;
    int         gravity_model_order;
    int         magnetic_model_degree;
    int         magnetic_model_order;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

class environment {
public:

    environment(const environment&)                    = delete;
    environment(environment&&)                         = delete;
    auto operator=(const environment&) -> environment& = delete;
    auto operator=(environment&&) -> environment&      = delete;

    explicit environment(const environment_model_properties& properties);
    virtual ~environment();

    /** compute environmental effects */
    [[nodiscard]] auto compute_effects(double t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_effects;

    [[nodiscard]] auto earth_mu() const -> double;

protected:

    // avoid re-allocation
    struct computation_cache {
        std::vector<double> rotation_matrix_buffer;

        // intermediate matrices
        mat3x3 R_ecef_to_eci;  // NOLINT(readability-identifier-naming)
        mat3x3 R_enu_to_ecef;  // NOLINT(readability-identifier-naming)
        mat3x3 R_enu_to_eci;   // NOLINT(readability-identifier-naming)

        // intermediate coordinates
        double lat_deg;
        double lon_deg;
        double alt_m;
        vec3   r_ecef_m;
        vec3   r_sun_eci;  // sun position in ECI

        // other
        double current_year;
    };

    /** Compute atmospheric density at cached transform */
    [[nodiscard]] auto atmospheric_density() const -> double;

    /** Compute magnetic fields at cached transform */
    [[nodiscard]] auto magnetic_field() const -> vec3;

    /** Compute gravitational fields at cached transform */
    [[nodiscard]] auto gravitational_field() const -> vec3;

    /** Cache coordinate transformation results and matrices */
    void cache_transform(double t_sec, const vec3& r_eci_m) const;

    [[nodiscard]] static auto earth_relative_v(const vec3& v_eci_m_s, const vec3& r_eci_m) -> vec3;

    [[nodiscard]] static auto sun_position_eci(double days_since_j2000) -> vec3;

    [[nodiscard]] static auto solar_perturbation(const vec3& r_sat_eci, const vec3& r_sun_eci) -> vec3;

    [[nodiscard]] static auto earth_shadow_factor(const vec3& r_sat, const vec3& r_sun) -> double;

private:

    double                       _start_year_decimal;
    mutable computation_cache    _cache;
    GeographicLib::Geocentric    _earth;
    GeographicLib::GravityModel  _gravity_model;
    GeographicLib::MagneticModel _magnetic_model;
    nrlmsise                     _atmospheric_model;
};

}  // namespace aos
