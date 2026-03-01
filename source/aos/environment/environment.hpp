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

struct environment_data {
    // NOLINTBEGIN(readability-identifier-naming)
    vec3   magnetic_field_eci_T;        // [Tesla] - B
    vec3   magnetic_field_dot_eci_T_s;  // [Tesla/s] - (dB/dt) Material Derivative
    vec3   gravity_eci_m_s2;            // [m/s^2] - Total Acceleration
    double atmospheric_density_kg_m3;   // [kg/m^3] - Density of Gasses including Anomalous Oxygen
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

class environment_model {
public:

    environment_model(const environment_model&)                    = delete;
    environment_model(environment_model&&)                         = delete;
    auto operator=(const environment_model&) -> environment_model& = delete;
    auto operator=(environment_model&&) -> environment_model&      = delete;

    explicit environment_model(const environment_model_properties& properties);
    virtual ~environment_model();

    /**
     * @brief Compute environmental effects
     * @param t_sec Global simulation time
     * @param r_eci_m Position vector in ECI frame
     * @param v_eci_m_s Velocity vector in ECI frame (required for magnetic gradient)
     */
    [[nodiscard]] auto calculate(double t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_data;

    [[nodiscard]] auto earth_mu() const -> double;

    /** @brief Get spacecraft velocity relative to earth rotation */
    [[nodiscard]] static auto earth_relative_v(const vec3& v_eci_m_s) -> vec3;

    [[nodiscard]] static auto sun_position_eci(double days_since_j2000) -> vec3;

    [[nodiscard]] static auto solar_perturbation(const vec3& r_sat_eci, const vec3& r_sun_eci) -> vec3;

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

private:

    double                       _start_year_decimal;
    mutable computation_cache    _cache;
    GeographicLib::Geocentric    _earth;
    GeographicLib::GravityModel  _gravity_model;
    GeographicLib::MagneticModel _magnetic_model;
    nrlmsise                     _atmospheric_model;
};

}  // namespace aos
