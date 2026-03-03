#pragma once

#include "aos/core/types.hpp"

#include <memory>
#include <string>

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
    double earth_mu;                    // [?] Earth mass
    // NOLINTEND(readability-identifier-naming)
};

struct environment_properties {
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

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class environment {
public:

    environment()                                      = default;
    environment(const environment&)                    = delete;
    environment(environment&&)                         = delete;
    auto operator=(const environment&) -> environment& = delete;
    auto operator=(environment&&) -> environment&      = delete;

    virtual ~environment();

    // compute environmental effects
    [[nodiscard]] virtual auto compute_effects(double t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_effects = 0;

    static auto create(const environment_properties& properties) -> std::shared_ptr<environment>;
};

}  // namespace aos
