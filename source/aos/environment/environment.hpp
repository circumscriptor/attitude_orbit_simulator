#pragma once

#include "aos/core/types.hpp"

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include <vector>

namespace aos {

struct geodetic_coords {
    double lat_deg{};
    double lon_deg{};
    double alt_m{};
};

struct environment_data {
    vec3 magnetic_field_eci_t;        // B (Tesla)
    vec3 magnetic_field_dot_eci_t_s;  // dB/dt (Tesla/s) - Material Derivative
    vec3 gravity_eci_m_s2;            // Total Gravity (m/s^2)
};

class environment_model {
public:

    environment_model(const environment_model&)                    = delete;
    environment_model(environment_model&&)                         = delete;
    auto operator=(const environment_model&) -> environment_model& = delete;
    auto operator=(environment_model&&) -> environment_model&      = delete;

    explicit environment_model(double start_year_decimal, int degree);
    virtual ~environment_model();

    /**
     * @brief Computes environmental vectors.
     * @param t_sec Global simulation time.
     * @param r_eci_m Position vector in ECI frame.
     * @param v_eci_m_s Velocity vector in ECI frame (required for magnetic gradient).
     */
    [[nodiscard]] auto calculate(double t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_data;

    [[nodiscard]] auto earth_mu() const -> double { return _gravity_model.MassConstant(); }

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
        double h_m;
        vec3   r_ecef_m;
    };

    struct field_at_point {
        vec3 b_eci;
        vec3 g_eci;
    };

    /** Compute magnetic and gravity fields at a specific spacetime point. */
    [[nodiscard]] auto compute_fields_at(double t_sec, const vec3& r_eci_m) const -> field_at_point;

    /** Compute Earth rotation and transform position ECI -> ECEF. */
    void update_ecef_transform(double t_sec, const vec3& r_eci_m) const;

    /** Convert ECEF -> Geodetic and compute ENU basis. */
    void update_geodetic_conversion() const;

private:

    double                       _start_year_decimal;
    mutable computation_cache    _cache;
    GeographicLib::Geocentric    _earth;
    GeographicLib::GravityModel  _gravity_model;
    GeographicLib::MagneticModel _magnetic_model;
};

}  // namespace aos
