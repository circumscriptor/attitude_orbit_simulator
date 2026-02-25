#pragma once

#include "aos/core/types.hpp"

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

extern "C" {
#include <nrlmsise-00.h>
}

#include <vector>

namespace aos {

struct geodetic_coords {
    double lat_deg{};
    double lon_deg{};
    double alt_m{};
};

struct environment_data {
    // NOLINTBEGIN(readability-identifier-naming)
    vec3 magnetic_field_eci_T;        // B (Tesla)
    vec3 magnetic_field_dot_eci_T_s;  // dB/dt (Tesla/s) - Material Derivative
    vec3 gravity_eci_m_s2;            // Total Gravity (m/s^2)
    // NOLINTEND(readability-identifier-naming)
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
        double              current_year;
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

        // atmospheric model
        nrlmsise_input  am_input;
        nrlmsise_flags  am_flags;
        nrlmsise_output am_output;
    };

    /** Compute magnetic fields at cached transform */
    [[nodiscard]] auto magnetic_field() const -> vec3;

    /** Compute gravitational fields at cached transform */
    [[nodiscard]] auto gravitational_field() const -> vec3;

    /** Compute atmospheric density at a specific spacetime point */
    [[nodiscard]] auto density_at(double t_sec, double lat_deg, double lon_deg, double alt_m) const -> double;

    /** Cache coordinate transformation results and matrices */
    void cache_transform(double t_sec, const vec3& r_eci_m) const;

private:

    double                       _start_year_decimal;
    mutable computation_cache    _cache;
    GeographicLib::Geocentric    _earth;
    GeographicLib::GravityModel  _gravity_model;
    GeographicLib::MagneticModel _magnetic_model;
};

}  // namespace aos
