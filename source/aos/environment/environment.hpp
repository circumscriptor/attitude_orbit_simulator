#pragma once

#include "aos/core/types.hpp"

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include <utility>
#include <vector>

namespace aos {

struct geodetic_coords {
    double lat_deg{};
    double lon_deg{};
    double alt_m{};
};

struct environment_data {
    vec3 magnetic_field_eci_t;          // Tesla
    vec3 gravity_disturbance_eci_m_s2;  // m/s^2
};

class environment_model {
public:

    environment_model(const environment_model&)            = delete;
    environment_model(environment_model&&)                 = delete;
    environment_model& operator=(const environment_model&) = delete;
    environment_model& operator=(environment_model&&)      = delete;

    explicit environment_model(double start_year_decimal);
    virtual ~environment_model();

    [[nodiscard]]
    environment_data calculate(double t_sec, const vec3& r_eci_m) const;

    [[nodiscard]]
    double earth_mu() const {
        return _gravity_model.MassConstant();
    }

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

    /** Compute Earth rotation and transform position ECI -> ECEF. */
    void step_eci_to_ecef(double t_sec, const vec3& r_eci_m) const;

    /**  Convert ECEF -> Geodetic and compute ENU basis. */
    void step_geodetic_conversion() const;

    /** Query GeographicLib models for ENU vectors. Returns { Magnetic_ENU, Gravity_ENU } */
    [[nodiscard]] std::pair<vec3, vec3> step_compute_physics_enu(double t_sec) const;

    /**  Combine rotations and transform ENU vectors to ECI. */
    [[nodiscard]] environment_data step_rotate_to_eci(const vec3& mag_enu, const vec3& grav_enu) const;

private:

    double                       _start_year_decimal;
    mutable computation_cache    _cache;
    GeographicLib::Geocentric    _earth;
    GeographicLib::GravityModel  _gravity_model;
    GeographicLib::MagneticModel _magnetic_model;
};

}  // namespace aos
