#pragma once

#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/MagneticModel.hpp>

namespace aos {

class environment {
public:

    environment()                              = default;
    environment(const environment&)            = delete;
    environment(environment&&)                 = delete;
    environment& operator=(const environment&) = delete;
    environment& operator=(environment&&)      = delete;

    virtual ~environment();

    [[nodiscard]]
    virtual vec3 inertial_magnetic_field_at(double t_sec) const = 0;
};

/**
 * @class Environment
 * @brief Manages all external environmental models for the simulation.
 *
 * This class is responsible for calculating time-dependent environmental
 * factors, such as the spacecraft's orbital position and the Earth's
 * magnetic field at that location.
 */
class wmm2025_environment : public environment {
public:

    struct geodetic_coords {
        double lat_deg{};
        double lon_deg{};
        double alt_m{};
    };

    struct properties {
        double orbit_altitude_km{};
        double orbit_inclination_deg{};
    };

    wmm2025_environment(const wmm2025_environment&)            = delete;
    wmm2025_environment(wmm2025_environment&&)                 = delete;
    wmm2025_environment& operator=(const wmm2025_environment&) = delete;
    wmm2025_environment& operator=(wmm2025_environment&&)      = delete;

    /**
     * @brief Constructs the wmm2020_environment model.
     * @param props The properties struct containing configuration.
     */
    explicit wmm2025_environment(const properties& props);

    ~wmm2025_environment() override;

    /**
     * @brief Calculates the Earth's magnetic field in the inertial frame (ECI).
     * @param t_sec The current simulation time in seconds.
     * @return The magnetic field vector in Teslas.
     */
    [[nodiscard]]
    vec3 inertial_magnetic_field_at(double t_sec) const override;

protected:

    [[nodiscard]]
    geodetic_coords propagate_orbit_to_ecef(double t_sec) const;

    [[nodiscard]]
    vec3 get_magnetic_field_in_ned(double t_sec, const geodetic_coords& coords) const;

    [[nodiscard]]
    static vec3 rotate_ned_to_eci(const vec3& b_ned, const geodetic_coords& coords, double t_sec);

    static double normalize_longitude_deg(double lon_deg);

private:

    double                       _orbit_altitude_m;
    double                       _orbit_inclination_rad;
    double                       _orbit_radius_m;
    double                       _orbit_period_s;
    GeographicLib::MagneticModel _magnetic_model;
};

}  // namespace aos
