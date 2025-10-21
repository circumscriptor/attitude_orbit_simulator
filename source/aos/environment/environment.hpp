#pragma once

#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/MagneticModel.hpp>

namespace aos {

class magnetic_model {
public:

    magnetic_model()                                 = default;
    magnetic_model(const magnetic_model&)            = delete;
    magnetic_model(magnetic_model&&)                 = delete;
    magnetic_model& operator=(const magnetic_model&) = delete;
    magnetic_model& operator=(magnetic_model&&)      = delete;

    virtual ~magnetic_model();

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
class wmm2025_magnetic_model : public magnetic_model {
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

    wmm2025_magnetic_model(const wmm2025_magnetic_model&)            = delete;
    wmm2025_magnetic_model(wmm2025_magnetic_model&&)                 = delete;
    wmm2025_magnetic_model& operator=(const wmm2025_magnetic_model&) = delete;
    wmm2025_magnetic_model& operator=(wmm2025_magnetic_model&&)      = delete;

    /**
     * @brief Constructs the wmm2020_environment model.
     * @param props The properties struct containing configuration.
     */
    explicit wmm2025_magnetic_model(const properties& props);

    ~wmm2025_magnetic_model() override;

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
