#pragma once

#include "aos/core/types.hpp"
#include "aos/simulation/config.hpp"

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/MagneticModel.hpp>

namespace aos::environment {

using core::mat3x3;
using core::vec3;
using simulation::simulation_parameters;

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
class wmm2020_environment : public environment {
public:

    wmm2020_environment(const wmm2020_environment&)            = delete;
    wmm2020_environment(wmm2020_environment&&)                 = delete;
    wmm2020_environment& operator=(const wmm2020_environment&) = delete;
    wmm2020_environment& operator=(wmm2020_environment&&)      = delete;

    /**
     * @brief Constructs the wmm2020_environment model.
     * @param params The simulation_parameters struct containing configuration.
     */
    explicit wmm2020_environment(const simulation_parameters& params)
        : m_orbit_altitude_km(params.orbit_altitude_km), m_orbit_inclination_deg(params.orbit_inclination_deg), m_magnetic_model("wmm2020") {}

    ~wmm2020_environment() override;

    /**
     * @brief Calculates the Earth's magnetic field in the inertial frame (ECI).
     * @param t_sec The current simulation time in seconds.
     * @return The magnetic field vector in Teslas.
     */
    [[nodiscard]]
    vec3 inertial_magnetic_field_at(double t_sec) const override;

private:

    double                       m_orbit_altitude_km{};
    double                       m_orbit_inclination_deg{};
    GeographicLib::MagneticModel m_magnetic_model;
};

}  // namespace aos::environment
