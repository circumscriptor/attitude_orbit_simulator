#pragma once
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include "../core/types.hpp"

#include "../simulation/config.hpp"

namespace aos::environment {

using core::vec3;
using simulation::SimulationParameters;

// Holds all external environmental models
class environment {
public:

    explicit environment(const SimulationParameters& params) : params_(params), m_magnetic_model("wmm2020") {}

    // Calculate the state of the environment at a given time
    [[nodiscard]]
    vec3 inertial_magnetic_field_at(double t_sec) const {
        // ... (Orbit propagation logic from the old dynamics class is moved here) ...
        double alt_m = params_.orbit_altitude_km * 1000.0;
        // ... calculation of lat, lon ...

        double bx, by, bz;  // nT
        m_magnetic_model(t_sec / 31536000.0 + 2025.0, lat, lon, alt_m, bx, by, bz);
        return vec3(bx * 1e-9, by * 1e-9, bz * 1e-9);  // Return in Tesla
    }

private:

    const SimulationParameters&  params_;
    GeographicLib::MagneticModel m_magnetic_model;
};

}  // namespace aos::environment
