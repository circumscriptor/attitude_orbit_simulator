#include "environment.hpp"

#include <GeographicLib/Constants.hpp>

#include <cmath>
#include <numbers>

aos::environment::environment::~environment() = default;

aos::environment::wmm2020_environment::~wmm2020_environment() = default;

aos::environment::vec3 aos::environment::wmm2020_environment::inertial_magnetic_field_at(double t_sec) const {
    // --- 1. Orbit Propagation ---
    // This section implements a simple circular orbit propagator to find the
    // spacecraft's position in the Earth-Centered, Earth-Fixed (ECEF) frame.

    const double earth_radius_m            = GeographicLib::Constants::WGS84_a();
    const double earth_mu                  = 3.986004418e14;  // Earth's gravitational parameter [m^3/s^2]
    const double earth_rotation_rate_rad_s = 7.2921150e-5;    // [rad/s]

    const double alt_m          = m_orbit_altitude_km * 1000.0;
    const double orbit_radius_m = earth_radius_m + alt_m;

    // Calculate the orbital period for a circular orbit
    const double orbit_period_s = 2.0 * std::numbers::pi * std::sqrt(std::pow(orbit_radius_m, 3) / earth_mu);

    // The angle of the satellite in its orbital plane (mean anomaly)
    const double orbit_angle_rad = 2.0 * std::numbers::pi * t_sec / orbit_period_s;

    const double inclination_rad = m_orbit_inclination_deg * std::numbers::pi / 180.0;

    // --- 2. Coordinate Conversion (Orbital to ECEF) ---
    // Calculate position in the ECEF frame, which GeographicLib expects.
    // This includes the effect of Earth's rotation.

    // Latitude calculation for an inclined circular orbit
    const double lat_rad = std::asin(std::sin(inclination_rad) * std::sin(orbit_angle_rad));
    const double lat_deg = lat_rad * 180.0 / std::numbers::pi;

    // Longitude calculation must account for both orbital motion and Earth's rotation
    const double lon_orbital_rad  = std::atan2(std::cos(inclination_rad) * std::sin(orbit_angle_rad), std::cos(orbit_angle_rad));
    const double lon_rotation_rad = earth_rotation_rate_rad_s * t_sec;

    // NOLINTBEGIN(readability-magic-numbers)
    // Combine and wrap longitude to the [-180, 180] degree range
    double lon_deg = (lon_orbital_rad - lon_rotation_rad) * 180.0 / std::numbers::pi;
    lon_deg        = std::fmod(lon_deg + 180.0, 360.0) - 180.0;
    if (lon_deg < -180.0) {
        lon_deg += 360.0;
    }
    // NOLINTEND(readability-magic-numbers)

    // --- 3. Magnetic Field Calculation ---
    // GeographicLib requires the year as a decimal.
    // We'll use a fixed epoch for simplicity. A more advanced model could
    // get this from the simulation start time.
    const double year_decimal = 2025.0 + (t_sec / (365.25 * 24 * 3600));

    double bx_nt{};
    double by_nt{};
    double bz_nt{};  // Output from the model is in nanoTeslas (nT)
    m_magnetic_model(year_decimal, lat_deg, lon_deg, alt_m, bx_nt, by_nt, bz_nt);

    // Convert from ECEF field (North, East, Down) to the standard ECI-like
    // cartesian frame (X, Y, Z) and convert units from nT to T.
    // NOTE: For this simplified orbit, we approximate the ECEF frame of the magnetic
    // field model as being aligned with the ECI frame at t=0. A high-fidelity
    // simulation would require a full ECI to ECEF rotation matrix.
    const double nt_to_t = 1e-9;  // nanotesla to tesla
    return {bx_nt * nt_to_t, by_nt * nt_to_t, bz_nt * nt_to_t};
}
