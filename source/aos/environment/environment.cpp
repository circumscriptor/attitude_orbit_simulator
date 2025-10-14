#include "environment.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>

#include <cmath>
#include <numbers>

namespace aos {

environment::~environment() = default;

wmm2025_environment::wmm2025_environment(const properties& props)
    : _orbit_altitude_m(props.orbit_altitude_km * km_to_m),
      _orbit_inclination_rad(props.orbit_inclination_deg * deg_to_rad),
      _orbit_radius_m(GeographicLib::Constants::WGS84_a() + _orbit_altitude_m),
      _orbit_period_s(2.0 * std::numbers::pi * std::sqrt(std::pow(_orbit_radius_m, 3) / earth_mu_m3_s2)),
      _magnetic_model("wmm2025") {}

wmm2025_environment::~wmm2025_environment() = default;

vec3 wmm2025_environment::inertial_magnetic_field_at(double t_sec) const {
    // 1. Calculate spacecraft position in ECEF frame (lat, lon, alt)
    const geodetic_coords coords_ecef = propagate_orbit_to_ecef(t_sec);

    // 2. Get the magnetic field in the local North-East-Down (NED) frame.
    // The NED frame is a local tangent frame where:
    // +X (North) is tangent to the meridian, pointing towards the North Pole.
    // +Y (East) is tangent to the parallel, pointing East.
    // +Z (Down) is normal to the ellipsoid, pointing towards the Earth's center.
    const vec3 b_ned = get_magnetic_field_in_ned(t_sec, coords_ecef);

    // 3. Rotate the NED vector to the ECI frame to get the inertial field.
    return rotate_ned_to_eci(b_ned, coords_ecef, t_sec);
}

wmm2025_environment::geodetic_coords wmm2025_environment::propagate_orbit_to_ecef(double t_sec) const {
    const double orbit_angle_rad = two_pi * t_sec / _orbit_period_s;

    // Note: For near-polar orbits, the asin() function can have numerical sensitivity
    // at the poles (latitudes near +/- 90 deg). For this model, it is assumed to be acceptable.
    const double lat_rad          = std::asin(std::sin(_orbit_inclination_rad) * std::sin(orbit_angle_rad));
    const double lon_orbital_rad  = std::atan2(std::cos(_orbit_inclination_rad) * std::sin(orbit_angle_rad), std::cos(orbit_angle_rad));
    const double lon_rotation_rad = earth_rotation_rate_rad_s * t_sec;

    const double lon_deg = (lon_orbital_rad - lon_rotation_rad) * rad_to_deg;

    return {.lat_deg = lat_rad * rad_to_deg, .lon_deg = normalize_longitude_deg(lon_deg), .alt_m = _orbit_altitude_m};
}

vec3 wmm2025_environment::get_magnetic_field_in_ned(double t_sec, const geodetic_coords& coords) const {
    const double year_decimal = simulation_start_year + (t_sec / seconds_per_year);

    // GeographicLib MagneticModel operator() returns components as:
    // Bx = easterly component (East)
    // By = northerly component (North)
    // Bz = vertical (up) component (Up)
    double bx_nt{};  // East component
    double by_nt{};  // North component
    double bz_nt{};  // Up component

    // Call the magnetic model function
    _magnetic_model(year_decimal, coords.lat_deg, coords.lon_deg, coords.alt_m, bx_nt, by_nt, bz_nt);

    // Construct the NED vector: (North, East, Down).
    // Note that the Down component is the negation of the Up component.
    return {
        by_nt * nt_to_t,   // North
        bx_nt * nt_to_t,   // East
        -bz_nt * nt_to_t,  // Down
    };
}

vec3 wmm2025_environment::rotate_ned_to_eci(const vec3& b_ned, const geodetic_coords& coords, double t_sec) {
    // Step 1: Rotate from North-East-Down (NED) to Earth-Centered, Earth-Fixed (ECEF)
    const double lat     = coords.lat_deg * deg_to_rad;
    const double lon     = coords.lon_deg * deg_to_rad;
    const double cos_lat = std::cos(lat);
    const double sin_lat = std::sin(lat);
    const double cos_lon = std::cos(lon);
    const double sin_lon = std::sin(lon);

    // The columns of the NED-to-ECEF matrix are the NED basis vectors expressed in ECEF coordinates.
    vec3 north_vec = {-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat};
    vec3 east_vec  = {-sin_lon, cos_lon, 0.0};
    vec3 down_vec  = {-cos_lat * cos_lon, -cos_lat * sin_lon, -sin_lat};

    // Construct the matrix from the column vectors. This is the robust Eigen way.
    mat3x3 r_ned_to_ecef;
    r_ned_to_ecef.col(0) = north_vec;
    r_ned_to_ecef.col(1) = east_vec;
    r_ned_to_ecef.col(2) = down_vec;

    const vec3 b_ecef = r_ned_to_ecef * b_ned;

    // Step 2: Rotate from ECEF to Earth-Centered Inertial (ECI)
    const double earth_rotation_angle = earth_rotation_rate_rad_s * t_sec;
    const double cos_rot              = std::cos(earth_rotation_angle);
    const double sin_rot              = std::sin(earth_rotation_angle);

    mat3x3 r_ecef_to_eci;
    r_ecef_to_eci << cos_rot, -sin_rot, 0.0, sin_rot, cos_rot, 0.0, 0.0, 0.0, 1.0;

    return r_ecef_to_eci * b_ecef;
}

double wmm2025_environment::normalize_longitude_deg(double lon_deg) {
    // NOLINTBEGIN(readability-magic-numbers)
    double lon = std::fmod(lon_deg + 180.0, 360.0);
    if (lon < 0) {
        lon += 360.0;
    }
    return lon - 180.0;
    // NOLINTEND(readability-magic-numbers)
}

}  // namespace aos
