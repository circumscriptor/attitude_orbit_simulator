#include "orbital_mechanics.hpp"

#include "aos/core/types.hpp"

#include <GeographicLib/Constants.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>

namespace aos {

std::pair<vec3, vec3> orbital_converter::to_cartesian(const keplerian_elements& el) {
    // Kepler's Equation for Eccentric Anomaly (E)
    // M = E - e*sin(E)
    double e = el.mean_anomaly_rad;  // initial guess: E ~ M

    // Newton-Raphson iteration
    for (int i = 0; i < max_iter; ++i) {
        const double delta = e - (el.eccentricity * std::sin(e)) - el.mean_anomaly_rad;
        if (std::abs(delta) < epsilon) {
            break;
        }
        // Derivative dM/dE = 1 - e*cos(E)
        e = e - delta / (1.0 - el.eccentricity * std::cos(e));
    }

    // true anomaly (nu)
    // tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    const double sqrt_factor = std::sqrt((1.0 + el.eccentricity) / (1.0 - el.eccentricity));
    const double tan_e_2     = std::tan(e / 2.0);
    const double nu          = 2.0 * std::atan(sqrt_factor * tan_e_2);

    // perifocal (PQW) coordinates
    const double p        = el.semi_major_axis_m * (1.0 - el.eccentricity * el.eccentricity);  // semi-latus rectum
    const double r        = p / (1.0 + el.eccentricity * std::cos(nu));                        // radial distance
    const double h_factor = std::sqrt(earth_mu / p);                                           // flight path velocity factor

    const vec3 r_pqw(r * std::cos(nu), r * std::sin(nu), 0.0);
    const vec3 v_pqw(-h_factor * std::sin(nu), h_factor * (el.eccentricity + std::cos(nu)), 0.0);

    // perifocal (PQW) -> inertial (ECI)
    // standard orbital transformation is a 3-1-3 Euler sequence:
    // 1. rotate by argument of periapsis (omega) around Z
    // 2. rotate by inclination (i) around X
    // 3. rotate by RAAN (Omega) around Z

    const aaxis  rot_omega(el.arg_of_periapsis_rad, vec3::UnitZ());
    const aaxis  rot_inc(el.inclination_rad, vec3::UnitX());
    const aaxis  rot_raan(el.raan_rad, vec3::UnitZ());
    const mat3x3 pqw_to_eci = (rot_raan * rot_inc * rot_omega).toRotationMatrix();

    return {
        pqw_to_eci * r_pqw,
        pqw_to_eci * v_pqw,
    };
}

void keplerian_elements::debug_print() const {
    std::cout << "-- orbit properties --"                                           //
              << "\n  semi-major axis:                   " << semi_major_axis_m     //
              << "\n  eccentricity:                      " << eccentricity          //
              << "\n  inclination:                       " << inclination_rad       //
              << "\n  right ascension of ascending node: " << raan_rad              //
              << "\n  argument of perapsis:              " << arg_of_periapsis_rad  //
              << "\n  mean anomaly:                      " << mean_anomaly_rad      //
              << '\n';
}

}  // namespace aos
