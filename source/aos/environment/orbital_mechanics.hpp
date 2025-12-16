#pragma once

#include "aos/core/types.hpp"  // Assuming vec3 is here

#include <utility>

namespace aos {

struct keplerian_elements {
    double semi_major_axis_m;     // a
    double eccentricity;          // e
    double inclination_rad;       // i
    double raan_rad;              // Omega (Right Ascension of Ascending Node)
    double arg_of_periapsis_rad;  // omega (Argument of Periapsis)
    double mean_anomaly_rad;      // M (Initial placement)

    void debug_print() const;
};

class orbital_converter {
public:

    static constexpr double earth_mu = 3.986004418e14;  // m^3/s^2
    static constexpr double epsilon  = 1e-9;
    static constexpr int    max_iter = 100;

    /** convert Keplerian elements to ECI cartesian state vectors, return {position, velocity} in ECI */
    static std::pair<vec3, vec3> to_cartesian(const keplerian_elements& el);
};

}  // namespace aos
