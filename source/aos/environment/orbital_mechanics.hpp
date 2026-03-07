#pragma once

#include "aos/core/types.hpp"  // Assuming vec3 is here

#include <utility>

namespace aos {

struct keplerian_elements {
    real semi_major_axis_m{};     // a
    real eccentricity{};          // e
    real inclination_rad{};       // i
    real raan_rad{};              // Omega (Right Ascension of Ascending Node)
    real arg_of_periapsis_rad{};  // omega (Argument of Periapsis)
    real mean_anomaly_rad{};      // M (Initial placement)

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class orbital_converter {
public:

    static constexpr real earth_mu = 3.986004418e14;  // m^3/s^2
    static constexpr real epsilon  = 1e-9;
    static constexpr int  max_iter = 100;

    /** convert Keplerian elements to ECI cartesian state vectors, return {position, velocity} in ECI */
    static auto to_cartesian(const keplerian_elements& el) -> std::pair<vec3, vec3>;
};

}  // namespace aos
