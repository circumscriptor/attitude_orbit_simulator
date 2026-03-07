#pragma once

#include "aos/core/types.hpp"  // Assuming vec3 is here

#include <utility>

namespace aos {

struct keplerian_elements {
    real_t semi_major_axis_m{};     // a
    real_t eccentricity{};          // e
    real_t inclination_rad{};       // i
    real_t raan_rad{};              // Omega (Right Ascension of Ascending Node)
    real_t arg_of_periapsis_rad{};  // omega (Argument of Periapsis)
    real_t mean_anomaly_rad{};      // M (Initial placement)

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class orbital_converter {
public:

    static constexpr real_t earth_mu = 3.986004418e14;  // m^3/s^2
    static constexpr real_t epsilon  = 1e-9;
    static constexpr int    max_iter = 100;

    /** convert Keplerian elements to ECI cartesian state vectors, return {position, velocity} in ECI */
    static auto to_cartesian(const keplerian_elements& el) -> std::pair<vec3, vec3>;
};

}  // namespace aos
