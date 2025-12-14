#pragma once

#include "aos/core/types.hpp"

namespace aos {

class permanent_magnet {
public:

    static constexpr double default_temp_coeff = -0.0002;
    static constexpr double default_ref_temp   = 20.0;

    /**
     * @brief Construct from generic volume and Remanence.
     */
    permanent_magnet(double remanence_t, double volume_m3, const vec3& orientation);

    /**
     * @brief Static factory for a Cylindrical Magnet.
     */
    static permanent_magnet cylindrical(double remanence_t, double length_m, double diameter_m, const vec3& orientation);

    /**
     * @brief Static factory for a Rectangular Bar Magnet.
     */
    static permanent_magnet rectangular(double remanence_t, double width_m, double height_m, double length_m, const vec3& orientation);

    [[nodiscard]] vec3   magnetic_moment() const { return _magnetic_moment_body; }
    [[nodiscard]] double remanence() const { return _remanence; }
    [[nodiscard]] double volume() const { return _volume; }

    /**
     * @brief Re-calculates moment using a temperature coefficient (e.g., -0.11% per C for Alnico).
     */
    void update_temperature(double temp_celsius, double temp_coeff = default_temp_coeff, double ref_temp = default_ref_temp);

private:

    double _remanence;             // [Tesla]
    double _volume;                // [m^3]
    vec3   _orientation_body;      // Unit vector
    vec3   _magnetic_moment_body;  // [A·m^2]
};

}  // namespace aos
