#pragma once

#include "aos/core/types.hpp"

namespace aos::components {

class permanent_magnet {
public:

    permanent_magnet(double remanence, double length_m, double diameter_m, const vec3& orientation);
    permanent_magnet(double br_t, double volume_m3, const vec3& orientation);

    [[nodiscard]] vec3 magnetic_moment() const { return _magnetic_moment_body; }

private:

    vec3 _magnetic_moment_body;  // [A·m^2]
};

}  // namespace aos::components
