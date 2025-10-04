#pragma once

#include "aos/core/types.hpp"

namespace aos::components {

using core::vec3;

class permanent_magnet {
public:

    permanent_magnet(double remanence, double length_m, double diameter_m, const core::vec3& orientation);

    [[nodiscard]] vec3 magnetic_moment() const { return m_magnetic_moment_body; }

private:

    core::vec3 m_magnetic_moment_body;  // [A·m^2]
};

}  // namespace aos::components
