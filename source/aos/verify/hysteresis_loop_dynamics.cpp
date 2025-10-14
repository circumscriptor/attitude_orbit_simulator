#include "hysteresis_loop_dynamics.hpp"

#include <cmath>
#include <numbers>

namespace aos {

void hysteresis_loop_dynamics::operator()(const hysteresis_state_type& m, hysteresis_state_type& dm_dt, double t) const {
    // Applied H-field and its time derivative
    const double h     = h_max * sin(2.0 * std::numbers::pi * frequency * t);
    const double dh_dt = h_max * 2.0 * std::numbers::pi * frequency * cos(2.0 * std::numbers::pi * frequency * t);
    dm_dt              = _rod->magnetization_derivative_from_h(m, h, dh_dt);
}

}  // namespace aos
