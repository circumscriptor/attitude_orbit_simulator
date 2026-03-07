#include "hysteresis_loop_dynamics.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <cmath>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
hysteresis_loop_dynamics::hysteresis_loop_dynamics(const hysteresis_rod& rod, real_t h_max, real_t frequency)
    : _rod(&rod), _h_max(h_max), _frequency(frequency) {}

auto hysteresis_loop_dynamics::h_max() const noexcept -> real_t {
    return _h_max;
}

auto hysteresis_loop_dynamics::frequency() const noexcept -> real_t {
    return _frequency;
}

void hysteresis_loop_dynamics::operator()(const hysteresis_state_type& m, hysteresis_state_type& dm_dt, real_t t) const {
    // Applied H-field and its time derivative (Sinusoidal Drive)
    const real_t h     = _h_max * std::sin(2.0 * pi * _frequency * t);
    const real_t dh_dt = _h_max * 2.0 * pi * _frequency * std::cos(2.0 * pi * _frequency * t);

    // TODO: Verify that magnetization_derivative handles the 'delta' parameter (the sign of dh_dt). If dh_dt is 0 at the peak of the sine wave, ensure the
    // model doesn't produce a singularity or NaN. Also, check if the model accounts for eddy current losses at 'frequency' > 400Hz.
    dm_dt = _rod->magnetization_derivative(m, h, dh_dt);
}

}  // namespace aos
