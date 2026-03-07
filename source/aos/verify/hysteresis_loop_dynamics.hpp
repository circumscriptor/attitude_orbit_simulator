#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/types.hpp"

namespace aos {

using hysteresis_state_type = real;

// Dynamics functor for the B-H loop test
class hysteresis_loop_dynamics {
public:

    static constexpr const real default_h_max     = 1000.0;  // [A/m] Max applied field. Must be > coercivity to see the full loop.
    static constexpr const real default_frequency = 1.0;     // [Hz] Frequency of the applied H-field.

    explicit hysteresis_loop_dynamics(const hysteresis_rod& rod, real h_max = default_h_max, real frequency = default_frequency);

    [[nodiscard]] auto h_max() const noexcept -> real;
    [[nodiscard]] auto frequency() const noexcept -> real;

    void operator()(const hysteresis_state_type& m, hysteresis_state_type& dm_dt, real t) const;

private:

    const hysteresis_rod* _rod;
    real                  _h_max;
    real                  _frequency;
};

}  // namespace aos
