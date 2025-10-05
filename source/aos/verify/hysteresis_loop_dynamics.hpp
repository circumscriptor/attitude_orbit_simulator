#pragma once

#include "aos/components/hysteresis_rod.hpp"

namespace aos::verify {

using components::hysteresis_rod;
using hysteresis_state_type = double;

// Dynamics functor for the B-H loop test
class hysteresis_loop_dynamics {
public:

    static constexpr const double h_max     = 100.0;  // [A/m] Max applied field. Must be > coercivity to see the full loop.
    static constexpr const double frequency = 1.0;    // [Hz] Frequency of the applied H-field.

    explicit hysteresis_loop_dynamics(const hysteresis_rod& rod) : _rod(&rod) {}

    void operator()(const hysteresis_state_type& m, hysteresis_state_type& dm_dt, double t) const;

private:

    const hysteresis_rod* _rod;
};

}  // namespace aos::verify
