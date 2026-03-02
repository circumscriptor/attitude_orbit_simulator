#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/types.hpp"

#include <span>
#include <vector>

namespace aos {

using hysteresis_rods_properties = std::vector<hysteresis_rod_properties>;

class hysteresis_rods {
public:

    hysteresis_rods(const hysteresis_rods&)                    = delete;
    hysteresis_rods(hysteresis_rods&&)                         = delete;
    auto operator=(const hysteresis_rods&) -> hysteresis_rods& = delete;
    auto operator=(hysteresis_rods&&) -> hysteresis_rods&      = delete;

    hysteresis_rods(const hysteresis_rods_properties& properties, const hysteresis_parameters& params);
    ~hysteresis_rods() = default;

    [[nodiscard]] auto rods() const -> std::span<const hysteresis_rod>;

    // compute total rod torque exerted by all rods
    [[nodiscard]] auto compute_rod_torques(const vecX& rod_magnetizations, const vec3& b_body) const -> vec3;

    // compute dM/dt for each rod, write dM/dt values into the dm_dt_out
    void compute_rod_derivatives(const vecX& rod_magnetizations, const vec3& b_body, const vec3& b_dot_body, vecX& dm_dt_out) const;

private:

    std::vector<hysteresis_rod> _rods;
};

}  // namespace aos
