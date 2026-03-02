#include "hysteresis_rods.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/types.hpp"

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

namespace aos {

hysteresis_rods::hysteresis_rods(const hysteresis_rods_properties& properties, const hysteresis_parameters& params) {
    _rods.reserve(properties.size());
    for (const auto& rod : properties) {
        _rods.emplace_back(rod, params);
    }
}

auto hysteresis_rods::rods() const -> std::span<const hysteresis_rod> {
    return _rods;
}

auto hysteresis_rods::compute_rod_torques(const vecX& rod_magnetizations, const vec3& b_body) const -> vec3 {
    const auto num_rods = static_cast<std::ptrdiff_t>(_rods.size());
    assert(rod_magnetizations.size() == num_rods);

    vec3 torque_sum = vec3::Zero();
    for (std::ptrdiff_t i = 0; i < num_rods; ++i) {
        torque_sum += _rods[i].magnetic_moment(rod_magnetizations(i), b_body).cross(b_body);
    }
    return torque_sum;
}

void hysteresis_rods::compute_rod_derivatives(const vecX& rod_magnetizations, const vec3& b_body, const vec3& b_dot_body, vecX& dm_dt_out) const {
    const auto num_rods = static_cast<std::ptrdiff_t>(_rods.size());
    assert(rod_magnetizations.size() == num_rods);
    assert(dm_dt_out.size() == num_rods);

    for (std::ptrdiff_t i = 0; i < num_rods; ++i) {
        dm_dt_out(i) = _rods[i].magnetization_derivative(rod_magnetizations(i), b_body, b_dot_body);
    }
}

}  // namespace aos
