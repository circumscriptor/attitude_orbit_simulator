#include "permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace aos {

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
permanent_magnet::permanent_magnet(const permanent_magnet_properties& properties)
    : _remanence_t(properties.remanence_t), _relative_permeability(properties.relative_permeability), _orientation_body(properties.orientation.normalized()) {
    std::visit(
        [&](auto&& shape) {
            using shape_type = std::decay_t<decltype(shape)>;

            if constexpr (std::is_same_v<shape_type, permanent_magnet_cylindrical>) {
                const double r = shape.radius_m;
                const double h = shape.length_m;
                _volume_m3     = std::numbers::pi * (r * r) * h;

                const double area_face  = std::numbers::pi * (r * r);
                const double area_side  = 2.0 * std::numbers::pi * r * h;
                _demagnetization_factor = area_face / (area_face + area_side);
            } else if constexpr (std::is_same_v<shape_type, permanent_magnet_rectangular>) {
                _volume_m3 = shape.width_m * shape.height_m * shape.length_m;

                const double area_z     = shape.width_m * shape.height_m;
                const double area_x     = shape.height_m * shape.length_m;
                const double area_y     = shape.width_m * shape.length_m;
                _demagnetization_factor = area_z / (area_z + area_x + area_y);
            }
        },
        properties.shape);

    if (_volume_m3 <= 0.0) {
        throw std::invalid_argument("Magnet volume must be positive.");
    }
}

auto permanent_magnet::magnetic_moment() const -> vec3 {
    return compute_magnetic_moment(_remanence_t, _volume_m3, _demagnetization_factor, _relative_permeability, _orientation_body);
}

auto permanent_magnet::magnetic_moment_at_temperature(double temp_celsius, double temp_coeff, double ref_temp) -> vec3 {
    const double temp_factor        = 1.0 + (temp_coeff * (temp_celsius - ref_temp));
    const double adjusted_remanence = _remanence_t * temp_factor;
    return compute_magnetic_moment(adjusted_remanence, _volume_m3, _demagnetization_factor, _relative_permeability, _orientation_body);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto permanent_magnet::compute_magnetic_moment(double remanence, double vol, double demag, double mu_r, const vec3& orientation) -> vec3 {
    // M_eff = (Br / mu0) / (1 + N * (mu_r - 1))
    const double magnetization_raw = remanence / vacuum_permeability;
    const double chi_m             = mu_r - 1.0;
    const double magnetization_eff = magnetization_raw / (1.0 + demag * chi_m);
    const double moment_mag        = magnetization_eff * vol;
    return moment_mag * orientation;
}

}  // namespace aos
