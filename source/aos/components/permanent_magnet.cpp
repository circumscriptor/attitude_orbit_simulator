#include "permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/impl/table.hpp>

#include <iostream>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace aos {

// NOLINTBEGIN(readability-magic-numbers)
void permanent_magnet_cylindrical::from_toml(const toml::table& table) {
    length_m = table["length_m"].value_or(0.05);
    radius_m = table["radius_m"].value_or(0.005);
}

void permanent_magnet_cylindrical::debug_print() const {
    std::cout << "--  permanent magnet (cylindrical) shape  --"  //
              << "\n  length: " << length_m                      //
              << "\n  radius: " << radius_m                      //
              << '\n';
}

void permanent_magnet_rectangular::from_toml(const toml::table& table) {
    width_m  = table["width_m"].value_or(0.02);
    height_m = table["height_m"].value_or(0.02);
    length_m = table["length_m"].value_or(0.02);
}

void permanent_magnet_rectangular::debug_print() const {
    std::cout << "--  permanent magnet (rectangular) shape  --"  //
              << "\n  width:  " << width_m                       //
              << "\n  height: " << height_m                      //
              << "\n  length: " << length_m                      //
              << '\n';
}

void permanent_magnet_properties::from_toml(const toml::table& table) {
    remanence_t           = table["remanence_t"].value_or(1.21);
    relative_permeability = table["relative_permeability"].value_or(1.0);

    if (const auto* vec = table["orientation"].as_array()) {
        orientation <<                   //
            vec->get(0)->value_or(0.0),  //
            vec->get(1)->value_or(0.0),  //
            vec->get(2)->value_or(1.0);
    }

    if (const auto* cyl = table["cylindrical"].as_table()) {
        shape.emplace<permanent_magnet_cylindrical>().from_toml(*cyl);
    } else if (const auto* rect = table["rectangular"].as_table()) {
        shape.emplace<permanent_magnet_rectangular>().from_toml(*rect);
    }
}

void permanent_magnet_properties::debug_print() const {
    std::cout << "--  permanent magnet properties  --"                                                                 //
              << "\n  remanence:             " << remanence_t                                                          //
              << "\n  relative permeability: " << relative_permeability                                                //
              << "\n  orientation:           " << orientation.x() << ' ' << orientation.y() << ' ' << orientation.z()  //
              << '\n';

    std::visit([](auto&& s) { s.debug_print(); }, shape);
}
// NOLINTEND(readability-magic-numbers)

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
