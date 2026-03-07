#include "permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/toml.hpp>

#include <iostream>
#include <print>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace aos {

// NOLINTBEGIN(readability-magic-numbers)
void permanent_magnet_cylindrical::from_toml(const toml_table& table) {
    length_m = table["length"].value_or(0.05);
    radius_m = table["radius"].value_or(0.005);
}

void permanent_magnet_cylindrical::debug_print() const {
    std::cout << "--  permanent magnet (cylindrical) shape  --"  //
              << "\n  length: " << length_m                      //
              << "\n  radius: " << radius_m                      //
              << '\n';
}

void permanent_magnet_rectangular::from_toml(const toml_table& table) {
    width_m  = table["width"].value_or(0.02);
    height_m = table["height"].value_or(0.02);
    length_m = table["length"].value_or(0.02);
}

void permanent_magnet_rectangular::debug_print() const {
    std::cout << "--  permanent magnet (rectangular) shape  --"  //
              << "\n  width:  " << width_m                       //
              << "\n  height: " << height_m                      //
              << "\n  length: " << length_m                      //
              << '\n';
}

void permanent_magnet_properties::from_toml(const toml_table& table) {
    remanence_t           = table["remanence"].value_or(1.21);
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
                set_cylindrical(shape);
            } else if constexpr (std::is_same_v<shape_type, permanent_magnet_rectangular>) {
                set_rectangular(shape);
            }
        },
        properties.shape);

    if (_volume_m3 <= 0.0) {
        throw std::invalid_argument("Magnet volume must be positive.");
    }

    _magnetic_moment_body = compute_magnetic_moment(_remanence_t, _volume_m3, _demagnetization_factor, _relative_permeability, _orientation_body);

    std::println("Permanent magnet remanence: {}", _remanence_t);
    std::println("Permanent magnet demagnitization factor: {}", _demagnetization_factor);
    std::println("Permanent magnet magnetic moment (body frame): {}, {}, {}", _magnetic_moment_body.x(), _magnetic_moment_body.y(), _magnetic_moment_body.z());
}

auto permanent_magnet::magnetic_moment() const -> vec3 {
    return _magnetic_moment_body;
}

auto permanent_magnet::compute_torque(const vec3& b_field_body) const -> vec3 {
    return _magnetic_moment_body.cross(b_field_body);
}

auto permanent_magnet::compute_potential_energy(const vec3& b_field_body) const -> real_t {
    return -_magnetic_moment_body.dot(b_field_body);
}

auto permanent_magnet::compute_magnetic_moment_at_temperature(real_t temp_celsius, real_t temp_coeff, real_t temp_ref) -> vec3 {
    const real_t temp_factor        = 1.0 + (temp_coeff * (temp_celsius - temp_ref));
    const real_t adjusted_remanence = _remanence_t * temp_factor;
    return compute_magnetic_moment(adjusted_remanence, _volume_m3, _demagnetization_factor, _relative_permeability, _orientation_body);
}

void permanent_magnet::set_rectangular(const permanent_magnet_rectangular& shape) {
    _volume_m3 = shape.width_m * shape.height_m * shape.length_m;

    const real_t area_z     = shape.width_m * shape.height_m;
    const real_t area_x     = shape.height_m * shape.length_m;
    const real_t area_y     = shape.width_m * shape.length_m;
    _demagnetization_factor = area_z / (area_z + area_x + area_y);
}

void permanent_magnet::set_cylindrical(const permanent_magnet_cylindrical& shape) {
    const real_t r = shape.radius_m;
    const real_t h = shape.length_m;
    _volume_m3     = pi * (r * r) * h;

    const real_t area_face  = pi * (r * r);
    const real_t area_side  = 2.0 * pi * r * h;
    _demagnetization_factor = area_face / (area_face + area_side);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto permanent_magnet::compute_magnetic_moment(real_t      remanence,
                                               real_t      volume,
                                               real_t      demagnitization_factor,
                                               real_t      relative_permeability,
                                               const vec3& orientation) -> vec3 {
    // M_eff = (Br / mu0) / (1 + N * (mu_r - 1))
    const real_t magnetization_raw = remanence / vacuum_permeability;
    const real_t chi_m             = relative_permeability - 1.0;
    const real_t magnetization_eff = magnetization_raw / (1.0 + demagnitization_factor * chi_m);
    const real_t moment_mag        = magnetization_eff * volume;
    return moment_mag * orientation;
}

}  // namespace aos
