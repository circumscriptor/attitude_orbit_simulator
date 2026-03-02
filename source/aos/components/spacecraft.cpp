#include "spacecraft.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft_face.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/impl/table.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <type_traits>
#include <variant>

namespace aos {

// NOLINTBEGIN(readability-magic-numbers)

void spacecraft_uniform::from_toml(const toml::table& table) {
    if (const auto* vec = table["dimensions_m"].as_array()) {
        dimensions_m <<                                      //
            vec->get(0)->value_or(default_spacecraft_size),  //
            vec->get(1)->value_or(default_spacecraft_size),  //
            vec->get(2)->value_or(default_spacecraft_size);
    }

    drag_coefficient = table["drag_coeffcient"].value_or(default_drag_coefficient);
    reflectivity     = table["reflectivity"].value_or(default_reflectivity);
}

void spacecraft_uniform::debug_print() const {
    std::cout << "--  spacecraft (uniform) shape  --"                                                                //
              << "\n  dimensions:       " << dimensions_m.x() << ' ' << dimensions_m.y() << ' ' << dimensions_m.z()  //
              << "\n  drag coefficient: " << drag_coefficient                                                        //
              << '\n';
}

void spacecraft_custom::from_toml(const toml::table& table) {
    if (const auto* vec = table["inertia"].as_array()) {
        inertia <<                                           //
            vec->get(0)->value_or(default_spacecraft_size),  //
            vec->get(1)->value_or(default_spacecraft_size),  //
            vec->get(2)->value_or(default_spacecraft_size),  //
            vec->get(3)->value_or(default_spacecraft_size),  //
            vec->get(4)->value_or(default_spacecraft_size),  //
            vec->get(5)->value_or(default_spacecraft_size),  //
            vec->get(6)->value_or(default_spacecraft_size),  //
            vec->get(7)->value_or(default_spacecraft_size),  //
            vec->get(8)->value_or(default_spacecraft_size);
    }

    if (const auto* arr = table["faces"].as_array()) {
        for (size_t i = 0; i < faces.size(); ++i) {
            if (const auto* face = arr->get(i)->as_table()) {
                faces.at(i).from_toml(*face);
            }
        }
    }
}

void spacecraft_custom::debug_print() const {
    std::cout << "--  spacecraft (custom) shape  --"                           //
              << "\n  inertia: "                                               //
              << inertia(0) << ' ' << inertia(1) << ' ' << inertia(2) << "; "  //
              << inertia(3) << ' ' << inertia(4) << ' ' << inertia(5) << "; "  //
              << inertia(6) << ' ' << inertia(7) << ' ' << inertia(8)          //
              << '\n';

    for (const auto& face : faces) {
        face.debug_print();
    }
}

void spacecraft_properties::from_toml(const toml::table& table) {
    mass_kg = table["mass_kg"].value_or(default_spacecraft_mass);

    if (const auto* uniform = table["uniform"].as_table()) {
        shape.emplace<spacecraft_uniform>().from_toml(*uniform);
    } else if (const auto* custom = table["custom"].as_table()) {
        shape.emplace<spacecraft_custom>().from_toml(*custom);
    }

    if (const auto* hyst = table["hysteresis"].as_table()) {
        hysteresis.from_toml(*hyst);
    }

    if (const auto* mag = table["magnet"].as_table()) {
        magnet.from_toml(*mag);
    }

    if (const auto* arr = table["rods"].as_array()) {
        rods.clear();
        rods.reserve(arr->size());
        for (size_t i = 0; i < arr->size(); ++i) {
            if (const auto* rod = arr->get(i)->as_table()) {
                rods.emplace_back().from_toml(*rod);
            }
        }
    }
}

void spacecraft_properties::debug_print() const {
    std::cout << "--  spacecraft properties  --"  //
              << "\n  mass: " << mass_kg          //
              << '\n';

    std::visit([](auto&& s) { s.debug_print(); }, shape);

    hysteresis.debug_print();
    magnet.debug_print();

    for (const auto& hysteresis_rod : rods) {
        hysteresis_rod.debug_print();
    }
}

// NOLINTEND(readability-magic-numbers)

spacecraft::spacecraft(const spacecraft_properties& properties) : _mass_kg(properties.mass_kg), _magnet(properties.magnet) {
    _rods.reserve(properties.rods.size());
    for (const auto& rod : properties.rods) {
        _rods.emplace_back(rod, properties.hysteresis);
    }

    std::visit(
        [&](auto&& shape) {
            using shape_type = std::decay_t<decltype(shape)>;

            if constexpr (std::is_same_v<shape_type, spacecraft_uniform>) {
                uniform(_mass_kg, shape);
            } else if constexpr (std::is_same_v<shape_type, spacecraft_custom>) {
                _inertia_tensor_kg_m2 = shape.inertia;
                _faces                = shape.faces;
            }
        },
        properties.shape);

    _inertia_tensor_kg_m2_inverse = _inertia_tensor_kg_m2.inverse();
}

auto spacecraft::mass_kg() const -> double {
    return _mass_kg;
}

auto spacecraft::inertia_tensor_kg_m2() const -> const mat3x3& {
    return _inertia_tensor_kg_m2;
}

auto spacecraft::inertia_tensor_kg_m2_inverse() const -> const mat3x3& {
    return _inertia_tensor_kg_m2_inverse;
}

auto spacecraft::faces() const -> std::span<const spacecraft_face> {
    return _faces;
}

auto spacecraft::magnet() const -> const permanent_magnet& {
    return _magnet;
}

auto spacecraft::rods() const -> std::span<const hysteresis_rod> {
    return _rods;
}

auto spacecraft::compute_gyroscopic_torque(const vec3& omega) const -> vec3 {
    return -omega.cross(_inertia_tensor_kg_m2 * omega);
}

auto spacecraft::compute_gravity_gradient_torque(const vec3& r_body, double earth_mu) const -> vec3 {
    const double r_sq  = r_body.squaredNorm();
    const double r_mag = std::sqrt(r_sq);

    // tau_gg = (3 * mu / r^5) * (r_body × (I * r_body))
    const double coef = (3.0 * earth_mu) / (r_sq * r_sq * r_mag);
    return coef * r_body.cross(_inertia_tensor_kg_m2 * r_body);
}

void spacecraft::uniform(double mass_kg, const spacecraft_uniform& shape) {
    const auto& dim_m               = shape.dimensions_m;
    const auto& drag_coefficient    = shape.drag_coefficient;
    const auto& reflectivity        = shape.reflectivity;
    const auto  face_surface_area_x = dim_m.y() * dim_m.z();
    const auto  face_surface_area_y = dim_m.x() * dim_m.z();
    const auto  face_surface_area_z = dim_m.x() * dim_m.y();

    // NOLINTBEGIN(readability-magic-numbers)
    _faces[0].surface_area_m2      = face_surface_area_x;
    _faces[0].surface_normal       = vec3(1, 0, 0);
    _faces[0].center_of_pressure_m = _faces[0].surface_normal * dim_m.x() * 0.5;
    _faces[0].drag_coefficient     = drag_coefficient;
    _faces[0].reflectivity         = reflectivity;
    _faces[1].surface_area_m2      = face_surface_area_x;
    _faces[1].surface_normal       = vec3(-1, 0, 0);
    _faces[1].center_of_pressure_m = _faces[1].surface_normal * dim_m.x() * 0.5;
    _faces[1].drag_coefficient     = drag_coefficient;
    _faces[1].reflectivity         = reflectivity;
    _faces[2].surface_area_m2      = face_surface_area_y;
    _faces[2].surface_normal       = vec3(0, 1, 0);
    _faces[2].center_of_pressure_m = _faces[2].surface_normal * dim_m.y() * 0.5;
    _faces[2].drag_coefficient     = drag_coefficient;
    _faces[2].reflectivity         = reflectivity;
    _faces[3].surface_area_m2      = face_surface_area_y;
    _faces[3].surface_normal       = vec3(0, -1, 0);
    _faces[3].center_of_pressure_m = _faces[3].surface_normal * dim_m.y() * 0.5;
    _faces[3].drag_coefficient     = drag_coefficient;
    _faces[3].reflectivity         = reflectivity;
    _faces[4].surface_area_m2      = face_surface_area_z;
    _faces[4].surface_normal       = vec3(0, 0, 1);
    _faces[4].center_of_pressure_m = _faces[4].surface_normal * dim_m.z() * 0.5;
    _faces[4].drag_coefficient     = drag_coefficient;
    _faces[4].reflectivity         = reflectivity;
    _faces[5].surface_area_m2      = face_surface_area_z;
    _faces[5].surface_normal       = vec3(0, 0, -1);
    _faces[5].center_of_pressure_m = _faces[5].surface_normal * dim_m.z() * 0.5;
    _faces[5].drag_coefficient     = drag_coefficient;
    _faces[5].reflectivity         = reflectivity;
    // NOLINTEND(readability-magic-numbers)

    _inertia_tensor_kg_m2 = compute_inertia_tensor(mass_kg, dim_m.x(), dim_m.y(), dim_m.z());
}

auto spacecraft::compute_inertia_tensor(double m_kg, double a, double b, double c) -> mat3x3 {
    const double i_x = (1. / 12.) * m_kg * (b * b + c * c);
    const double i_y = (1. / 12.) * m_kg * (a * a + c * c);
    const double i_z = (1. / 12.) * m_kg * (a * a + b * b);

    mat3x3 tensor = mat3x3::Zero();
    tensor(0, 0)  = i_x;
    tensor(1, 1)  = i_y;
    tensor(2, 2)  = i_z;
    return tensor;
}

}  // namespace aos
