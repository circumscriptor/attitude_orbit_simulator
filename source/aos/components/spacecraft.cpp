#include "spacecraft.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/types.hpp"

#include <iostream>
#include <span>
#include <type_traits>
#include <variant>

namespace aos {

void spacecraft_properties::debug_print() const {
    // std::cout << "--  spacecraft properties  --"                                                  //
    //           << "\n  mass:               " << mass_g                                             //
    //           << "\n  dimensions:         " << dim_m.x() << ' ' << dim_m.y() << ' ' << dim_m.z()  //
    //           << '\n';

    //   magnet.debug_print();

    for (const auto& hysteresis_rod : rods) {
        hysteresis_rod.debug_print();
    }
    std::cout << '\n';
}

spacecraft::spacecraft(const spacecraft_properties& properties) : _magnet(properties.magnet) {
    _rods.reserve(properties.rods.size());
    for (const auto& rod : properties.rods) {
        _rods.emplace_back(rod);
    }

    std::visit(
        [&](auto&& shape) {
            using shape_type = std::decay_t<decltype(shape)>;

            if constexpr (std::is_same_v<shape_type, spacecraft_uniform>) {
                uniform(shape.mass_g, shape.dimensions_m, shape.drag_coefficient);
            } else if constexpr (std::is_same_v<shape_type, spacecraft_custom>) {
                _inertia_tensor = shape.inertia;
                _faces          = shape.faces;
            }
        },
        properties.shape);

    _inertia_tensor_inverse = _inertia_tensor.inverse();
}

auto spacecraft::inertia_tensor() const -> const mat3x3& {
    return _inertia_tensor;
}

auto spacecraft::inertia_tensor_inverse() const -> const mat3x3& {
    return _inertia_tensor_inverse;
}

auto spacecraft::magnet() const -> const permanent_magnet& {
    return _magnet;
}

auto spacecraft::rods() const -> std::span<const hysteresis_rod> {
    return _rods;
}

void spacecraft::uniform(double mass, const vec3& dim_m, double drag_coefficient) {
    const auto face_surface_area_x = dim_m.y() * dim_m.z();
    const auto face_surface_area_y = dim_m.x() * dim_m.z();
    const auto face_surface_area_z = dim_m.x() * dim_m.y();

    // NOLINTBEGIN(readability-magic-numbers)
    _faces[0].surface_area_m2      = face_surface_area_x;
    _faces[0].surface_normal       = vec3(1, 0, 0);
    _faces[0].center_of_pressure_m = _faces[0].surface_normal * dim_m.x() * 0.5;
    _faces[0].drag_coefficient     = drag_coefficient;
    _faces[1].surface_area_m2      = face_surface_area_x;
    _faces[1].surface_normal       = vec3(-1, 0, 0);
    _faces[1].center_of_pressure_m = _faces[1].surface_normal * dim_m.x() * 0.5;
    _faces[1].drag_coefficient     = drag_coefficient;
    _faces[2].surface_area_m2      = face_surface_area_y;
    _faces[2].surface_normal       = vec3(0, 1, 0);
    _faces[2].center_of_pressure_m = _faces[2].surface_normal * dim_m.y() * 0.5;
    _faces[2].drag_coefficient     = drag_coefficient;
    _faces[3].surface_area_m2      = face_surface_area_y;
    _faces[3].surface_normal       = vec3(0, -1, 0);
    _faces[3].center_of_pressure_m = _faces[3].surface_normal * dim_m.y() * 0.5;
    _faces[3].drag_coefficient     = drag_coefficient;
    _faces[4].surface_area_m2      = face_surface_area_z;
    _faces[4].surface_normal       = vec3(0, 0, 1);
    _faces[4].center_of_pressure_m = _faces[4].surface_normal * dim_m.z() * 0.5;
    _faces[4].drag_coefficient     = drag_coefficient;
    _faces[5].surface_area_m2      = face_surface_area_z;
    _faces[5].surface_normal       = vec3(0, 0, -1);
    _faces[5].center_of_pressure_m = _faces[5].surface_normal * dim_m.z() * 0.5;
    _faces[5].drag_coefficient     = drag_coefficient;
    // NOLINTEND(readability-magic-numbers)

    _inertia_tensor = compute_inertia_tensor(mass, dim_m.x(), dim_m.y(), dim_m.z());
}

auto spacecraft::compute_inertia_tensor(double m, double a, double b, double c) -> mat3x3 {
    const double i_x = (1. / 12.) * m * (b * b + c * c);
    const double i_y = (1. / 12.) * m * (a * a + c * c);
    const double i_z = (1. / 12.) * m * (a * a + b * b);

    mat3x3 tensor = mat3x3::Zero();
    tensor(0, 0)  = i_x;
    tensor(1, 1)  = i_y;
    tensor(2, 2)  = i_z;
    return tensor;
}

}  // namespace aos
