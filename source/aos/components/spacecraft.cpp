#include "spacecraft.hpp"

#include "aos/components/permanent_magnet.hpp"
#include "aos/core/types.hpp"

#include <array>
#include <iostream>

namespace aos {

void spacecraft::properties::debug_print() const {
    std::cout << "--  spacecraft properties  --"                                                                                                              //
              << "\n  mass:                                    " << mass_g                                                                                    //
              << "\n  dimensions:                              " << dim_m.x() << ' ' << dim_m.y() << ' ' << dim_m.z()                                         //
              << "\n  magnet orientation:                      " << magnet_orientation.x() << ' ' << magnet_orientation.y() << ' ' << magnet_orientation.z()  //
              << "\n  magnet remanence:                        " << magnet_remanence                                                                          //
              << "\n  magnet length:                           " << magnet_length                                                                             //
              << "\n  magnet diameter:                         " << magnet_diameter                                                                           //
              << "\n  hysteresis rod volume:                   " << hysteresis_rod_volume                                                                     //
              << "\n  hysteresis rod orientations:             ";                                                                                             //

    for (const auto& orientation : hysteresis_rod_orientations) {
        std::cout << '[' << orientation.x() << ' ' << orientation.y() << ' ' << orientation.z() << "] ";
    }
    std::cout << '\n';

    hysteresis_params.debug_print();
}

spacecraft::spacecraft(const mat3x3& inertia, const properties& properties)
    : _inertia_tensor(inertia),
      _inertia_tensor_inverse(inertia.inverse()),
      _magnet(permanent_magnet::cylindrical(properties.magnet_remanence, properties.magnet_length, properties.magnet_diameter, properties.magnet_orientation)) {
    _rods.reserve(properties.hysteresis_rod_orientations.size());
    for (const auto& orientation : properties.hysteresis_rod_orientations) {
        _rods.emplace_back(properties.hysteresis_rod_volume, orientation, properties.hysteresis_params);
    }
    set_faces(properties);
}

auto spacecraft::get_inertia_tensor(double m, double a, double b, double c) -> mat3x3 {
    const double i_x = (1. / 12.) * m * (b * b + c * c);
    const double i_y = (1. / 12.) * m * (a * a + c * c);
    const double i_z = (1. / 12.) * m * (a * a + b * b);

    mat3x3 tensor = mat3x3::Zero();
    tensor(0, 0)  = i_x;
    tensor(1, 1)  = i_y;
    tensor(2, 2)  = i_z;
    return tensor;
}

void spacecraft::set_faces(const properties& properties) {
    const auto face_surface_area_x = properties.dim_m.y() * properties.dim_m.z();
    const auto face_surface_area_y = properties.dim_m.x() * properties.dim_m.z();
    const auto face_surface_area_z = properties.dim_m.x() * properties.dim_m.y();

    // NOLINTBEGIN(readability-magic-numbers)
    _faces[0].surface_area       = face_surface_area_x;
    _faces[0].surface_normal     = vec3(1, 0, 0);
    _faces[0].center_of_pressure = _faces[0].surface_normal * properties.dim_m.x() * 0.5;
    _faces[1].surface_area       = face_surface_area_x;
    _faces[1].surface_normal     = vec3(-1, 0, 0);
    _faces[1].center_of_pressure = _faces[1].surface_normal * properties.dim_m.x() * 0.5;
    _faces[2].surface_area       = face_surface_area_y;
    _faces[2].surface_normal     = vec3(0, 1, 0);
    _faces[2].center_of_pressure = _faces[2].surface_normal * properties.dim_m.y() * 0.5;
    _faces[3].surface_area       = face_surface_area_y;
    _faces[3].surface_normal     = vec3(0, -1, 0);
    _faces[3].center_of_pressure = _faces[3].surface_normal * properties.dim_m.y() * 0.5;
    _faces[4].surface_area       = face_surface_area_z;
    _faces[4].surface_normal     = vec3(0, 0, 1);
    _faces[4].center_of_pressure = _faces[4].surface_normal * properties.dim_m.z() * 0.5;
    _faces[5].surface_area       = face_surface_area_z;
    _faces[5].surface_normal     = vec3(0, 0, -1);
    _faces[5].center_of_pressure = _faces[5].surface_normal * properties.dim_m.z() * 0.5;
    // NOLINTEND(readability-magic-numbers)
}

}  // namespace aos
