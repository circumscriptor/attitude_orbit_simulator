#include "spacecraft_shape.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/toml.hpp>

#include <cstddef>
#include <iostream>

namespace aos {

// NOLINTBEGIN(readability-magic-numbers)

void spacecraft_uniform::from_toml(const toml_table& table) {
    if (const auto* vec = table["dimensions"].as_array()) {
        dimensions_m <<                                      //
            vec->get(0)->value_or(default_spacecraft_size),  //
            vec->get(1)->value_or(default_spacecraft_size),  //
            vec->get(2)->value_or(default_spacecraft_size);
    }

    drag_coefficient                = table["drag_coeffcient"].value_or(default_drag_coefficient);
    specular_reflection_coefficient = table["specular_reflection_coefficient"].value_or(0.0);
    diffuse_reflection_coefficient  = table["diffuse_reflection_coefficient"].value_or(0.0);
}

void spacecraft_uniform::debug_print() const {
    std::cout << "--  spacecraft (uniform) shape  --"                                                                               //
              << "\n  dimensions:                      " << dimensions_m.x() << ' ' << dimensions_m.y() << ' ' << dimensions_m.z()  //
              << "\n  drag coefficient:                " << drag_coefficient                                                        //
              << "\n  specular reflection coefficient: " << specular_reflection_coefficient                                         //
              << "\n  diffuse reflection coefficient:  " << diffuse_reflection_coefficient                                          //
              << '\n';
}

void spacecraft_custom::from_toml(const toml_table& table) {
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

// NOLINTEND(readability-magic-numbers)

}  // namespace aos
