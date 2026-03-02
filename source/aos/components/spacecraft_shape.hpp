#pragma once

#include "aos/components/spacecraft_face.hpp"
#include "aos/core/types.hpp"

#include <variant>

namespace aos {

struct spacecraft_uniform {
    vec3   dimensions_m;
    double drag_coefficient{};
    double specular_reflection_coefficient{};
    double diffuse_reflection_coefficient{};

    void from_toml(const toml::table& table);
    void debug_print() const;
};

struct spacecraft_custom {
    mat3x3                inertia;
    spacecraft_face_array faces;

    void from_toml(const toml::table& table);
    void debug_print() const;
};

using spacecraft_shape = std::variant<spacecraft_uniform, spacecraft_custom>;

}  // namespace aos
