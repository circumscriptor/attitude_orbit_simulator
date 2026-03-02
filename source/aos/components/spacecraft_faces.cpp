#include "spacecraft_faces.hpp"

#include "aos/components/spacecraft_face.hpp"
#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <cstddef>
#include <type_traits>
#include <variant>

namespace aos {

spacecraft_faces::spacecraft_faces(const spacecraft_shape& shape) {
    std::visit(
        [&](auto&& shape) {
            using shape_type = std::decay_t<decltype(shape)>;

            if constexpr (std::is_same_v<shape_type, spacecraft_uniform>) {
                uniform(shape);
            } else if constexpr (std::is_same_v<shape_type, spacecraft_custom>) {
                _faces = shape.faces;
            }
        },
        shape);
}

spacecraft_faces::spacecraft_faces(const spacecraft_custom& shape) : _faces(shape.faces) {}

spacecraft_faces::spacecraft_faces(const spacecraft_uniform& shape) {
    uniform(shape);
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
auto spacecraft_faces::compute_face_effects(const environment_effects& env, const quat& q_att, const quat& q_inv, const vec3& omega_body) const
    -> face_effects {
    vec3 torque_body_sum = vec3::Zero();
    vec3 force_body_sum  = vec3::Zero();

    const vec3 s_body = (q_inv * env.r_sun_eci).normalized();
    const vec3 v_body = q_inv * env.v_earth_rel;
    for (const auto& face : _faces) {
        const vec3 f_body = face.compute_force(env, v_body, s_body, omega_body);  // drag + srp
        force_body_sum += f_body;
        torque_body_sum += face.center_of_pressure_m.cross(f_body);
    }

    return {
        .torque_body = torque_body_sum,
        .force_eci   = q_att * force_body_sum,
    };
}

auto spacecraft_faces::compute_faces_effects_with_forces(const environment_effects& env, const quat& q_inv, const vec3& omega_body) const
    -> face_effects_with_forces {
    face_effects_with_forces result;

    const vec3 s_body = (q_inv * env.r_sun_eci).normalized();
    const vec3 v_body = q_inv * env.v_earth_rel;
    for (size_t i = 0; i < _faces.size(); ++i) {
        const auto& face   = _faces.at(i);
        const auto  f      = face.compute_forces(env, v_body, s_body, omega_body);  // drag + srp
        const vec3  f_body = f.force_drag_body + f.force_srp_body;
        const vec3  t_body = face.center_of_pressure_m.cross(f_body);

        result.force_body += f_body;
        result.torque_body += t_body;
        result.forces_body.at(i) = f;
    }

    return result;
}
// NOLINTEND(bugprone-easily-swappable-parameters)

void spacecraft_faces::uniform(const spacecraft_uniform& shape) {
    const auto& dim_m                           = shape.dimensions_m;
    const auto& drag_coefficient                = shape.drag_coefficient;
    const auto& specular_reflection_coefficient = shape.specular_reflection_coefficient;
    const auto& diffuse_reflection_coefficient  = shape.diffuse_reflection_coefficient;
    const auto  face_surface_area_x             = dim_m.y() * dim_m.z();
    const auto  face_surface_area_y             = dim_m.x() * dim_m.z();
    const auto  face_surface_area_z             = dim_m.x() * dim_m.y();

    // NOLINTBEGIN(readability-magic-numbers)
    _faces[0].surface_area_m2                 = face_surface_area_x;
    _faces[0].surface_normal                  = vec3(1, 0, 0);
    _faces[0].center_of_pressure_m            = _faces[0].surface_normal * dim_m.x() * 0.5;
    _faces[0].drag_coefficient                = drag_coefficient;
    _faces[0].specular_reflection_coefficient = specular_reflection_coefficient;
    _faces[0].diffuse_reflection_coefficient  = diffuse_reflection_coefficient;
    _faces[1].surface_area_m2                 = face_surface_area_x;
    _faces[1].surface_normal                  = vec3(-1, 0, 0);
    _faces[1].center_of_pressure_m            = _faces[1].surface_normal * dim_m.x() * 0.5;
    _faces[1].drag_coefficient                = drag_coefficient;
    _faces[1].specular_reflection_coefficient = specular_reflection_coefficient;
    _faces[1].diffuse_reflection_coefficient  = diffuse_reflection_coefficient;
    _faces[2].surface_area_m2                 = face_surface_area_y;
    _faces[2].surface_normal                  = vec3(0, 1, 0);
    _faces[2].center_of_pressure_m            = _faces[2].surface_normal * dim_m.y() * 0.5;
    _faces[2].drag_coefficient                = drag_coefficient;
    _faces[2].specular_reflection_coefficient = specular_reflection_coefficient;
    _faces[2].diffuse_reflection_coefficient  = diffuse_reflection_coefficient;
    _faces[3].surface_area_m2                 = face_surface_area_y;
    _faces[3].surface_normal                  = vec3(0, -1, 0);
    _faces[3].center_of_pressure_m            = _faces[3].surface_normal * dim_m.y() * 0.5;
    _faces[3].drag_coefficient                = drag_coefficient;
    _faces[3].specular_reflection_coefficient = specular_reflection_coefficient;
    _faces[3].diffuse_reflection_coefficient  = diffuse_reflection_coefficient;
    _faces[4].surface_area_m2                 = face_surface_area_z;
    _faces[4].surface_normal                  = vec3(0, 0, 1);
    _faces[4].center_of_pressure_m            = _faces[4].surface_normal * dim_m.z() * 0.5;
    _faces[4].drag_coefficient                = drag_coefficient;
    _faces[4].specular_reflection_coefficient = specular_reflection_coefficient;
    _faces[4].diffuse_reflection_coefficient  = diffuse_reflection_coefficient;
    _faces[5].surface_area_m2                 = face_surface_area_z;
    _faces[5].surface_normal                  = vec3(0, 0, -1);
    _faces[5].center_of_pressure_m            = _faces[5].surface_normal * dim_m.z() * 0.5;
    _faces[5].drag_coefficient                = drag_coefficient;
    _faces[5].specular_reflection_coefficient = specular_reflection_coefficient;
    _faces[5].diffuse_reflection_coefficient  = diffuse_reflection_coefficient;
    // NOLINTEND(readability-magic-numbers)
}

}  // namespace aos
