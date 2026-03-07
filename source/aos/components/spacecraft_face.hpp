#pragma once

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <array>
#include <cstddef>

namespace aos {

static constexpr size_t spacecraft_num_faces = 6;  // Cube

struct face_effects {
    vec3 torque_body;  //!< [Nm] Torque in body frame
    vec3 force_eci;    //!< [N] Force applied to satellite in ECI
};

struct face_forces {
    vec3 force_drag_body;
    vec3 force_srp_body;
};

struct face_effects_with_forces {
    vec3 torque_body;  //!< [Nm] Torque in body frame
    vec3 force_body;   //!< [N] Force applied to satellite in body frame

    std::array<face_forces, spacecraft_num_faces> forces_body;  // Per-face forces in body frame
};

struct spacecraft_face {
    vec3 center_of_pressure_m;                        //!< [m] Face center relative to spacecraft center
    vec3 surface_normal;                              //!< [-] Face normal (outward)
    real surface_area_m2{};                           //!< [m^2] Face surface area
    real drag_coefficient{default_drag_coefficient};  //!< [-] Face drag coefficient
    real specular_reflection_coefficient{};           //!< [-] Specular reflection coefficient (0..1)
    real diffuse_reflection_coefficient{};            //!< [-] Diffuse reflection coefficient (0..1)

    void from_toml(const toml_table& table);

    void debug_print() const;

    [[nodiscard]] auto compute_force(const environment_effects& data, const vec3& v_body, const vec3& s_body, const vec3& omega_body) const -> vec3;
    [[nodiscard]] auto compute_forces(const environment_effects& data, const vec3& v_body, const vec3& s_body, const vec3& omega_body) const -> face_forces;
    [[nodiscard]] auto compute_force_srp_body(real pressure, const vec3& s_body, real shadow_factor) const -> vec3;
    [[nodiscard]] auto compute_force_drag_body(real density, const vec3& v_rel_body) const -> vec3;
    [[nodiscard]] auto compute_v_rel_body(const vec3& v_com_body, const vec3& omega_body) const -> vec3;
};

using spacecraft_face_array = std::array<spacecraft_face, spacecraft_num_faces>;

}  // namespace aos
