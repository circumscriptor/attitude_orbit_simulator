#pragma once

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <array>
#include <cstddef>

namespace aos {

struct spacecraft_face_effects {
    vec3 torque_body;
    vec3 force_eci;
};

struct spacecraft_face {
    static constexpr size_t num_faces = 6;  // Cube

    vec3   center_of_pressure_m;                        //!< [m] Face center relative to spacecraft center
    vec3   surface_normal;                              //!< [-] Face normal (outward)
    double surface_area_m2{};                           //!< [m^2] Face surface area
    double drag_coefficient{default_drag_coefficient};  //!< [-] Face drag coefficient
    double specular_reflection_coefficient{};           //!< [-] Specular reflection coefficient (0..1)
    double diffuse_reflection_coefficient{};            //!< [-] Diffuse reflection coefficient (0..1)

    void from_toml(const toml::table& table);

    void debug_print() const;

    [[nodiscard]] auto compute_force(const environment_data& data, const vec3& v_body, const vec3& s_body, const vec3& omega_body) const -> vec3;
    [[nodiscard]] auto compute_force_srp_body(double pressure, const vec3& s_body, double shadow_factor) const -> vec3;
    [[nodiscard]] auto compute_force_drag_body(double density, const vec3& v_rel_body) const -> vec3;
    [[nodiscard]] auto compute_v_rel_body(const vec3& v_com_body, const vec3& omega_body) const -> vec3;
};

using spacecraft_faces = std::array<spacecraft_face, spacecraft_face::num_faces>;

}  // namespace aos
