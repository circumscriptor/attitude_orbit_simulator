#include "spacecraft_face.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <cmath>
#include <iostream>
#include <limits>

namespace aos {

void spacecraft_face::from_toml(const toml::table& table) {
    if (const auto* vec = table["center_of_pressure_m"].as_array()) {
        center_of_pressure_m <<          //
            vec->get(0)->value_or(0.0),  //
            vec->get(1)->value_or(0.0),  //
            vec->get(2)->value_or(0.0);
    }

    if (const auto* vec = table["surface_normal"].as_array()) {
        surface_normal <<                //
            vec->get(0)->value_or(0.0),  //
            vec->get(1)->value_or(0.0),  //
            vec->get(2)->value_or(0.0);
    }

    surface_area_m2                 = table["surface_area_m2"].value_or(0.0);
    drag_coefficient                = table["drag_coefficient"].value_or(default_drag_coefficient);
    specular_reflection_coefficient = table["specular_reflection_coefficient"].value_or(0.0);
    diffuse_reflection_coefficient  = table["diffuse_reflection_coefficient"].value_or(0.0);
}

void spacecraft_face::debug_print() const {
    std::cout << "--  spacecraft face  --"                                                                                                                  //
              << "\n  center of pressure:              " << center_of_pressure_m.x() << ' ' << center_of_pressure_m.y() << ' ' << center_of_pressure_m.z()  //
              << "\n  surface normal:                  " << surface_normal.x() << ' ' << surface_normal.y() << ' ' << surface_normal.z()                    //
              << "\n  surface area:                    " << surface_area_m2                                                                                 //
              << "\n  drag coefficient:                " << drag_coefficient                                                                                //
              << "\n  specular reflection coefficient: " << specular_reflection_coefficient                                                                 //
              << "\n  diffuse reflection coefficient:  " << diffuse_reflection_coefficient                                                                  //
              << '\n';
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto spacecraft_face::compute_force(const environment_effects& data, const vec3& v_body, const vec3& s_body, const vec3& omega_body) const -> vec3 {
    const auto [f_drag_body, f_srp_body] = compute_forces(data, v_body, s_body, omega_body);
    return f_drag_body + f_srp_body;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto spacecraft_face::compute_forces(const environment_effects& data, const vec3& v_body, const vec3& s_body, const vec3& omega_body) const -> face_forces {
    const vec3 v_rel_body  = compute_v_rel_body(v_body, omega_body);
    const vec3 f_drag_body = compute_force_drag_body(data.atmospheric_density_kg_m3, v_rel_body);
    const vec3 f_srp_body  = compute_force_srp_body(data.solar_pressure_Pa, s_body, data.shadow_factor);
    return {
        .force_drag_body = f_drag_body,
        .force_srp_body  = f_srp_body,
    };
}

auto spacecraft_face::compute_force_srp_body(double pressure, const vec3& s_body, double shadow_factor) const -> vec3 {
    const double cos_alpha = surface_normal.dot(s_body);
    if (cos_alpha <= 0.0 || shadow_factor <= 0.0) {
        return vec3::Zero();
    }

    const double scalar        = pressure * surface_area_m2 * cos_alpha * shadow_factor;
    const vec3   f_sun_dir     = -s_body * scalar * (1.0 - specular_reflection_coefficient);
    const double normal_scalar = scalar * (2.0 * specular_reflection_coefficient * cos_alpha + (2.0 / 3.0) * diffuse_reflection_coefficient);
    const vec3   f_normal_dir  = -surface_normal * normal_scalar;
    return f_sun_dir + f_normal_dir;
}

auto spacecraft_face::compute_force_drag_body(double density, const vec3& v_rel_body) const -> vec3 {
    const double     v_sq = v_rel_body.squaredNorm();
    constexpr double eps  = std::numeric_limits<double>::epsilon();
    if (v_sq <= (eps * eps)) {
        return vec3::Zero();
    }

    const double v_mag     = std::sqrt(v_sq);
    const double cos_theta = surface_normal.dot(-v_rel_body / v_mag);
    if (!(cos_theta > 0.0)) {
        return vec3::Zero();
    }

    const double scalar = 0.5 * density * drag_coefficient * surface_area_m2 * cos_theta * v_mag;
    return scalar * v_rel_body;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto spacecraft_face::compute_v_rel_body(const vec3& v_com_body, const vec3& omega_body) const -> vec3 {
    return v_com_body + omega_body.cross(center_of_pressure_m);
}

}  // namespace aos
