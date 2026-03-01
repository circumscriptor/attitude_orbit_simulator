#include "spacecraft_face.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/impl/table.hpp>

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

    surface_area_m2  = table["surface_area_m2"].value_or(0.0);
    drag_coefficient = table["drag_coefficient"].value_or(default_drag_coefficient);
    reflectivity     = table["reflectivity"].value_or(default_reflectivity);
}

void spacecraft_face::debug_print() const {
    std::cout << "--  spacecraft face  --"                                                                                                     //
              << "\n  center of pressure: " << center_of_pressure_m.x() << ' ' << center_of_pressure_m.y() << ' ' << center_of_pressure_m.z()  //
              << "\n  surface normal:     " << surface_normal.x() << ' ' << surface_normal.y() << ' ' << surface_normal.z()                    //
              << "\n  surface area:       " << surface_area_m2                                                                                 //
              << "\n  drag coefficient:   " << drag_coefficient                                                                                //
              << "\n  reflectivity:       " << reflectivity                                                                                    //
              << '\n';
}

auto spacecraft_face::compute_force_srp_body(double pressure, const vec3& s_body, double shadow_factor) const -> vec3 {
    const double cos_alpha = surface_normal.dot(s_body);
    if (cos_alpha <= 0.0 || shadow_factor <= 0.0) {
        return vec3::Zero();
    }
    return s_body * (pressure * surface_area_m2 * cos_alpha * shadow_factor * reflectivity);
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
