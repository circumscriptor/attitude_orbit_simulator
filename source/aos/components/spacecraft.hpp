#pragma once

#include "spacecraft_face.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <span>
#include <variant>
#include <vector>

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
    mat3x3           inertia;
    spacecraft_faces faces;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

using spacecraft_shape = std::variant<spacecraft_uniform, spacecraft_custom>;

struct spacecraft_properties {
    double                                 mass_kg;
    spacecraft_shape                       shape;
    hysteresis_parameters                  hysteresis;
    permanent_magnet_properties            magnet;
    std::vector<hysteresis_rod_properties> rods;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

class spacecraft {
public:

    explicit spacecraft(const spacecraft_properties& properties);

    [[nodiscard]] auto mass_kg() const -> double;
    [[nodiscard]] auto inertia_tensor_kg_m2() const -> const mat3x3&;
    [[nodiscard]] auto inertia_tensor_kg_m2_inverse() const -> const mat3x3&;
    [[nodiscard]] auto faces() const -> std::span<const spacecraft_face>;
    [[nodiscard]] auto magnet() const -> const permanent_magnet&;
    [[nodiscard]] auto rods() const -> std::span<const hysteresis_rod>;

    void derivative(const environment_effects& env, double earth_mu, const system_state& current_state, system_state& state_derivative) const;

    // sums permanent magnet, gyroscopic, and gravity gradient torques
    [[nodiscard]] auto compute_torques(const vec3& omega, const vec3& b_body, const vec3& r_body, double earth_mu) const -> vec3;

    // compute gyroscopic torque [-omega × (I * ω)]
    [[nodiscard]] auto compute_gyroscopic_torque(const vec3& omega) const -> vec3;

    // compute gravity gradient torque
    [[nodiscard]] auto compute_gravity_gradient_torque(const vec3& r_body, double earth_mu) const -> vec3;

    // compute drag and srp torque+force
    [[nodiscard]] auto compute_face_effects(const environment_effects& env, const quat& q_att, const quat& q_inv, const vec3& omega_body) const
        -> spacecraft_face_effects;

    // compute total rod torque and dM/dt for each rod, return the total torque exerted by all rods, write dM/dt values into the dm_dt_out
    [[nodiscard]] auto compute_rod_effects(const vecX& rod_magnetizations, const vec3& b_body, const vec3& b_dot_body, vecX& dm_dt_out) const -> vec3;

protected:

    void uniform(double mass_kg, const spacecraft_uniform& shape);

    [[nodiscard]] static auto compute_inertia_tensor(double m_kg, double a, double b, double c) -> mat3x3;

private:

    double                      _mass_kg;                       // [kg] Mass
    mat3x3                      _inertia_tensor_kg_m2;          // [kg*m^2] Inertia
    mat3x3                      _inertia_tensor_kg_m2_inverse;  // [1/(kg*m^2)] Inverse inertia
    spacecraft_faces            _faces;
    permanent_magnet            _magnet;
    std::vector<hysteresis_rod> _rods;
};

}  // namespace aos
