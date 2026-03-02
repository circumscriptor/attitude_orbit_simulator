#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/hysteresis_rods.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft_faces.hpp"
#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

namespace aos {

struct spacecraft_properties {
    double                      mass_kg;
    spacecraft_shape            shape;
    hysteresis_parameters       hysteresis;
    permanent_magnet_properties magnet;
    hysteresis_rods_properties  rods;

    void from_toml(const toml::table& table);
    void debug_print() const;
};

class spacecraft {
public:

    explicit spacecraft(const spacecraft_properties& properties);

    [[nodiscard]] auto mass_kg() const -> double;
    [[nodiscard]] auto inertia_tensor_kg_m2() const -> const mat3x3&;
    [[nodiscard]] auto inertia_tensor_kg_m2_inverse() const -> const mat3x3&;
    [[nodiscard]] auto faces() const -> const spacecraft_faces&;
    [[nodiscard]] auto magnet() const -> const permanent_magnet&;
    [[nodiscard]] auto hystresis() const -> const hysteresis_rods&;

    void derivative(const environment_effects& env, double earth_mu, const system_state& current_state, system_state& state_derivative) const;

    // sums permanent magnet, gyroscopic, and gravity gradient torques
    [[nodiscard]] auto compute_torques(const vec3& omega, const vec3& b_body, const vec3& r_body, double earth_mu) const -> vec3;

    // compute gyroscopic torque [-omega × (I * ω)]
    [[nodiscard]] auto compute_gyroscopic_torque(const vec3& omega) const -> vec3;

    // compute gravity gradient torque
    [[nodiscard]] auto compute_gravity_gradient_torque(const vec3& r_body, double earth_mu) const -> vec3;

protected:

    [[nodiscard]] static auto compute_inertia_tensor(double m_kg, double a, double b, double c) -> mat3x3;

private:

    double           _mass_kg;                       // [kg] Mass
    mat3x3           _inertia_tensor_kg_m2;          // [kg*m^2] Inertia
    mat3x3           _inertia_tensor_kg_m2_inverse;  // [1/(kg*m^2)] Inverse inertia
    spacecraft_faces _faces;
    permanent_magnet _magnet;
    hysteresis_rods  _hystresis;
};

}  // namespace aos
