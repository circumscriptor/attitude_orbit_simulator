#pragma once

#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/types.hpp"

namespace aos {

class inertia_tensor {
public:

    inertia_tensor(real_t mass_kg, const spacecraft_shape& shape);
    inertia_tensor(real_t mass_kg, const spacecraft_uniform& shape);
    explicit inertia_tensor(const spacecraft_custom& shape);

    [[nodiscard]] auto value() const -> const mat3x3&;
    [[nodiscard]] auto inverse() const -> const mat3x3&;

    // compute gyroscopic torque [-omega × (I * ω)]
    [[nodiscard]] auto compute_gyroscopic_torque(const vec3& omega) const -> vec3;

    // compute gravity gradient torque
    [[nodiscard]] auto compute_gravity_gradient_torque(const vec3& r_body, real_t earth_mu) const -> vec3;

protected:

    [[nodiscard]] static auto compute_inertia_tensor(real_t m_kg, real_t a, real_t b, real_t c) -> mat3x3;

private:

    mat3x3 _inertia_tensor_kg_m2;          // [kg*m^2] Inertia
    mat3x3 _inertia_tensor_kg_m2_inverse;  // [1/(kg*m^2)] Inverse inertia
};

}  // namespace aos
