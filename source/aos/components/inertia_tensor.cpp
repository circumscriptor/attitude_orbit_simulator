#include "inertia_tensor.hpp"

#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/types.hpp"

#include <cmath>
#include <type_traits>
#include <variant>

namespace aos {

inertia_tensor::inertia_tensor(double mass_kg, const spacecraft_shape& shape) {
    std::visit(
        [&](auto&& shape) {
            using shape_type = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<shape_type, spacecraft_uniform>) {
                _inertia_tensor_kg_m2 = compute_inertia_tensor(mass_kg, shape.dimensions_m.x(), shape.dimensions_m.y(), shape.dimensions_m.z());
            } else if constexpr (std::is_same_v<shape_type, spacecraft_custom>) {
                _inertia_tensor_kg_m2 = shape.inertia;
            }
        },
        shape);
    _inertia_tensor_kg_m2_inverse = _inertia_tensor_kg_m2.inverse();
}

auto inertia_tensor::value() const -> const mat3x3& {
    return _inertia_tensor_kg_m2;
}

auto inertia_tensor::inverse() const -> const mat3x3& {
    return _inertia_tensor_kg_m2_inverse;
}

auto inertia_tensor::compute_gyroscopic_torque(const vec3& omega) const -> vec3 {
    return -omega.cross(_inertia_tensor_kg_m2 * omega);
}

auto inertia_tensor::compute_gravity_gradient_torque(const vec3& r_body, double earth_mu) const -> vec3 {
    const double r_sq  = r_body.squaredNorm();
    const double r_mag = std::sqrt(r_sq);

    // tau_gg = (3 * mu / r^5) * (r_body × (I * r_body))
    const double coef = (3.0 * earth_mu) / (r_sq * r_sq * r_mag);
    return coef * r_body.cross(_inertia_tensor_kg_m2 * r_body);
}

auto inertia_tensor::compute_inertia_tensor(double m_kg, double a, double b, double c) -> mat3x3 {
    const double i_x = (1. / 12.) * m_kg * (b * b + c * c);
    const double i_y = (1. / 12.) * m_kg * (a * a + c * c);
    const double i_z = (1. / 12.) * m_kg * (a * a + b * b);

    mat3x3 tensor = mat3x3::Zero();
    tensor(0, 0)  = i_x;
    tensor(1, 1)  = i_y;
    tensor(2, 2)  = i_z;
    return tensor;
}

}  // namespace aos
