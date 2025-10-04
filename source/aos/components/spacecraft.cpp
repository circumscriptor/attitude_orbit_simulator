#include "spacecraft.hpp"

#include "aos/core/types.hpp"

namespace aos::components {

spacecraft::spacecraft(const mat3x3& inertia, const properties& properties)
    : _inertia_tensor(inertia),
      _inertia_tensor_inverse(inertia.inverse()),
      _magnet(properties.magnet_remanence, properties.magnet_length, properties.magnet_diameter, properties.magnet_orientation) {
    _rods.reserve(properties.hysteresis_rod_orientations.size());
    for (const auto& orientation : properties.hysteresis_rod_orientations) {
        _rods.emplace_back(properties.hysteresis_rod_volume, orientation, properties.hysteresis_params);
    }
}

mat3x3 spacecraft::get_inertia_tensor(double m, double a, double b, double c) {
    const double i_x = (1. / 12.) * m * (b * b + c * c);
    const double i_y = (1. / 12.) * m * (a * a + c * c);
    const double i_z = (1. / 12.) * m * (a * a + b * b);

    mat3x3 tensor = mat3x3::Zero();
    tensor(0, 0)  = i_x;
    tensor(1, 1)  = i_y;
    tensor(2, 2)  = i_z;
    return tensor;
}

}  // namespace aos::components
