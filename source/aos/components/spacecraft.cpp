#include "spacecraft.hpp"

aos::components::spacecraft::spacecraft(const mat3x3& inertia, const properties& properties)
    : m_inertia_tensor(inertia),
      m_inertia_tensor_inverse(inertia.inverse()),
      m_magnet(properties.magnet_remanence, properties.magnet_length, properties.magnet_diameter, properties.magnet_orientation) {
    m_rods.reserve(properties.hysteresis_rod_orientations.size());
    for (const auto& orientation : properties.hysteresis_rod_orientations) {
        m_rods.emplace_back(properties.hysteresis_rod_volume, orientation, properties.hysteresis_params);
    }
}

aos::components::mat3x3 aos::components::spacecraft::get_inertia_tensor(double m, double a, double b, double c) {
    const double i_x = (1. / 12.) * m * (b * b + c * c);
    const double i_y = (1. / 12.) * m * (a * a + c * c);
    const double i_z = (1. / 12.) * m * (a * a + b * b);

    mat3x3 tensor = mat3x3::Zero();
    tensor(0, 0)  = i_x;
    tensor(1, 1)  = i_y;
    tensor(2, 2)  = i_z;
    return tensor;
}
