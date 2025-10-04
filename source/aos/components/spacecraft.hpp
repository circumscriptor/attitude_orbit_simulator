#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/types.hpp"

#include <span>
#include <vector>

namespace aos::components {

using core::mat3x3;
using core::vec3;

class spacecraft {
public:

    struct properties {
        vec3                          magnet_orientation;
        double                        magnet_remanence{};
        double                        magnet_length{};
        double                        magnet_diameter{};
        double                        hysteresis_rod_volume{};
        std::vector<vec3>             hysteresis_rod_orientations;
        hysteresis_rod::ja_parameters hysteresis_params;
    };

    spacecraft(double mass, double size, const properties& properties) : spacecraft(mass, size, size, size, properties) {}

    spacecraft(double mass, double width, double height, double length, const properties& properties)
        : spacecraft(get_inertia_tensor(mass, width, height, length), properties) {}

    explicit spacecraft(const mat3x3& inertia, const properties& properties);

    [[nodiscard]] const mat3x3&                   inertia_tensor() const { return m_inertia_tensor; }
    [[nodiscard]] const mat3x3&                   inertia_tensor_inverse() const { return m_inertia_tensor_inverse; }
    [[nodiscard]] const permanent_magnet&         magnet() const { return m_magnet; }
    [[nodiscard]] std::span<const hysteresis_rod> rods() const { return m_rods; }

protected:

    static mat3x3 get_inertia_tensor(double m, double a, double b, double c);

private:

    mat3x3                      m_inertia_tensor;
    mat3x3                      m_inertia_tensor_inverse;
    permanent_magnet            m_magnet;
    std::vector<hysteresis_rod> m_rods;
};

}  // namespace aos::components
