#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/types.hpp"

#include <span>
#include <vector>

namespace aos {

class spacecraft {
public:

    struct properties {
        double                        mass_g{};
        vec3                          dim_m;
        vec3                          magnet_orientation;
        double                        magnet_remanence{};
        double                        magnet_length{};
        double                        magnet_diameter{};
        double                        hysteresis_rod_volume{};
        std::vector<vec3>             hysteresis_rod_orientations;
        hysteresis_rod::ja_parameters hysteresis_params;

        void debug_print() const;
    };

    explicit spacecraft(const properties& properties)
        : spacecraft(get_inertia_tensor(properties.mass_g, properties.dim_m.x(), properties.dim_m.y(), properties.dim_m.z()), properties) {}

    spacecraft(const mat3x3& inertia, const properties& properties);

    [[nodiscard]] const mat3x3&                   inertia_tensor() const { return _inertia_tensor; }
    [[nodiscard]] const mat3x3&                   inertia_tensor_inverse() const { return _inertia_tensor_inverse; }
    [[nodiscard]] const permanent_magnet&         magnet() const { return _magnet; }
    [[nodiscard]] std::span<const hysteresis_rod> rods() const { return _rods; }

    static mat3x3 get_inertia_tensor(double m, double a, double b, double c);

private:

    mat3x3                      _inertia_tensor;
    mat3x3                      _inertia_tensor_inverse;
    permanent_magnet            _magnet;
    std::vector<hysteresis_rod> _rods;
};

}  // namespace aos
