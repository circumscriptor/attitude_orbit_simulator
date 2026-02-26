#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/types.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace aos {

struct spacecraft_face {
    vec3   center_of_pressure_m;  //!< Face center relative to spacecraft center
    vec3   surface_normal;        //!< Face normal (outward)
    double surface_area_m2{};     //!< Face surface area
};

class spacecraft {
public:

    static constexpr size_t num_faces = 6;  // Cube

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

    [[nodiscard]] auto inertia_tensor() const -> const mat3x3& { return _inertia_tensor; }
    [[nodiscard]] auto inertia_tensor_inverse() const -> const mat3x3& { return _inertia_tensor_inverse; }
    [[nodiscard]] auto magnet() const -> const permanent_magnet& { return _magnet; }
    [[nodiscard]] auto rods() const -> std::span<const hysteresis_rod> { return _rods; }

    static auto get_inertia_tensor(double m, double a, double b, double c) -> mat3x3;

protected:

    void set_faces(const properties& properties);

private:

    mat3x3                                 _inertia_tensor;          // I
    mat3x3                                 _inertia_tensor_inverse;  // I^(-1)
    permanent_magnet                       _magnet;
    std::vector<hysteresis_rod>            _rods;
    std::array<spacecraft_face, num_faces> _faces;
};

}  // namespace aos
