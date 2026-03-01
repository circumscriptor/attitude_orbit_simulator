#pragma once

#include "spacecraft_face.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/types.hpp"

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
    double drag_coefficient;

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
    double                                 mass_g;
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

    [[nodiscard]] auto mass_g() const -> double;
    [[nodiscard]] auto inertia_tensor_kg_m2() const -> const mat3x3&;
    [[nodiscard]] auto inertia_tensor_kg_m2_inverse() const -> const mat3x3&;
    [[nodiscard]] auto faces() const -> std::span<const spacecraft_face>;
    [[nodiscard]] auto magnet() const -> const permanent_magnet&;
    [[nodiscard]] auto rods() const -> std::span<const hysteresis_rod>;

protected:

    void uniform(double mass_kg, const vec3& dim_m, double drag_coefficient);

    [[nodiscard]] static auto compute_inertia_tensor(double m_kg, double a, double b, double c) -> mat3x3;

private:

    double                      _mass_g;                        // [g] Mass
    mat3x3                      _inertia_tensor_kg_m2;          // [kg*m^2] Pre-computed for physics
    mat3x3                      _inertia_tensor_kg_m2_inverse;  // [1/(kg*m^2)] Pre-computed
    spacecraft_faces            _faces;
    permanent_magnet            _magnet;
    std::vector<hysteresis_rod> _rods;
};

}  // namespace aos
