#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <array>
#include <cstddef>
#include <span>
#include <variant>
#include <vector>

namespace aos {

struct spacecraft_face {
    static constexpr size_t num_faces = 6;  // Cube

    vec3   center_of_pressure_m;                        //!< [m] Face center relative to spacecraft center
    vec3   surface_normal;                              //!< [-] Face normal (outward)
    double surface_area_m2{};                           //!< [m^2] Face surface area
    double drag_coefficient{default_drag_coefficient};  //!< [-] Face drag coefficient

    void from_toml(const toml::table& table);

    void debug_print() const;

    [[nodiscard]] auto compute_force_drag(double density, const quat& q_att, const vec3& v_eci, const vec3& omega_body) const -> vec3;

    // compute relative velocity of spacecraft face center of mass
    [[nodiscard]] auto compute_v_rel_atmosphere(const quat& q_att, const vec3& v_eci, const vec3& omega_body) const -> vec3;
};

using spacecraft_faces = std::array<spacecraft_face, spacecraft_face::num_faces>;

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

    [[nodiscard]] auto mass() const -> double;
    [[nodiscard]] auto inertia_tensor() const -> const mat3x3&;
    [[nodiscard]] auto inertia_tensor_inverse() const -> const mat3x3&;
    [[nodiscard]] auto faces() const -> std::span<const spacecraft_face>;
    [[nodiscard]] auto magnet() const -> const permanent_magnet&;
    [[nodiscard]] auto rods() const -> std::span<const hysteresis_rod>;

protected:

    void uniform(double mass, const vec3& dim_m, double drag_coefficient);

    [[nodiscard]] static auto compute_inertia_tensor(double m, double a, double b, double c) -> mat3x3;

private:

    double                      _mass_g;
    mat3x3                      _inertia_tensor;          // I
    mat3x3                      _inertia_tensor_inverse;  // I^(-1)
    spacecraft_faces            _faces;
    permanent_magnet            _magnet;
    std::vector<hysteresis_rod> _rods;
};

}  // namespace aos
