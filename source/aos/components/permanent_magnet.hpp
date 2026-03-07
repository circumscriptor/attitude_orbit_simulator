#pragma once

#include "aos/core/types.hpp"

#include <variant>

namespace aos {

struct permanent_magnet_cylindrical {
    real_t length_m;
    real_t radius_m;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

struct permanent_magnet_rectangular {
    real_t width_m;
    real_t height_m;
    real_t length_m;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

using permanent_magnet_shape = std::variant<permanent_magnet_cylindrical, permanent_magnet_rectangular>;

struct permanent_magnet_properties {
    real_t                 remanence_t;
    real_t                 relative_permeability;
    vec3                   orientation;
    permanent_magnet_shape shape;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class permanent_magnet {
public:

    static constexpr real_t default_temperature_coefficient = -0.0002;  // [-] Default temperature coefficient
    static constexpr real_t default_temperature_reference   = 20.0;     // [deg C] Default reference temperature

    explicit permanent_magnet(const permanent_magnet_properties& properties);

    [[nodiscard]] auto magnetic_moment() const -> vec3;

    // compute torque (m x B)
    [[nodiscard]] auto compute_torque(const vec3& b_field_body) const -> vec3;

    // compute energy
    [[nodiscard]] auto compute_potential_energy(const vec3& b_field_body) const -> real_t;

    // compute magnetic moment using temperature coefficient
    [[nodiscard]] auto compute_magnetic_moment_at_temperature(real_t temp_celsius,
                                                              real_t temp_coeff = default_temperature_coefficient,
                                                              real_t temp_ref   = default_temperature_reference) -> vec3;

protected:

    void set_rectangular(const permanent_magnet_rectangular& shape);
    void set_cylindrical(const permanent_magnet_cylindrical& shape);

    // compute magnetic moment from the state of permanent magnet
    [[nodiscard]] static auto compute_magnetic_moment(real_t      remanence,
                                                      real_t      volume,
                                                      real_t      demagnitization_factor,
                                                      real_t      relative_permeability,
                                                      const vec3& orientation) -> vec3;

private:

    real_t _remanence_t{};             // [T] Remanent magnetization (Br)
    real_t _volume_m3{};               // [m^3] Material volume
    real_t _demagnetization_factor{};  // [-] Geometry-dependent factor (N)
    real_t _relative_permeability{};   // [-] Material permeability (mu_r)
    vec3   _orientation_body;          // [-] Unit vector of magnetization direction
    vec3   _magnetic_moment_body;      // [Am^2] Total magnetic dipole moment (m)
};

}  // namespace aos
