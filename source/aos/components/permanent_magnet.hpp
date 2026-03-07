#pragma once

#include "aos/core/types.hpp"

#include <variant>

namespace aos {

struct permanent_magnet_cylindrical {
    real length_m;
    real radius_m;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

struct permanent_magnet_rectangular {
    real width_m;
    real height_m;
    real length_m;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

using permanent_magnet_shape = std::variant<permanent_magnet_cylindrical, permanent_magnet_rectangular>;

struct permanent_magnet_properties {
    real                   remanence_t;
    real                   relative_permeability;
    vec3                   orientation;
    permanent_magnet_shape shape;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class permanent_magnet {
public:

    static constexpr real default_temperature_coefficient = -0.0002;  // [-] Default temperature coefficient
    static constexpr real default_temperature_reference   = 20.0;     // [deg C] Default reference temperature

    explicit permanent_magnet(const permanent_magnet_properties& properties);

    [[nodiscard]] auto magnetic_moment() const -> vec3;

    // compute torque (m x B)
    [[nodiscard]] auto compute_torque(const vec3& b_field_body) const -> vec3;

    // compute energy
    [[nodiscard]] auto compute_potential_energy(const vec3& b_field_body) const -> real;

    // compute magnetic moment using temperature coefficient
    [[nodiscard]] auto compute_magnetic_moment_at_temperature(real temp_celsius,
                                                              real temp_coeff = default_temperature_coefficient,
                                                              real temp_ref   = default_temperature_reference) -> vec3;

protected:

    void set_rectangular(const permanent_magnet_rectangular& shape);
    void set_cylindrical(const permanent_magnet_cylindrical& shape);

    // compute magnetic moment from the state of permanent magnet
    [[nodiscard]] static auto compute_magnetic_moment(real        remanence,
                                                      real        volume,
                                                      real        demagnitization_factor,
                                                      real        relative_permeability,
                                                      const vec3& orientation) -> vec3;

private:

    real _remanence_t{};             // [T] Remanent magnetization (Br)
    real _volume_m3{};               // [m^3] Material volume
    real _demagnetization_factor{};  // [-] Geometry-dependent factor (N)
    real _relative_permeability{};   // [-] Material permeability (mu_r)
    vec3 _orientation_body;          // [-] Unit vector of magnetization direction
    vec3 _magnetic_moment_body;      // [Am^2] Total magnetic dipole moment (m)
};

}  // namespace aos
