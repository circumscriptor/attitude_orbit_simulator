#pragma once

#include "aos/core/types.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <variant>

namespace aos {

struct permanent_magnet_cylindrical {
    double length_m;
    double radius_m;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

struct permanent_magnet_rectangular {
    double width_m;
    double height_m;
    double length_m;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

using permanent_magnet_shape = std::variant<permanent_magnet_cylindrical, permanent_magnet_rectangular>;

struct permanent_magnet_properties {
    double                 remanence_t;
    double                 relative_permeability;
    vec3                   orientation;
    permanent_magnet_shape shape;

    void from_toml(const toml::table& table);

    void debug_print() const;
};

class permanent_magnet {
public:

    static constexpr double default_temp_coeff = -0.0002;  // [-] Default temperature coefficient
    static constexpr double default_temp_ref   = 20.0;     // [deg celsius] Default reference temperature

    explicit permanent_magnet(const permanent_magnet_properties& properties);

    // compute magnetic moment based on the current state of the magnet
    [[nodiscard]] auto magnetic_moment() const -> vec3;
    // compute magnetic moment using temperature coefficient
    [[nodiscard]] auto magnetic_moment_at_temperature(double temp_celsius, double temp_coeff = default_temp_coeff, double ref_temp = default_temp_ref) -> vec3;

protected:

    // compute magnetic moment from the state of permanent magnet
    [[nodiscard]] static auto compute_magnetic_moment(double remanence, double vol, double demag, double mu_r, const vec3& orientation) -> vec3;

private:

    double _remanence_t;             // [Tesla]
    double _volume_m3;               // [m^3]
    double _demagnetization_factor;  //
    double _relative_permeability;   //
    vec3   _orientation_body;        // Unit vector
    vec3   _magnetic_moment_body;    // [A·m^2]
};

}  // namespace aos
