#pragma once

#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/environment/nrlmsise.hpp"

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include <vector>

namespace aos {

class environment_impl : public environment {
public:

    environment_impl(const environment_impl&)                    = delete;
    environment_impl(environment_impl&&)                         = delete;
    auto operator=(const environment_impl&) -> environment_impl& = delete;
    auto operator=(environment_impl&&) -> environment_impl&      = delete;

    explicit environment_impl(const environment_properties& properties);
    ~environment_impl() override;

    [[nodiscard]] auto compute_effects(double t_sec, const vec3& r_eci_m, const vec3& v_eci_m_s) const -> environment_effects override;
    [[nodiscard]] auto earth_mu() const -> double override;

protected:

    // avoid re-allocation
    struct computation_cache {
        std::vector<double> rotation_matrix_buffer;

        // intermediate matrices
        mat3x3 R_ecef_to_eci;  // NOLINT(readability-identifier-naming)
        mat3x3 R_enu_to_ecef;  // NOLINT(readability-identifier-naming)
        mat3x3 R_enu_to_eci;   // NOLINT(readability-identifier-naming)

        // intermediate coordinates
        double lat_deg{};
        double lon_deg{};
        double alt_m{};
        vec3   r_ecef_m;
        vec3   r_sun_eci;  // sun position in ECI

        // other
        double current_year{};
    };

    /** Compute atmospheric density at cached transform */
    [[nodiscard]] auto atmospheric_density() const -> double;

    /** Compute magnetic fields at cached transform */
    [[nodiscard]] auto magnetic_field() const -> vec3;

    /** Compute gravitational fields at cached transform */
    [[nodiscard]] auto gravitational_field() const -> vec3;

    /** Cache coordinate transformation results and matrices */
    void cache_transform(double t_sec, const vec3& r_eci_m) const;

    [[nodiscard]] static auto earth_relative_v(const vec3& v_eci_m_s, const vec3& r_eci_m) -> vec3;

    [[nodiscard]] static auto sun_position_eci(double days_since_j2000) -> vec3;

    [[nodiscard]] static auto solar_perturbation(const vec3& r_sat_eci, const vec3& r_sun_eci) -> vec3;

    [[nodiscard]] static auto earth_shadow_factor(const vec3& r_sat, const vec3& r_sun) -> double;

private:

    double                       _start_year_decimal;
    mutable computation_cache    _cache;
    GeographicLib::Geocentric    _earth;
    GeographicLib::GravityModel  _gravity_model;
    GeographicLib::MagneticModel _magnetic_model;
    nrlmsise                     _atmospheric_model;
};

}  // namespace aos
