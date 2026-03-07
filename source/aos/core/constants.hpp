#pragma once

#include "aos/core/types.hpp"

#include <numbers>

namespace aos {

static constexpr real default_spacecraft_size   = 0.1;    // 100 mm
static constexpr real default_spacecraft_mass   = 1300.;  // 1.3 kg
static constexpr real default_drag_coefficient  = 2.2;
static constexpr real default_reflectivity      = 1.3;
static constexpr real default_absolute_error    = 1e-10;
static constexpr real default_relative_error    = 1e-10;
static constexpr real kilometer_to_meter        = 1000.;
static constexpr real meter_to_kilometer        = 0.001;
static constexpr real au_to_m                   = 149597870700.;
inline constexpr real au_to_m_2                 = au_to_m * au_to_m;
static constexpr real gram_to_kilogram          = 0.001;
static constexpr real pi                        = std::numbers::pi_v<real>;
static constexpr real deg_to_rad                = pi / 180.;
static constexpr real rad_to_deg                = 180.0 / pi;
static constexpr real nanotesla_to_tesla        = 1e-9;
static constexpr real year_to_days              = 365.25;
static constexpr real day_to_seconds            = 24. * 60. * 60.;
static constexpr real earth_mu_m3_s2            = 3.986004418e14;
static constexpr real sun_mu_m3_s2              = 1.32712440018e20;
static constexpr real earth_radius_km           = 6371.;
static constexpr real earth_radius_m            = earth_radius_km * kilometer_to_meter;
static constexpr real reentry_altitude_km       = 100.0;  // Kármán line
static constexpr real reentry_altitude_m        = reentry_altitude_km * kilometer_to_meter;
static constexpr real earth_rotation_rate_rad_s = 7.2921150e-5;
inline constexpr real sun_radius_km             = 696340.;
inline constexpr real sun_radius_m              = sun_radius_km * kilometer_to_meter;
inline constexpr real solar_pressure_1au        = 1367. / 299792458.;  // P/c ≈ 4.56e-6 N/m^2
static constexpr real seconds_per_year          = 365.25 * 24. * 3600.;
static constexpr real simulation_start_year     = 2025.;
static constexpr real two_pi                    = 2. * pi;
static constexpr real vacuum_permeability       = 4. * pi * 1e-7;
static constexpr real dt_gradient               = 0.01;

}  // namespace aos
