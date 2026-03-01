#pragma once

#include <numbers>

namespace aos {

static constexpr double default_spacecraft_size   = 0.1;    // 100 mm
static constexpr double default_spacecraft_mass   = 1300.;  // 1.3 kg
static constexpr double default_drag_coefficient  = 2.2;
static constexpr double default_reflectivity      = 1.3;
static constexpr double default_absolute_error    = 1e-10;
static constexpr double default_relative_error    = 1e-10;
static constexpr double kilometer_to_meter        = 1000.;
static constexpr double meter_to_kilometer        = 0.001;
static constexpr double au_to_m                   = 149597870700.;
inline constexpr double au_to_m_2                 = au_to_m * au_to_m;
static constexpr double gram_to_kilogram          = 0.001;
static constexpr double deg_to_rad                = std::numbers::pi / 180.;
static constexpr double rad_to_deg                = 180.0 / std::numbers::pi;
static constexpr double nanotesla_to_tesla        = 1e-9;
static constexpr double year_to_days              = 365.25;
static constexpr double day_to_seconds            = 24. * 60. * 60.;
static constexpr double earth_mu_m3_s2            = 3.986004418e14;
static constexpr double sun_mu_m3_s2              = 1.32712440018e20;
static constexpr double earth_radius_km           = 6371.;
static constexpr double earth_radius_m            = earth_radius_km * kilometer_to_meter;
static constexpr double earth_rotation_rate_rad_s = 7.2921150e-5;
inline constexpr double sun_radius_km             = 696340.;
inline constexpr double sun_radius_m              = sun_radius_km * kilometer_to_meter;
inline constexpr double solar_pressure_1au        = 1367. / 299792458.;  // P/c ≈ 4.56e-6 N/m^2
static constexpr double seconds_per_year          = 365.25 * 24. * 3600.;
static constexpr double simulation_start_year     = 2025.;
static constexpr double two_pi                    = 2. * std::numbers::pi;
static constexpr double vacuum_permeability       = 4. * std::numbers::pi * 1e-7;
static constexpr double dt_gradient               = 0.01;

}  // namespace aos
