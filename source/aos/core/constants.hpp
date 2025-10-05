#pragma once

#include <numbers>

namespace aos {

static constexpr double earth_mu_m3_s2            = 3.986004418e14;
static constexpr double earth_radius_km           = 6371.;
static constexpr double earth_rotation_rate_rad_s = 7.2921150e-5;
static constexpr double epsilon                   = 1e-6;
static constexpr double km_to_m                   = 1000.0;
static constexpr double nt_to_t                   = 1e-9;  // nanoTesla to Tesla
static constexpr double deg_to_rad                = std::numbers::pi / 180.0;
static constexpr double rad_to_deg                = 180.0 / std::numbers::pi;
static constexpr double seconds_per_year          = 365.25 * 24.0 * 3600.0;
static constexpr double simulation_start_year     = 2025.0;
static constexpr double two_pi                    = 2.0 * std::numbers::pi;
static constexpr double vacuum_permeability       = 1.25663706212e-6;

}  // namespace aos
