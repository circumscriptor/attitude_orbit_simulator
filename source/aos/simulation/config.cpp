#include "config.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/constants.hpp"

#include <numbers>

namespace aos::simulation {

// {{"N35", 1.21},  // Using nominal Br in Tesla
//  {"N42", 1.32},
//  {"N52", 1.45},
//  {"N35SH", 1.19}};

simulation_parameters simulation_parameters::get_default() {
    // NOLINTBEGIN(readability-magic-numbers)
    return {
        .spacecraft =
            {
                .mass_g                = 1300.0,                                  // 1.3 kg
                .dim_m                 = {0.1, 0.1, 0.1},                         // 10x10x10 cm
                .magnet_orientation    = {0.0, 0.0, 1.0},                         // Permanent Magnet: A Grade N35 NdFeB magnet
                .magnet_remanence      = 1.21,                                    // [T] for Grade N35
                .magnet_length         = 0.05,                                    // 5 cm long
                .magnet_diameter       = 0.01,                                    // 1 cm diameter
                .hysteresis_rod_volume = 0.005 * 0.005 * std::numbers::pi * 0.1,  // 0.5 cm radius, 10 cm length
                .hysteresis_rod_orientations =
                    {
                        {1.0, 0.0, 0.0},
                        {-1.0, 0.0, 0.0},  // One pair along the X-axis
                        {0.0, 1.0, 0.0},
                        {0.0, -1.0, 0.0}  // One pair along the Y-axis
                    },
                .hysteresis_params = components::hysteresis_rod::ja_parameters::hymu80(),
            },
        .environment =
            {
                .orbit_altitude_km     = earth_radius_km + 650.,  // 650 km above surface
                .orbit_inclination_deg = 51.,                     // 51 degree inclicanation
            },
        .initial_angular_velocity = {0.1, -0.05, 0.08},              // rad/s
        .t_start                  = 0.0,                             //
        .t_end                    = 2.0 * 7.0 * 24.0 * 60.0 * 60.0,  // 2 weeks
        .dt_initial               = 0.1,                             //
    };
    // NOLINTEND(readability-magic-numbers)
}

}  // namespace aos::simulation
