#include "config.hpp"

#include "aos/components/hysteresis_rod.hpp"
#include "aos/core/constants.hpp"

#include <iostream>
#include <numbers>

namespace aos {

// {{"N35", 1.21},  // Using nominal Br in Tesla
//  {"N42", 1.32},
//  {"N52", 1.45},
//  {"N35SH", 1.19}};

void simulation_parameters::debug_print() const {
    std::cout << std::boolalpha;

    spacecraft.debug_print();

    std::cout << "-- simulation properties --";
    std::cout << "\nangular velocity:                        " << angular_velocity.x() << ' ' << angular_velocity.y() << ' ' << angular_velocity.z()  //
              << "\nsimulation year:                         " << simulation_year                                                                     //
              << "\ngravity model degree:                    " << gravity_model_degree                                                                //
              << "\nt start:                                 " << t_start                                                                             //
              << "\nt end:                                   " << t_end                                                                               //
              << "\ninitial dt:                              " << dt_initial                                                                          //
              << "\nabsolute error:                          " << absolute_error                                                                      //
              << "\nrelative error:                          " << relative_error                                                                      //
              << "\nhigher order:                            " << higher_order                                                                        //
              << "\ncheckpoint interval:                     " << checkpoint_interval                                                                 //
              << '\n';

    std::cout << "----\n";
}

simulation_parameters simulation_parameters::get_default() {
    // NOLINTBEGIN(readability-magic-numbers)
    return {
        .spacecraft =
            {.mass_g                      = 1300.0,                                  // 1.3 kg
             .dim_m                       = {0.1, 0.1, 0.1},                         // 10x10x10 cm
             .magnet_orientation          = {0.0, 0.0, 1.0},                         // Permanent Magnet: A Grade N35 NdFeB magnet
             .magnet_remanence            = 1.21,                                    // [T] for Grade N35
             .magnet_length               = 0.05,                                    // 5 cm long
             .magnet_diameter             = 0.01,                                    // 1 cm diameter
             .hysteresis_rod_volume       = 0.005 * 0.005 * std::numbers::pi * 0.1,  // 0.5 cm radius, 10 cm length
             .hysteresis_rod_orientations = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0}},  // Pairs along the X-axis and Y-axis
             .hysteresis_params           = hysteresis_rod::ja_parameters::hymu80()},
        .orbit                = {.semi_major_axis_m    = 6818137.0,
                                 .eccentricity         = 0.001,
                                 .inclination_rad      = 1.396263,
                                 .raan_rad             = 0.0,
                                 .arg_of_periapsis_rad = 0.0,
                                 .mean_anomaly_rad     = 0.0},  // 440 km above surface, 80 degree inclicanation
        .observer             = {.exclude_elements = false, .exclude_magnitudes = false},
        .angular_velocity     = {0.1, -0.05, 0.08},  // rad/s
        .simulation_year      = 2026.0,
        .gravity_model_degree = 12,
        .t_start              = 0.0,                             //
        .t_end                = 2.0 * 7.0 * 24.0 * 60.0 * 60.0,  // 2 weeks
        .dt_initial           = 0.1,                             //
        .absolute_error       = default_absolute_error,
        .relative_error       = default_relative_error,
        .higher_order         = false,
        .checkpoint_interval  = 0.0,
    };
    // NOLINTEND(readability-magic-numbers)
}

}  // namespace aos
