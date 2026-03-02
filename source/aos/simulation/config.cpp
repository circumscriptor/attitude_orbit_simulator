#include "config.hpp"

#include "aos/core/constants.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <iostream>

namespace aos {

// {{"N35", 1.21},  // Using nominal Br in Tesla
//  {"N42", 1.32},
//  {"N52", 1.45},
//  {"N35SH", 1.19}};

void simulation_properties::from_toml(const toml::table& table) {
    // NOLINTBEGIN(readability-magic-numbers)

    if (const auto* sat = table["satellite"].as_table()) {
        satellite.from_toml(*sat);
    }

    if (const auto* orb = table["orbit"].as_table()) {
        orbit.from_toml(*orb);
    }

    if (const auto* obs = table["observer"].as_table()) {
        observer.from_toml(*obs);
    }

    if (const auto* env = table["environment"].as_table()) {
        environment.from_toml(*env);
    }

    if (const auto* vec = table["angular_velocity"].as_array()) {
        angular_velocity <<              //
            vec->get(0)->value_or(0.0),  //
            vec->get(1)->value_or(0.0),  //
            vec->get(2)->value_or(0.0);
    }

    t_start             = table["t_start"].value_or(0.0);
    t_end               = table["t_end"].value_or(2.0 * 7.0 * 24.0 * 60.0 * 60.0);  // 2 weeks
    dt_initial          = table["dt_initial"].value_or(0.1);
    absolute_error      = table["absolute_error"].value_or(default_absolute_error);
    relative_error      = table["relative_error"].value_or(default_relative_error);
    stepper_function    = table["stepper_function"].value_or(0);
    checkpoint_interval = table["checkpoint_interval"].value_or(0.0);

    // NOLINTEND(readability-magic-numbers)
}

void simulation_properties::debug_print() const {
    std::cout << std::boolalpha;

    satellite.debug_print();
    orbit.debug_print();
    environment.debug_print();

    std::cout << "-- simulation properties --";
    std::cout << "\n  angular velocity:    " << angular_velocity.x() << ' ' << angular_velocity.y() << ' ' << angular_velocity.z()  //
              << "\n  t start:             " << t_start                                                                             //
              << "\n  t end:               " << t_end                                                                               //
              << "\n  initial dt:          " << dt_initial                                                                          //
              << "\n  absolute error:      " << absolute_error                                                                      //
              << "\n  relative error:      " << relative_error                                                                      //
              << "\n  stepper function:    " << stepper_function                                                                    //
              << "\n  checkpoint interval: " << checkpoint_interval                                                                 //
              << '\n';

    std::cout << "----\n";
}

}  // namespace aos
