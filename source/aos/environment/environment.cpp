#include "environment.hpp"

#include "aos/core/types.hpp"
#include "aos/environment/details/environment_impl.hpp"

#include <toml++/toml.hpp>

#include <iostream>
#include <memory>
#include <string>

namespace aos {

void environment_properties::from_toml(const toml_table& table) {
    // NOLINTBEGIN(readability-magic-numbers)

    start_year_decimal    = table["start_year_decimal"].value_or(2026.0);
    gravity_model_name    = table["gravity_model_name"].value_or<std::string>("egm2008");
    gravity_model_path    = table["gravity_model_path"].value_or<std::string>("");
    magnetic_model_name   = table["magnetic_model_name"].value_or<std::string>("wmm2025");
    magnetic_model_path   = table["magnetic_model_path"].value_or<std::string>("");
    weather_data_path     = table["weather_data_path"].value_or<std::string>(WEATHER_DATA_PATH);
    gravity_model_degree  = table["gravity_model_degree"].value_or(-1);
    gravity_model_order   = table["gravity_model_order"].value_or(-1);
    magnetic_model_degree = table["magnetic_model_degree"].value_or(-1);
    magnetic_model_order  = table["magnetic_model_order"].value_or(-1);

    // NOLINTEND(readability-magic-numbers)
}

void environment_properties::debug_print() const {
    std::cout << "--  environmnet properties  --"                        //
              << "\n  start year:            " << start_year_decimal     //
              << "\n  gravity model name:    " << gravity_model_name     //
              << "\n  gravity model path:    " << gravity_model_path     //
              << "\n  magnetic model name:   " << magnetic_model_name    //
              << "\n  magnetic model path:   " << magnetic_model_path    //
              << "\n  weather data path:     " << weather_data_path      //
              << "\n  gravity model degree:  " << gravity_model_degree   //
              << "\n  gravity model order:   " << gravity_model_order    //
              << "\n  magnetic model degree: " << magnetic_model_degree  //
              << "\n  magnetic model order:  " << magnetic_model_order   //
              << '\n';
}

environment::~environment() = default;

auto environment::create(const environment_properties& properties) -> std::shared_ptr<environment> {
    return std::make_shared<environment_impl>(properties);
}

}  // namespace aos
