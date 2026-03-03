#include "spacecraft.hpp"

#include "aos/components/details/spacecraft_impl.hpp"
#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/hysteresis_rods.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <toml++/toml.hpp>

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <variant>

namespace aos {

void spacecraft_properties::from_toml(const toml_table& table) {
    mass_kg = table["mass"].value_or(default_spacecraft_mass);

    if (const auto* uniform = table["uniform"].as_table()) {
        shape.emplace<spacecraft_uniform>().from_toml(*uniform);
    } else if (const auto* custom = table["custom"].as_table()) {
        shape.emplace<spacecraft_custom>().from_toml(*custom);
    }

    if (const auto* hyst = table["hysteresis"].as_table()) {
        hysteresis.from_toml(*hyst);
    }

    if (const auto* mag = table["magnet"].as_table()) {
        magnet.from_toml(*mag);
    }

    if (const auto* arr = table["rods"].as_array()) {
        rods.clear();
        rods.reserve(arr->size());
        for (size_t i = 0; i < arr->size(); ++i) {
            if (const auto* rod = arr->get(i)->as_table()) {
                rods.emplace_back().from_toml(*rod);
            }
        }
    }
}

void spacecraft_properties::debug_print() const {
    std::cout << "--  spacecraft properties  --"  //
              << "\n  mass: " << mass_kg          //
              << '\n';

    std::visit([](auto&& s) { s.debug_print(); }, shape);

    hysteresis.debug_print();
    magnet.debug_print();

    for (const auto& hysteresis_rod : rods) {
        hysteresis_rod.debug_print();
    }
}

spacecraft::~spacecraft() = default;

auto spacecraft::create(const spacecraft_properties& properties) -> std::shared_ptr<spacecraft> {
    return std::make_shared<spacecraft_impl>(properties);
}

}  // namespace aos
