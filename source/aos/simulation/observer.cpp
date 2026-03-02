#include "observer.hpp"

#include "aos/simulation/details/observer_impl.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
#include <cstddef>
#include <memory>
#include <string>
// clang-format on

namespace aos {

void observer_properties::from_toml(const toml::table& table) {
    exclude_elements   = table["exclude_elements"].value_or(false);
    exclude_magnitudes = table["exclude_magnitudes"].value_or(false);
    precission         = table["precission"].value_or(default_precission);
}

observer::~observer() = default;

auto observer::create(const std::string& filename, std::size_t num_rods, const observer_properties& properties) -> std::shared_ptr<observer> {
    return std::make_shared<observer_impl>(filename, num_rods, properties);
}

}  // namespace aos
