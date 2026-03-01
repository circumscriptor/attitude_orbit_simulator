#pragma once

#include "aos/core/state.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>

namespace aos {

struct state_observer_properties {
    bool exclude_elements;    // per-element entries
    bool exclude_magnitudes;  // magnitude (vector length) entries

    void from_toml(const toml::table& table);
};

class state_observer {
public:

    explicit state_observer(const std::string& filename, std::size_t num_rods, const state_observer_properties& props);

    void operator()(const system_state& state, double time) const;

private:

    std::shared_ptr<std::ofstream> _file;
    std::size_t                    _num_rods;
    bool                           _include_elements;
    bool                           _include_magnitudes;
};

}  // namespace aos
