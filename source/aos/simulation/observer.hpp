#pragma once

#include "aos/core/state.hpp"

// clang-format off
#include <toml++/toml.hpp>
#include <toml++/impl/table.hpp>
// clang-format on

#include <cstddef>
#include <fstream>
#include <ostream>
#include <string>

namespace aos {

struct observer_properties {
    static constexpr int default_precission = 5;

    bool exclude_elements{};              // per-element entries
    bool exclude_magnitudes{};            // magnitude (vector length) entries
    int  precission{default_precission};  // output number decimal precision

    void from_toml(const toml::table& table);
};

class observer {
public:

    observer(const observer&)                    = delete;
    observer(observer&&)                         = delete;
    auto operator=(const observer&) -> observer& = delete;
    auto operator=(observer&&) -> observer&      = delete;

    explicit observer(const std::string& filename, std::size_t num_rods, const observer_properties& properties);
    virtual ~observer();

    virtual auto write_header() -> std::ostream&;
    virtual auto write(const system_state& state, double time) -> std::ostream&;

protected:

    auto file() -> std::ofstream&;

private:

    std::ofstream _file;
    std::size_t   _num_rods;
    bool          _include_elements;
    bool          _include_magnitudes;
};

}  // namespace aos
