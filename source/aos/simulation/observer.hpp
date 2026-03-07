#pragma once

#include "aos/core/state.hpp"
#include "aos/core/types.hpp"

#include <cstddef>
#include <memory>
#include <ostream>
#include <string>

namespace aos {

struct observer_properties {
    static constexpr int default_precission = 5;

    bool exclude_elements{};              // per-element entries
    bool exclude_magnitudes{};            // magnitude (vector length) entries
    int  precission{default_precission};  // output number decimal precision

    void from_toml(const toml_table& table);
};

class observer {
public:

    observer()                                   = default;
    observer(const observer&)                    = delete;
    observer(observer&&)                         = delete;
    auto operator=(const observer&) -> observer& = delete;
    auto operator=(observer&&) -> observer&      = delete;

    virtual ~observer();

    virtual auto write_header() -> std::ostream&                              = 0;
    virtual auto write(const system_state& state, real time) -> std::ostream& = 0;

    static auto create(const std::string& filename, std::size_t num_rods, const observer_properties& properties) -> std::shared_ptr<observer>;
};

}  // namespace aos
