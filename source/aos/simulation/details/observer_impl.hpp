#pragma once

#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/simulation/observer.hpp"

#include <cstddef>
#include <fstream>
#include <ostream>
#include <string>

namespace aos {

class observer_impl : public observer {
public:

    observer_impl(const observer_impl&)                    = delete;
    observer_impl(observer_impl&&)                         = delete;
    auto operator=(const observer_impl&) -> observer_impl& = delete;
    auto operator=(observer_impl&&) -> observer_impl&      = delete;

    explicit observer_impl(const std::string& filename, std::size_t num_rods, const observer_properties& properties);
    ~observer_impl() override;

    auto write_header() -> std::ostream& override;
    auto write(const system_state& state, real time) -> std::ostream& override;

private:

    std::ofstream _file;
    std::size_t   _num_rods;
    bool          _include_elements;
    bool          _include_magnitudes;
};

}  // namespace aos
