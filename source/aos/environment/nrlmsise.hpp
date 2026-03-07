#pragma once

#include "aos/core/types.hpp"
#include "aos/environment/space_weather.hpp"

#include <cstddef>
#include <filesystem>

extern "C" {
#include <nrlmsise-00.h>
}

namespace aos {

struct nrlmsise {
    mutable nrlmsise_input  input{};
    mutable nrlmsise_flags  flags{};
    mutable nrlmsise_output output{};
    space_weather           weather{};
    mutable size_t          hint{};

    explicit nrlmsise(const std::filesystem::path& filepath);

    // compute atmospheric density at a specific spacetime point
    [[nodiscard]] auto density_at(real_t year_decimal, real_t lat_deg, real_t lon_deg, real_t alt_m) const -> real_t;
};

}  // namespace aos
