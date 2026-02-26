#pragma once

#include "aos/environment/space_weather.hpp"

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

    explicit nrlmsise(const std::filesystem::path& filepath);

    /** Compute atmospheric density at a specific spacetime point */
    [[nodiscard]] auto density_at(double year_decimal, double lat_deg, double lon_deg, double alt_m) const -> double;
};

}  // namespace aos
