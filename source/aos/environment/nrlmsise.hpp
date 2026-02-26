#pragma once

extern "C" {
#include <nrlmsise-00.h>
}

namespace aos {

struct nrlmsise {
    nrlmsise_input  input{};
    nrlmsise_flags  flags{};
    nrlmsise_output output{};

    explicit nrlmsise();

    /** Compute atmospheric density at a specific spacetime point */
    [[nodiscard]] auto density_at(double year_decimal, double lat_deg, double lon_deg, double alt_m) -> double;
};

}  // namespace aos
