#include "nrlmsise.hpp"

#include "aos/core/constants.hpp"

#include <cmath>

namespace aos {

nrlmsise::nrlmsise() {
    for (int& fswitch : flags.switches) {
        fswitch = 1;
    }
    // flags.switches[0] = 0;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto nrlmsise::density_at(double year_decimal, double lat_deg, double lon_deg, double alt_m) -> double {
    const int    year                = static_cast<int>(std::floor(year_decimal));
    const double day_of_year_decimal = (year_decimal - year) * year_to_days;
    const int    day_of_year         = static_cast<int>(std::floor(day_of_year_decimal));
    const double seconds             = (day_of_year_decimal - day_of_year) * day_to_seconds;

    const double alt_km           = alt_m * 0.001;
    const double local_solar_time = (seconds / 3600.0) + (lon_deg / 15.0);

    input.alt    = alt_km;
    input.g_lat  = lat_deg;
    input.g_long = lon_deg;
    input.doy    = day_of_year;
    input.sec    = seconds;
    input.lst    = local_solar_time;
    input.f107A  = 80.0;     // TODO: use real values
    input.f107   = 150.0;    // TODO: use real values
    input.ap     = 0.0;      // TODO: what to do with this?
    input.ap_a   = nullptr;  // TODO: what to do with this?
    gtd7d(&input, &flags, &output);
    return output.d[5];  // NOLINT(readability-magic-numbers)
}

}  // namespace aos
