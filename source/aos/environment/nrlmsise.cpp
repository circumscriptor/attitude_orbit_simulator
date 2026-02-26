#include "nrlmsise.hpp"

#include "aos/core/constants.hpp"
#include "aos/environment/space_weather.hpp"

#include <cmath>
#include <filesystem>

namespace aos {

nrlmsise::nrlmsise(const std::filesystem::path& filepath) : weather(space_weather_parser{}.parse(filepath)) {
    for (int& fswitch : flags.switches) {
        fswitch = 1;
    }
    // flags.switches[0] = 0; // Uncomment to use g/cm3 instead of kg/m3
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto nrlmsise::density_at(double year_decimal, double lat_deg, double lon_deg, double alt_m) const -> double {
    const int    year                = static_cast<int>(std::floor(year_decimal));
    const double day_of_year_decimal = (year_decimal - year) * year_to_days;
    const int    day_of_year         = static_cast<int>(std::floor(day_of_year_decimal));
    const double seconds             = (day_of_year_decimal - day_of_year) * day_to_seconds;

    const double alt_km           = alt_m * 0.001;
    const double local_solar_time = (seconds / 3600.0) + (lon_deg / 15.0);

    const auto& data = weather.get_at(year_decimal);

    input.alt    = alt_km;
    input.g_lat  = lat_deg;
    input.g_long = lon_deg;
    input.doy    = day_of_year;
    input.sec    = seconds;
    input.lst    = local_solar_time;
    input.f107A  = data.f107a;
    input.f107   = data.f107;
    input.ap     = 0.0;      // TODO: what to do with this?
    input.ap_a   = nullptr;  // TODO: what to do with this?
    gtd7d(&input, &flags, &output);
    return output.d[5];  // NOLINT(readability-magic-numbers)
}

}  // namespace aos
