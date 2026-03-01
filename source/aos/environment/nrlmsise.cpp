#include "nrlmsise.hpp"

#include "aos/core/constants.hpp"
#include "aos/environment/space_weather.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>

namespace aos {

nrlmsise::nrlmsise(const std::filesystem::path& filepath) : weather(space_weather_parser{}.parse(filepath)) {
    // NOLINTBEGIN
    for (size_t i = 0; i < 24; ++i) {
        flags.switches[i] = 1;
        flags.sw[i]       = 0.0;
        flags.swc[i]      = 0.0;
    }
    // flags.switches[9] = 0;
    // flags.switches[0] = 0; // Uncomment to use g/cm3 instead of kg/m3
    // NOLINTEND
    for (double& d : output.d) {
        d = 0.0;
    }
    for (double& t : output.t) {
        t = 0.0;
    }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto nrlmsise::density_at(double year_decimal, double lat_deg, double lon_deg, double alt_m) const -> double {
    const int    year                  = static_cast<int>(std::floor(year_decimal));
    const bool   is_leap               = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    const double days_in_year          = is_leap ? 366.0 : 365.0;
    const double total_seconds_in_year = (year_decimal - year) * days_in_year * 86400.0;
    const auto&  data                  = weather.get_linear_month_at(year_decimal);

    // NOLINTBEGIN(readability-magic-numbers)
    input.doy    = static_cast<int>(std::floor(total_seconds_in_year / 86400.0)) + 1;
    input.sec    = std::fmod(total_seconds_in_year, 86400.0);
    input.alt    = alt_m * meter_to_kilometer;
    input.g_lat  = lat_deg;
    input.g_long = lon_deg;
    input.g_lat  = lat_deg;
    input.g_long = lon_deg;
    input.lst    = (input.sec / 3600.0) + (lon_deg / 15.0);

    if (input.alt < 80.0) {
        input.f107A = 150.0;
        input.f107  = 150.0;
        input.ap    = 4.0;
    } else {
        input.f107A = data.f107a;
        input.f107  = data.f107;
        input.ap    = 4.0;  // Default quiet-time magnetic index
    }
    // NOLINTEND(readability-magic-numbers)

    input.ap_a = nullptr;
    gtd7d(&input, &flags, &output);
    return output.d[5];  // NOLINT(readability-magic-numbers)
}

}  // namespace aos
