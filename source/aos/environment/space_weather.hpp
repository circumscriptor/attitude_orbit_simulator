#pragma once

#include "aos/core/types.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <vector>

namespace aos {

enum space_weather_type : uint8_t {
    space_weather_type_obs,
    space_weather_type_int,
    space_weather_type_prd,
    space_weather_type_prm,
    space_weather_type_unknown,
};

struct space_weather_data {
    real               year_decimal;
    real               f107;
    real               f107a;
    space_weather_type type;
};

using space_weather_data_vec = std::vector<space_weather_data>;

struct space_weather {
    space_weather_data_vec observed;
    space_weather_data_vec predicted_days;
    space_weather_data_vec predicted_months;

    real observed_cutoff_year_decimal;
    real predicted_days_cutoff_year_decimal;
    real predicted_months_cutoff_year_decimal;

    [[nodiscard]] auto get_at(real year_decimal, size_t& hint) const -> const space_weather_data&;
    [[nodiscard]] auto get_month_at(real year_decimal) const -> const space_weather_data&;
    [[nodiscard]] auto get_linear_month_at(real year_decimal) const -> space_weather_data;
};

class space_weather_parser {
public:

    auto parse(const std::filesystem::path& filepath) -> space_weather;

protected:

    [[nodiscard]] auto is_valid() const noexcept -> bool;
    [[nodiscard]] auto is_target_column(int index) const -> bool;

    void map_headers(std::string_view header_line);

    void row(space_weather& result, std::string_view line);
    void column(space_weather_data& data, int index, std::string_view cell) const;

    static void store(space_weather& result, const space_weather_data& data);

    static auto convert_type(std::string_view type) -> space_weather_type;
    static auto convert_date(std::string_view date) -> real;

private:

    int _column_date{-1};
    int _column_f107{-1};
    int _column_f107a{-1};
    int _column_type{-1};
};

}  // namespace aos
