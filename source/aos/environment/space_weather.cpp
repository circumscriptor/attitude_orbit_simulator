#include "space_weather.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>

namespace aos {

auto space_weather::get_at(double year_decimal) const -> const space_weather_data& {
    const space_weather_data_vec* target = nullptr;

    if (year_decimal <= observed_cutoff_year_decimal) {
        target = &observed;
    } else if (year_decimal <= predicted_days_cutoff_year_decimal) {
        target = &predicted_days;
    } else if (year_decimal <= predicted_months_cutoff_year_decimal) {
        target = &predicted_months;
    } else {
        throw std::out_of_range("requested date is out of space weather data range");
    }

    if (target->empty()) {
        throw std::runtime_error("space weather data are empty");
    }

    // TODO: Optimize in the future to access by day or month without iterating
    auto it = std::ranges::lower_bound(*target, year_decimal, {}, &space_weather_data::year_decimal);
    if (it == target->end()) {
        return target->back();
    }
    return *it;
}

auto space_weather_parser::parse(const std::filesystem::path& filepath) -> space_weather {
    std::ifstream file(filepath);
    if (not file.is_open()) {
        throw std::runtime_error("could not open space weather file");
    }

    std::string line;
    if (not std::getline(file, line)) {
        throw std::runtime_error("space weather file is empty");
    }

    map_headers(line);
    if (not is_valid()) {
        throw std::runtime_error("space weather file has unexpected header");
    }

    space_weather result;
    while (std::getline(file, line)) {
        try {
            row(result, line);
        } catch (const std::exception& ex) {
            throw std::runtime_error(std::string("space weather file data are malformed: ") + ex.what());
        }
    }

    result.observed_cutoff_year_decimal =
        result.observed.empty() ? 0.0 : std::ranges::max(result.observed | std::views::transform(&space_weather_data::year_decimal));

    result.predicted_days_cutoff_year_decimal =
        result.predicted_days.empty() ? 0.0 : std::ranges::max(result.predicted_days | std::views::transform(&space_weather_data::year_decimal));

    result.predicted_months_cutoff_year_decimal =
        result.predicted_months.empty() ? 0.0 : std::ranges::max(result.predicted_months | std::views::transform(&space_weather_data::year_decimal));

    return result;
}

auto space_weather_parser::is_valid() const noexcept -> bool {
    return _column_date >= 0 && _column_f107 >= 0 && _column_f107a >= 0 && _column_type >= 0;
}

auto space_weather_parser::is_target_column(int index) const -> bool {
    return index == _column_date || index == _column_f107 || index == _column_f107a || index == _column_type;
}

void space_weather_parser::map_headers(std::string_view line) {
    _column_date = _column_f107 = _column_f107a = _column_type = -1;

    if (line.ends_with('\r')) {
        line.remove_suffix(1);
    }

    int column_counter = 0;
    for (auto part : std::views::split(line, ',')) {
        const std::string_view column_name(part.begin(), part.end());
        if (column_name == "DATE") {
            _column_date = column_counter;
        } else if (column_name == "F10.7_ADJ") {
            _column_f107 = column_counter;
        } else if (column_name == "F10.7_ADJ_LAST81") {
            _column_f107a = column_counter;
        } else if (column_name == "F10.7_DATA_TYPE") {
            _column_type = column_counter;
        }
        ++column_counter;
    }
}

void space_weather_parser::row(space_weather& result, std::string_view line) {
    space_weather_data data{};

    if (line.ends_with('\r')) {
        line.remove_suffix(1);
    }

    int column_counter = 0;
    for (auto part : std::views::split(line, ',')) {
        if (is_target_column(column_counter)) {
            std::string_view cell{part.begin(), part.end()};
            column(data, column_counter, cell);
        }
        ++column_counter;
    }

    store(result, data);
}

void space_weather_parser::column(space_weather_data& data, int index, std::string_view cell) const {
    if (index == _column_date) {
        data.year_decimal = convert_date(cell);
    } else if (index == _column_f107) {
        std::from_chars(cell.begin(), cell.end(), data.f107);
    } else if (index == _column_f107a) {
        std::from_chars(cell.begin(), cell.end(), data.f107a);
    } else if (index == _column_type) {
        data.type = convert_type(cell);
    }
}

auto space_weather_parser::convert_type(std::string_view type) -> space_weather_type {
    if (type == "OBS") {
        return space_weather_type_obs;
    }
    if (type == "INT") {
        return space_weather_type_int;
    }
    if (type == "PRD") {
        return space_weather_type_prd;
    }
    if (type == "PRM") {
        return space_weather_type_prm;
    }
    throw std::runtime_error("unknown data type");
}

void space_weather_parser::store(space_weather& result, const space_weather_data& data) {
    switch (data.type) {
        case space_weather_type_obs:
        case space_weather_type_int:
            result.observed.push_back(data);
            break;
        case space_weather_type_prd:
            result.predicted_days.push_back(data);
            break;
        case space_weather_type_prm:
            result.predicted_months.push_back(data);
            break;
        case space_weather_type_unknown:
            throw std::runtime_error("unknown or missing data type");
    }
}

auto space_weather_parser::convert_date(std::string_view date) -> double {
    static constexpr std::array days_before_month = std::to_array({0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334});

    // Format: YYYY-MM-DD
    int year{};
    int month{};
    int day{};

    {
        const auto date_year = date.substr(0, 4);
        std::from_chars(date_year.begin(), date_year.end(), year);
    }
    {
        const auto date_month = date.substr(5, 2);
        std::from_chars(date_month.begin(), date_month.end(), month);
    }
    {
        const auto date_day = date.substr(8, 2);
        std::from_chars(date_day.begin(), date_day.end(), day);
    }

    const bool   is_leap     = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    const int    day_of_year = days_before_month.at(month - 1) + day + static_cast<int>(is_leap && month > 2);
    const double total_days  = is_leap ? 366. : 365.;
    return static_cast<double>(year) + (static_cast<double>(day_of_year - 1) / total_days);
}

}  // namespace aos
