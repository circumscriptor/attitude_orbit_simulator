#include "dynamics.hpp"

namespace aos {

dynamics::dynamics()  = default;
dynamics::~dynamics() = default;

auto dynamics::get_time_offset() const noexcept -> double {
    return _time_offset;
}

void dynamics::set_time_offset(double offset_s) {
    _time_offset = offset_s;
}

}  // namespace aos
