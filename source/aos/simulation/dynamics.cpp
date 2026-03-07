#include "dynamics.hpp"

#include "aos/core/types.hpp"
#include "aos/simulation/details/dynamics_impl.hpp"

#include <memory>
#include <utility>

namespace aos {

dynamics::dynamics()  = default;
dynamics::~dynamics() = default;

auto dynamics::get_time_offset() const noexcept -> real {
    return _time_offset;
}

void dynamics::set_time_offset(real offset_s) {
    _time_offset = offset_s;
}

auto dynamics::create(std::shared_ptr<const spacecraft> spacecraft, std::shared_ptr<const environment> environment) -> std::shared_ptr<dynamics> {
    return std::make_shared<dynamics_impl>(std::move(spacecraft), std::move(environment));
}

}  // namespace aos
