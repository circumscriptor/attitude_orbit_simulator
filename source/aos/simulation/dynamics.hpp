#pragma once

#include "aos/core/state.hpp"

namespace aos {

class dynamics {
public:

    dynamics(const dynamics&)                    = delete;
    dynamics(dynamics&&)                         = delete;
    auto operator=(const dynamics&) -> dynamics& = delete;
    auto operator=(dynamics&&) -> dynamics&      = delete;

    dynamics();
    virtual ~dynamics();

    virtual void step(const system_state& current_state, system_state& state_derivative, double t_sec) const = 0;

    [[nodiscard]]
    auto get_time_offset() const noexcept -> double;
    void set_time_offset(double offset_s);

private:

    double _time_offset{};
};

}  // namespace aos
