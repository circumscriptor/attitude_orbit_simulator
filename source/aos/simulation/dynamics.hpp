#pragma once

#include "aos/core/state.hpp"
#include "aos/core/types.hpp"

#include <memory>

namespace aos {

class spacecraft;
class environment;

class dynamics {
public:

    dynamics(const dynamics&)                    = delete;
    dynamics(dynamics&&)                         = delete;
    auto operator=(const dynamics&) -> dynamics& = delete;
    auto operator=(dynamics&&) -> dynamics&      = delete;

    dynamics();
    virtual ~dynamics();

    virtual void step(const system_state& current_state, system_state& state_derivative, real t_sec) const = 0;

    [[nodiscard]]
    auto get_time_offset() const noexcept -> real;
    void set_time_offset(real offset_s);

    static auto create(std::shared_ptr<const spacecraft> spacecraft, std::shared_ptr<const environment> environment) -> std::shared_ptr<dynamics>;

private:

    real _time_offset{};
};

}  // namespace aos
