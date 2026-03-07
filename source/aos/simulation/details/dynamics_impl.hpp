#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/dynamics.hpp"

#include <memory>

namespace aos {

class dynamics_impl : public dynamics {
public:

    dynamics_impl(const dynamics_impl&)                    = delete;
    dynamics_impl(dynamics_impl&&)                         = delete;
    auto operator=(const dynamics_impl&) -> dynamics_impl& = delete;
    auto operator=(dynamics_impl&&) -> dynamics_impl&      = delete;

    dynamics_impl(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment> environment_model);
    ~dynamics_impl() override;

    void step(const system_state& current_state, system_state& state_derivative, real t_sec) const override;

    [[nodiscard]] auto get_spacecraft() const -> const spacecraft&;
    [[nodiscard]] auto get_environment() const -> const environment&;

private:

    std::shared_ptr<const spacecraft>  _spacecraft;
    std::shared_ptr<const environment> _environment;
    real                               _time_offset{};
};

}  // namespace aos
