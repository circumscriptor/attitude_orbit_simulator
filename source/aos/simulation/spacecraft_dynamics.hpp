#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/dynamics.hpp"

#include <memory>

namespace aos {

class spacecraft_dynamics : public dynamics {
public:

    spacecraft_dynamics(const spacecraft_dynamics&)                    = delete;
    spacecraft_dynamics(spacecraft_dynamics&&)                         = delete;
    auto operator=(const spacecraft_dynamics&) -> spacecraft_dynamics& = delete;
    auto operator=(spacecraft_dynamics&&) -> spacecraft_dynamics&      = delete;

    spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment> environment_model);
    ~spacecraft_dynamics() override;

    void step(const system_state& current_state, system_state& state_derivative, double t_sec) const override;

    [[nodiscard]] auto get_spacecraft() const -> const spacecraft&;
    [[nodiscard]] auto get_environment() const -> const environment&;

private:

    std::shared_ptr<const spacecraft>  _spacecraft;
    std::shared_ptr<const environment> _environment;
    double                             _time_offset{};
};

}  // namespace aos
