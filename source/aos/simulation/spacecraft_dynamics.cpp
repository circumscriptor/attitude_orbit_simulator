#include "spacecraft_dynamics.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"

#include <memory>
#include <utility>

namespace aos {

spacecraft_dynamics::spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment> environment_model)
    : _spacecraft(std::move(spacecraft_model)), _environment(std::move(environment_model)) {}

spacecraft_dynamics::~spacecraft_dynamics() = default;

void spacecraft_dynamics::step(const system_state& current_state, system_state& state_derivative, double t_sec) const {
    const auto env = _environment->compute_effects(_time_offset + t_sec, current_state.position_m, current_state.velocity_m_s);
    _spacecraft->derivative(env, _environment->earth_mu(), current_state, state_derivative);
}

auto spacecraft_dynamics::get_spacecraft() const -> const spacecraft& {
    return *_spacecraft;
}

auto spacecraft_dynamics::get_environment() const -> const environment& {
    return *_environment;
}

}  // namespace aos
