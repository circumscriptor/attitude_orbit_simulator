#include "dynamics.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"

#include <memory>
#include <utility>

namespace aos {

spacecraft_dynamics::spacecraft_dynamics(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment_model> environment_model)
    : _spacecraft(std::move(spacecraft_model)), _environment(std::move(environment_model)) {}

void spacecraft_dynamics::operator()(const system_state& current_state, system_state& state_derivative, double t_sec) const {
    const auto env_data = _environment->calculate(_global_time_offset + t_sec, current_state.position_m, current_state.velocity_m_s);
    _spacecraft->derivative(env_data, _environment->earth_mu(), current_state, state_derivative);
}

void spacecraft_dynamics::set_global_time_offset(double offset_s) {
    _global_time_offset = offset_s;
}

}  // namespace aos
