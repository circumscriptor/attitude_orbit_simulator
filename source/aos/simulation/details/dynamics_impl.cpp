#include "dynamics_impl.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <memory>
#include <utility>

namespace aos {

dynamics_impl::dynamics_impl(std::shared_ptr<const spacecraft> spacecraft_model, std::shared_ptr<const environment> environment_model)
    : _spacecraft(std::move(spacecraft_model)), _environment(std::move(environment_model)) {}

dynamics_impl::~dynamics_impl() = default;

void dynamics_impl::step(const system_state& current_state, system_state& state_derivative, real t_sec) const {
    const auto env = _environment->compute_effects(_time_offset + t_sec, current_state.position_m, current_state.velocity_m_s);
    _spacecraft->derivative(env, current_state, state_derivative);
}

auto dynamics_impl::get_spacecraft() const -> const spacecraft& {
    return *_spacecraft;
}

auto dynamics_impl::get_environment() const -> const environment& {
    return *_environment;
}

}  // namespace aos
