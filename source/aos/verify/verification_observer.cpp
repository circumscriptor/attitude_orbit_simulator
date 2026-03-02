#include "verification_observer.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/observer.hpp"

#include <memory>
#include <ostream>
#include <print>
#include <string>
#include <utility>

namespace aos {

verification_observer::verification_observer(const std::string&                 filename,
                                             std::shared_ptr<const spacecraft>  sat,
                                             std::shared_ptr<const environment> env,
                                             const observer_properties&         props)
    : observer(filename, sat->rods().size(), props), _sat(std::move(sat)), _env(std::move(env)) {}

auto verification_observer::write_header() -> std::ostream& {
    return observer::write_header() << ",sun_x,sun_y,sun_z,shadow,mag_x,mag_y,mag_z,t_mag_x,t_mag_y,t_mag_z,t_grav_x,t_grav_y,t_grav_z";
}

auto verification_observer::write(const system_state& state, double time) -> std::ostream& {
    const auto env    = _env->compute_effects(time, state.position_m, state.velocity_m_s);
    const quat q_inv  = state.attitude.conjugate();
    const vec3 b_body = q_inv * env.magnetic_field_eci_T;
    const vec3 r_body = q_inv * state.position_m;
    const vec3 t_mag  = _sat->magnet().compute_torque(b_body);
    const vec3 t_grav = _sat->compute_gravity_gradient_torque(r_body, _env->earth_mu());

    auto& f = observer::write(state, time);
    std::print(f, ",{},{},{},{},{},{},{},{},{},{},{},{},{}",                                              //
               env.r_sun_eci.x(), env.r_sun_eci.y(), env.r_sun_eci.z(),                                   //
               env.shadow_factor,                                                                         //
               env.magnetic_field_eci_T.x(), env.magnetic_field_eci_T.y(), env.magnetic_field_eci_T.z(),  //
               t_mag.x(), t_mag.y(), t_mag.z(),                                                           //
               t_grav.x(), t_grav.y(), t_grav.z());
    return f;
}

}  // namespace aos
