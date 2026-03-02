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
    return observer::write_header() << ",sun_x,sun_y,sun_z,shadow,mag_x,mag_y,mag_z,mag_dot_x,mag_dot_y,mag_dot_z,"
                                       "t_mag_x,t_mag_y,t_mag_z,t_grav_x,t_grav_y,t_grav_z,t_gyro_x,t_gyro_y,t_gyro_z,"
                                       "t_rods_x,t_rods_y,t_rods_z,t_face_x,t_face_y,t_face_z,f_face_x,f_face_y,f_face_z,"
                                       "rho,solar_p";
}

auto verification_observer::write(const system_state& state, double time) -> std::ostream& {
    const auto env        = _env->compute_effects(time, state.position_m, state.velocity_m_s);
    const quat q_inv      = state.attitude.conjugate();
    const vec3 b_body     = q_inv * env.magnetic_field_eci_T;
    const vec3 r_body     = q_inv * state.position_m;
    const vec3 b_dot_orb  = q_inv * env.magnetic_field_dot_eci_T_s;
    const vec3 b_dot_body = b_dot_orb - state.angular_velocity_m_s.cross(b_body);
    const vec3 t_mag      = _sat->magnet().compute_torque(b_body);
    const vec3 t_grav     = _sat->compute_gravity_gradient_torque(r_body, _env->earth_mu());
    const vec3 t_gyro     = _sat->compute_gyroscopic_torque(state.angular_velocity_m_s);
    const vec3 t_rods     = _sat->compute_rod_torques(state.rod_magnetizations, b_body);
    const auto face_eff   = _sat->compute_face_effects(env, state.attitude, q_inv, state.angular_velocity_m_s);

    auto& f = observer::write(state, time);
    std::print(f,
               ",{},{},{},{},{},{},{},{},{},{},"  // Sun + Mag
               "{},{},{},{},{},{},{},{},{},"      // Mag + Grav + Gyro Torques
               "{},{},{},{},{},{},{},{},{},"      // Rods + Face Torques + Face Forces
               "{},{}",                           // Atmosphere + Solar Pressure
               env.r_sun_eci.x(), env.r_sun_eci.y(), env.r_sun_eci.z(), env.shadow_factor, env.magnetic_field_eci_T.x(), env.magnetic_field_eci_T.y(),
               env.magnetic_field_eci_T.z(), env.magnetic_field_dot_eci_T_s.x(), env.magnetic_field_dot_eci_T_s.y(), env.magnetic_field_dot_eci_T_s.z(),
               t_mag.x(), t_mag.y(), t_mag.z(), t_grav.x(), t_grav.y(), t_grav.z(), t_gyro.x(), t_gyro.y(), t_gyro.z(), t_rods.x(), t_rods.y(), t_rods.z(),
               face_eff.torque_body.x(), face_eff.torque_body.y(), face_eff.torque_body.z(), face_eff.force_eci.x(), face_eff.force_eci.y(),
               face_eff.force_eci.z(), env.atmospheric_density_kg_m3, env.solar_pressure_Pa);
    return f;
}

}  // namespace aos
