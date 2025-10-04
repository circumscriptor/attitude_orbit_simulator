#pragma once

#include <Eigen/Eigen>

namespace aos {

using Vector3d = Eigen::Matrix<double, 3, 1>;
using Packed6d = Eigen::Matrix<double, 6, 1>;  // packed two vectors 3d

static constexpr double standard_gravitational_parameter_earth = 3.986004418e14;

/**
 * @brief Calculates the time derivatives for a two-body orbit around Earth.
 *
 * This function implements the differential equations for the two-body problem.
 * The state vector is composed of the position and velocity vectors. The derivative
 * of the state is composed of the velocity and acceleration vectors.
 *
 * d/dt [r] = [v]
 *      [v]   [a]
 *
 * where a = -μ * r / |r|^3
 *
 * @param t Current time (unused).
 * @param state The 6-element state vector [rx, ry, rz, vx, vy, vz] in meters and m/s.
 * @return The 6-element derivative vector [vx, vy, vz, ax, ay, az].
 */
constexpr Packed6d orbit_derivatives(double t, const Packed6d& state) {
    (void)t;

    // 1. Unpack the state vector into position and velocity
    Vector3d position = state.head<3>();
    Vector3d velocity = state.tail<3>();

    // 2. Calculate the magnitude of the position vector (distance from center of Earth)
    double r_norm = position.norm();

    // Safety check to prevent division by zero if the object is at the singularity
    if (r_norm < 1e-6)
        throw std::runtime_error("Position vector magnitude is near zero, cannot compute gravitational acceleration.");

    // 3. Calculate the gravitational acceleration vector (a = dv/dt)
    double   r_norm_cubed = r_norm * r_norm * r_norm;
    Vector3d acceleration = -(standard_gravitational_parameter_earth / r_norm_cubed) * position;

    // 4. Pack the derivatives into a new 6D vector
    // The derivative of position is velocity (dr/dt = v)
    // The derivative of velocity is acceleration (dv/dt = a)
    Packed6d derivatives;
    derivatives.head<3>() = velocity;
    derivatives.tail<3>() = acceleration;
    return derivatives;
}

}  // namespace aos
