#pragma once

#include <Eigen/Eigen>
#include <expected>

namespace aos {

template <typename T>
struct is_eigen_matrix : std::false_type {};

template <typename S, int R, int C, int O, int MR, int MC>
struct is_eigen_matrix<Eigen::Matrix<S, R, C, O, MR, MC>> : std::true_type {};

template <typename T>
concept IsEigenMatrix = is_eigen_matrix<T>::value;

template <IsEigenMatrix StateType>
using Derivative = std::function<StateType(double, const StateType&)>;

template <IsEigenMatrix StateType>
using Step = std::function<void(double, const StateType&)>;

enum class rfk45_status { success, max_step_reached, step_size_underflow };

/**
 * @brief A production-ready C++ template for the RKF45 method using the Eigen library.
 *
 * This version is enhanced with diagnostic outputs, robustness checks, and greater user control.
 *
 * @tparam StateType An Eigen Vector type (e.g., Eigen::VectorXd).
 * @tparam UseLocalExtrapolation A compile-time boolean to select the 5th order solution for state updates.
 * @param f The function defining the ODE system, f(t, y).
 * @param y The initial state vector; will be updated to the final state.
 * @param t The initial time; will be updated to the final time.
 * @param t_end The end time of the integration.
 * @param h The initial step size; will be adapted.
 * @param tol The desired tolerance for the local error.
 * @param accepted_steps A counter for accepted steps.
 * @param rejected_steps A counter for rejected steps.
 * @param safety_factor The safety factor for step size adjustment (typically 0.8 to 0.9).
 * @param max_steps The maximum number of steps to prevent infinite loops.
 * @return An rfk45_status indicating the outcome of the integration.
 */
template <IsEigenMatrix StateType, bool UseLocalExtrapolation = false>
constexpr rfk45_status rkf45(const Derivative<StateType>& f,
                             StateType&                   y,
                             double&                      t,
                             double                       t_end,
                             double&                      h,
                             double                       tol,
                             long&                        accepted_steps,
                             long&                        rejected_steps,
                             double                       safety_factor = 0.9,
                             long                         max_steps     = 100000) noexcept {
    // Butcher tableau coefficients
    const double c2 = 1. / 4., a21 = 1. / 4.;
    const double c3 = 3. / 8., a31 = 3. / 32., a32 = 9. / 32.;
    const double c4 = 12. / 13., a41 = 1932. / 2197., a42 = -7200. / 2197., a43 = 7296. / 2197.;
    const double c5 = 1., a51 = 439. / 216., a52 = -8., a53 = 3680. / 513., a54 = -845. / 4104.;
    const double c6 = 1. / 2., a61 = -8. / 27., a62 = 2., a63 = -3544. / 2565., a64 = 1859. / 4104., a65 = -11. / 40.;
    const double b1_4th = 25. / 216., b3_4th = 1408. / 2565., b4_4th = 2197. / 4104., b5_4th = -1. / 5.;
    const double b1_5th = 16. / 135., b3_5th = 6656. / 12825., b4_5th = 28561. / 56430., b5_5th = -9. / 50., b6_5th = 2. / 55.;

    const double max_step_increase = 5.0;
    const double min_step_decrease = 0.1;

    while (t < t_end) {
        if (accepted_steps + rejected_steps >= max_steps) {
            return rfk45_status::max_step_reached;
        }

        if (t + h > t_end) {
            h = t_end - t;
        }

        const StateType k1 = f(t, y);
        const StateType k2 = f(t + c2 * h, y + h * a21 * k1);
        const StateType k3 = f(t + c3 * h, y + h * (a31 * k1 + a32 * k2));
        const StateType k4 = f(t + c4 * h, y + h * (a41 * k1 + a42 * k2 + a43 * k3));
        const StateType k5 = f(t + c5 * h, y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
        const StateType k6 = f(t + c6 * h, y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5));

        const StateType y4 = y + h * (b1_4th * k1 + b3_4th * k3 + b4_4th * k4 + b5_4th * k5);
        const StateType y5 = y + h * (b1_5th * k1 + b3_5th * k3 + b4_5th * k4 + b5_5th * k5 + b6_5th * k6);

        const double error = (y5 - y4).norm();

        if (error <= tol) {
            t += h;
            accepted_steps++;

            if constexpr (UseLocalExtrapolation) {
                y = y5;
            } else {
                y = y4;
            }

            if (error < std::numeric_limits<double>::epsilon()) {
                h *= max_step_increase;
            } else {
                double h_new = safety_factor * h * std::pow(tol / error, 0.2);
                h            = std::min(h * max_step_increase, h_new);
            }
        } else {
            rejected_steps++;
            double h_new = safety_factor * h * std::pow(tol / error, 0.25);  // Use 5th order formula for rejected steps
            h            = std::max(h * min_step_decrease, h_new);
        }

        if (std::abs(h) < std::numeric_limits<double>::epsilon() * t) {
            return rfk45_status::step_size_underflow;
        }

        if (t >= t_end)
            break;
    }
    return rfk45_status::success;
}

/**
 * @brief A simplified wrapper for the RKF45 solver with minimal inputs.
 *
 * This function handles the setup of diagnostic counters and provides sensible defaults
 * for advanced parameters, returning the final state upon success.
 *
 * @tparam StateType An Eigen Vector type (e.g., Eigen::VectorXd).
 * @tparam UseLocalExtrapolation A compile-time boolean to select the 5th order solution for state updates.
 * @param f The function defining the ODE system, f(t, y).
 * @param y0 The initial state vector at t0.
 * @param t0 The initial time.
 * @param t_end The end time of the integration.
 * @param h0 Optional: A suggested initial step size. If not provided, a default is calculated.
 * @param tol The desired tolerance for the local error.
 * @return A std::expected containing the final state vector if integration is successful,
 *         otherwise std::unexpected containing status.
 */
template <IsEigenMatrix StateType, bool UseLocalExtrapolation = false>
constexpr std::expected<StateType, rfk45_status> rkf45(const Derivative<StateType>& f, const StateType& y0, double t0, double t_end, double h0, double tol) {
    StateType y              = y0;
    double    t              = t0;
    long      accepted_steps = 0;
    long      rejected_steps = 0;
    // Determine a reasonable initial step size if not provided by the user
    // A common heuristic is a fraction of the total integration interval.
    double       h             = (h0 > 0.0) ? h0 : (t_end - t0) / 100.0;
    const double safety_factor = 0.9;
    const long   max_steps     = (t_end - t0) * 1000;

    rfk45_status status = rkf45<StateType, UseLocalExtrapolation>(f, y, t, t_end, h, tol, accepted_steps, rejected_steps, safety_factor, max_steps);
    if (status == rfk45_status::success) {
        return y;
    }
    return std::unexpected(status);
}

}  // namespace aos
