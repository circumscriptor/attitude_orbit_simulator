#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/dynamics.hpp"

#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>

using aos::environment;
using aos::mat3x3;
using aos::permanent_magnet;
using aos::quat;
using aos::spacecraft;
using aos::spacecraft_dynamics;
using aos::system_state;
using aos::vec3;
using aos::wmm2025_environment;

// --- Mock Environments for Isolated Testing ---

class ZeroFieldEnvironment : public environment {
public:

    ZeroFieldEnvironment() = default;

    [[nodiscard]]
    vec3 inertial_magnetic_field_at([[maybe_unused]] double t_sec) const override {
        return {0, 0, 0};
    }
};

class ConstantFieldEnvironment : public environment {
public:

    ConstantFieldEnvironment() = default;

    [[nodiscard]]
    vec3 inertial_magnetic_field_at([[maybe_unused]] double t_sec) const override {
        return {0.0, 3e-5, 0.0};  // B field along Y-axis
    }
};

// --- Test Cases ---

TEST(SimulationPermanentMagnetTest, CorrectMomentCalculation) {
    // Grade N52 has Br = 1.45 T. M = Br / mu0. Moment = M * Volume.
    // Volume = pi * (0.005^2) * 0.05 = 3.927e-6 m^3
    // M = 1.45 / (1.256637e-6) = 1.154e6 A/m
    // Moment Magnitude = 1.154e6 * 3.927e-6 = 4.53 A*m^2
    permanent_magnet magnet(1.45, 0.05, 0.01, {0, 0, 1});
    EXPECT_NEAR(magnet.magnetic_moment().norm(), 4.53, 1e-2);
}

TEST(SpacecraftTest, InertiaTensorCalculation) {
    // For a 12kg, 2x2x2m cube
    double m    = 12.0;
    double side = 2.0;
    // I = (1/12) * m * (side^2 + side^2) = (1/12)*12*(4+4) = 8
    mat3x3 inertia = spacecraft::get_inertia_tensor(m, side, side, side);
    EXPECT_NEAR(inertia(0, 0), 8.0, 1e-9);
    EXPECT_NEAR(inertia(1, 1), 8.0, 1e-9);
    EXPECT_NEAR(inertia(2, 2), 8.0, 1e-9);
}

TEST(EnvironmentTest, MagneticFieldPlausibility) {
    wmm2025_environment::properties props{};
    props.orbit_altitude_km     = 500.0;
    props.orbit_inclination_deg = 51.6;
    wmm2025_environment env(props);  // Use the base environment class

    vec3 b0   = env.inertial_magnetic_field_at(0.0);
    vec3 b100 = env.inertial_magnetic_field_at(100.0);

    ASSERT_GT(b0.norm(), 1e-9);
    ASSERT_NE(b0, b100);
    ASSERT_GE(b0.norm(), 2.0e-5);
    ASSERT_LE(b0.norm(), 6.0e-5);
}

// --- Test Fixture for Dynamics ---
class DynamicsTest : public ::testing::Test {
protected:

    // This method is called before each test
    void SetUp() override {
        // Setup a simple diagonal inertia tensor for easy calculations
        inertia_ << 0.1, 0, 0, 0, 0.2, 0, 0, 0, 0.3;

        spacecraft::properties props;
        props.magnet_remanence   = 1.45;
        props.magnet_length      = 0.05;
        props.magnet_diameter    = 0.01;
        props.magnet_orientation = {0, 0, 1};  // Along Z
        // FIX: The properties struct must be fully initialized
        props.hysteresis_rod_orientations.emplace_back(1, 0, 0);  // Add at least one rod

        satellite_ = std::make_shared<spacecraft>(inertia_, props);

        // FIX: Initialize the state objects here to be used in all tests
        current_state.attitude = quat::Identity();
        current_state.rod_magnetizations.resize(static_cast<int>(satellite_->rods().size()));
        current_state.rod_magnetizations.setZero();

        derivative.rod_magnetizations.resize(static_cast<int>(satellite_->rods().size()));
    }

    mat3x3                      inertia_;
    std::shared_ptr<spacecraft> satellite_;
    system_state                current_state;
    system_state                derivative;
};

TEST_F(DynamicsTest, EquilibriumProducesNoAcceleration) {
    current_state.angular_velocity = {0, 0, 0};  // No rotation

    auto zero_env = std::make_shared<ZeroFieldEnvironment>();
    // REFINEMENT: Pass by reference, not shared_ptr
    spacecraft_dynamics zero_dynamics(satellite_, zero_env);

    zero_dynamics(current_state, derivative, 0.0);

    EXPECT_NEAR(derivative.angular_velocity.norm(), 0.0, 1e-12);
}

TEST_F(DynamicsTest, GyroscopicTorqueIsCorrect) {
    // Let's pick an interesting spin
    current_state.angular_velocity = {0.1, 0.5, 0.3};

    auto                zero_env = std::make_shared<ZeroFieldEnvironment>();
    spacecraft_dynamics zero_dynamics(satellite_, zero_env);
    zero_dynamics(current_state, derivative, 0.0);

    // Expected torque: T = -ω x (Iω)
    const vec3 w               = current_state.angular_velocity;
    const vec3 inertia_w       = inertia_ * w;
    const vec3 expected_torque = -w.cross(inertia_w);
    const vec3 expected_accel  = inertia_.inverse() * expected_torque;

    EXPECT_NEAR(derivative.angular_velocity.x(), expected_accel.x(), 1e-9);
    EXPECT_NEAR(derivative.angular_velocity.y(), expected_accel.y(), 1e-9);
    EXPECT_NEAR(derivative.angular_velocity.z(), expected_accel.z(), 1e-9);
}

TEST_F(DynamicsTest, MagneticTorqueIsCorrect) {
    current_state.angular_velocity = {0, 0, 0};  // No spin

    auto                const_env = std::make_shared<ConstantFieldEnvironment>();
    spacecraft_dynamics const_dynamics(satellite_, const_env);
    const_dynamics(current_state, derivative, 0.0);

    // Expected torque: T = m x B
    // m is along Z, B is along Y. T should be along -X.
    vec3 m               = satellite_->magnet().magnetic_moment();
    vec3 b               = {0.0, 3e-5, 0.0};
    vec3 expected_torque = m.cross(b);
    vec3 expected_accel  = inertia_.inverse() * expected_torque;

    EXPECT_NEAR(derivative.angular_velocity.x(), expected_accel.x(), 1e-9);
    EXPECT_NEAR(derivative.angular_velocity.y(), expected_accel.y(), 1e-9);
    EXPECT_NEAR(derivative.angular_velocity.z(), expected_accel.z(), 1e-9);
    EXPECT_GT(std::abs(derivative.angular_velocity.x()), 1e-9);
    EXPECT_NEAR(derivative.angular_velocity.y(), 0.0, 1e-9);
}
