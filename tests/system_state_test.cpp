#include "aos/core/state.hpp"
#include "aos/core/types.hpp"

#include <gtest/gtest.h>

using aos::core::quat;
using aos::core::system_state;

// Test fixture to set up a reusable state object
class SystemStateTest : public ::testing::Test {
protected:

    void SetUp() override {
        s1.attitude         = quat(1.0, 0.0, 0.0, 0.0);
        s1.angular_velocity = {1.0, -2.0, 3.0};
        s1.rod_magnetizations.resize(2);
        s1.rod_magnetizations << 10.0, -20.0;

        s2.attitude         = quat(0.0, 1.0, 0.0, 0.0);
        s2.angular_velocity = {0.5, 0.5, -0.5};
        s2.rod_magnetizations.resize(2);
        s2.rod_magnetizations << 5.0, 5.0;
    }

    system_state s1;
    system_state s2;
};

TEST_F(SystemStateTest, Addition) {
    system_state result = s1 + s2;
    EXPECT_NEAR(result.attitude.w(), 1.0, 1e-9);
    EXPECT_NEAR(result.attitude.x(), 1.0, 1e-9);
    EXPECT_NEAR(result.angular_velocity.x(), 1.5, 1e-9);
    EXPECT_NEAR(result.rod_magnetizations(1), -15.0, 1e-9);
}

TEST_F(SystemStateTest, ScalarMultiplication) {
    system_state result = s1 * 2.5;
    EXPECT_NEAR(result.attitude.w(), 2.5, 1e-9);
    EXPECT_NEAR(result.angular_velocity.y(), -5.0, 1e-9);
    EXPECT_NEAR(result.rod_magnetizations(0), 25.0, 1e-9);
}

TEST_F(SystemStateTest, CommutativeScalarMultiplication) {
    system_state result = 2.5 * s1;
    EXPECT_NEAR(result.attitude.w(), 2.5, 1e-9);
    EXPECT_NEAR(result.angular_velocity.y(), -5.0, 1e-9);
}

TEST_F(SystemStateTest, ScalarAddition) {
    system_state result = 10.0 + s1;
    EXPECT_NEAR(result.attitude.w(), 11.0, 1e-9);
    EXPECT_NEAR(result.angular_velocity.y(), 8.0, 1e-9);
    EXPECT_NEAR(result.rod_magnetizations(1), -10.0, 1e-9);
}

TEST_F(SystemStateTest, Abs) {
    system_state result = abs(s1);
    EXPECT_GE(result.attitude.w(), 0.0);
    EXPECT_NEAR(result.angular_velocity.x(), 1.0, 1e-9);
    EXPECT_NEAR(result.angular_velocity.y(), 2.0, 1e-9);
    EXPECT_NEAR(result.rod_magnetizations(1), 20.0, 1e-9);
}
