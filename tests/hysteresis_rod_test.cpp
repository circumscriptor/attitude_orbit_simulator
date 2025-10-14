#include "aos/components/hysteresis_rod.hpp"

#include "aos/core/types.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>

using aos::hysteresis_rod;
using aos::vec3;

class HysteresisRodTest : public ::testing::Test {
protected:

    void SetUp() override {
        hysteresis_rod::ja_parameters params{.ms = 1.4e5, .a = 2.0e3, .k = 1.0e3, .c = 0.2, .alpha = 1.0e-3};
        rod = std::make_shared<hysteresis_rod>(0.1, vec3(1, 0, 0), params);
    }

    std::shared_ptr<hysteresis_rod> rod;
};

TEST_F(HysteresisRodTest, ConstructorAndMagneticMoment) {
    const vec3 moment = rod->magnetic_moment(1.0e4);
    EXPECT_NEAR(moment.x(), 1.0e3, 1e-9);
    EXPECT_NEAR(moment.y(), 0.0, 1e-9);
    EXPECT_NEAR(moment.z(), 0.0, 1e-9);
}

TEST_F(HysteresisRodTest, MagnetizationDerivativeFromH) {
    const double m_scalar_am = 5.0e4;
    const double h_along_rod = 1.5e3;
    const double dh_dt       = 1.0e2;
    const double dm_dt       = rod->magnetization_derivative_from_h(m_scalar_am, h_along_rod, dh_dt);
    EXPECT_NEAR(dm_dt, -782.51825848, 1e-6);
}

TEST_F(HysteresisRodTest, MagnetizationDerivative) {
    double       m_scalar_am = 5.0e4;
    const vec3   b_body_t(0.002, 0.001, 0.0);
    const vec3   omega_rad_s(0.0, 0.0, 0.1);
    const double dm_dt = rod->magnetization_derivative(m_scalar_am, b_body_t, omega_rad_s);
    EXPECT_NEAR(dm_dt, -510.24644226, 1e-6);
}

TEST_F(HysteresisRodTest, NearZeroEffectiveField) {
    const double m_scalar_am = 0.0;
    const double h_along_rod = 0.0;
    const double dh_dt       = 1.0e-7;
    const double dm_dt       = rod->magnetization_derivative_from_h(m_scalar_am, h_along_rod, dh_dt);
    EXPECT_FALSE(std::isnan(dm_dt));
    EXPECT_FALSE(std::isinf(dm_dt));
}

TEST_F(HysteresisRodTest, NearZeroDenominator) {
    // These values are specifically chosen to make the denominator
    // in the dmirr_dh calculation close to zero.
    const double m_scalar_am = 1.3e5;
    const double h_along_rod = 2.5e3;
    const double dh_dt       = 1.0;
    const double dm_dt       = rod->magnetization_derivative_from_h(m_scalar_am, h_along_rod, dh_dt);
    EXPECT_FALSE(std::isnan(dm_dt));
    EXPECT_FALSE(std::isinf(dm_dt));
}

TEST_F(HysteresisRodTest, NegativeDhDt) {
    const double m_scalar_am = 5.0e4;
    const double h_along_rod = 1.5e3;
    const double dh_dt       = -1.0e2;  // Negative rate of change
    const double dm_dt       = rod->magnetization_derivative_from_h(m_scalar_am, h_along_rod, dh_dt);
    EXPECT_NEAR(dm_dt, -1650.58156137, 1e-6);
}
