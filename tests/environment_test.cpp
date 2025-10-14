#include "aos/environment/environment.hpp"

#include "aos/core/types.hpp"
#include "gtest/gtest.h"

#include <cmath>
#include <numbers>

using aos::vec3;
using aos::wmm2025_environment;

// Test fixture for the wmm2025_environment class
class Wmm2020EnvironmentTest : public ::testing::Test {
protected:

    void SetUp() override {
        props                       = {};
        props.orbit_altitude_km     = 500.0;
        props.orbit_inclination_deg = 45.0;
    }

    wmm2025_environment::properties props{};
};

TEST_F(Wmm2020EnvironmentTest, ConstructorInitializesProperties) {
    // This test primarily ensures that the constructor does not throw
    // an exception and can be instantiated.
    wmm2025_environment env(props);
    SUCCEED();
}

TEST_F(Wmm2020EnvironmentTest, EquatorialOrbitAtTimeZero) {
    props.orbit_inclination_deg = 0.0;
    wmm2025_environment env(props);

    vec3 b_inertial = env.inertial_magnetic_field_at(0.0);

    // Ground Truth for 2025.0 at lat=0, lon=0, alt=500km, rotated to ECI.
    // At t=0, the ECEF and ECI frames are aligned.
    // These values are computed from WMM2020 model via GeographicLib
    vec3 expected_b_inertial = {1.0906478300890021e-05, -2.1590470e-06, 2.1287368e-05};

    ASSERT_NEAR(b_inertial.x(), expected_b_inertial.x(), 1e-6);
    ASSERT_NEAR(b_inertial.y(), expected_b_inertial.y(), 1e-6);
    ASSERT_NEAR(b_inertial.z(), expected_b_inertial.z(), 1e-6);
}

TEST_F(Wmm2020EnvironmentTest, PolarOrbitAtTimeZero) {
    props.orbit_inclination_deg = 90.0;
    wmm2025_environment env(props);

    // At t=0, the satellite is at (lat=0, lon=0), so the result should be identical.
    vec3 b_inertial_t0          = env.inertial_magnetic_field_at(0.0);
    vec3 expected_b_inertial_t0 = {1.0906478300890021e-05, -2.1590470e-06, 2.1287368e-05};

    ASSERT_NEAR(b_inertial_t0.x(), expected_b_inertial_t0.x(), 1e-6);
    ASSERT_NEAR(b_inertial_t0.y(), expected_b_inertial_t0.y(), 1e-6);
    ASSERT_NEAR(b_inertial_t0.z(), expected_b_inertial_t0.z(), 1e-6);

    // Check that the field is different at a later time (over the pole)
    const double earth_radius_m = 6378137.0;
    const double alt_m          = props.orbit_altitude_km * 1000.0;
    const double orbit_radius_m = earth_radius_m + alt_m;
    const double orbit_period_s = 2.0 * std::numbers::pi * std::sqrt(std::pow(orbit_radius_m, 3) / 3.986004418e14);

    vec3 b_inertial_t_quarter = env.inertial_magnetic_field_at(orbit_period_s / 4.0);
    EXPECT_GT((b_inertial_t0 - b_inertial_t_quarter).norm(), 1e-6);
}

TEST_F(Wmm2020EnvironmentTest, FieldVariesWithTime) {
    wmm2025_environment env(props);

    vec3 b_t0   = env.inertial_magnetic_field_at(0.0);
    vec3 b_t100 = env.inertial_magnetic_field_at(100.0);

    // The field vector should change due to satellite motion and Earth rotation
    EXPECT_GT((b_t0 - b_t100).norm(), 1e-7);

    // Magnitudes should remain within a plausible range for LEO (20-60 microtesla)
    EXPECT_GT(b_t0.norm(), 2e-5);
    EXPECT_LT(b_t0.norm(), 6e-5);
    EXPECT_GT(b_t100.norm(), 2e-5);
    EXPECT_LT(b_t100.norm(), 6e-5);
}
