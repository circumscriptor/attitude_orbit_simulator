#include "aos/components/permanent_magnet.hpp"

#include "aos/core/constants.hpp"
#include "aos/core/types.hpp"

#include <gtest/gtest.h>

#include <numbers>

// Test fixture for permanent_magnet
using aos::vacuum_permeability;
using aos::vec3;
using aos::components::permanent_magnet;

class PermanentMagnetTest : public ::testing::Test {
protected:

    void SetUp() override {
        remanence   = 1.0;
        length_m    = 0.1;
        diameter_m  = 0.02;
        orientation = {1.0, 0.0, 2.0};
    }

    double remanence{};
    double length_m{};
    double diameter_m{};
    vec3   orientation;
};

// Test that constructor correctly initializes magnetic moment orientation
TEST_F(PermanentMagnetTest, MagneticMomentOrientation) {
    permanent_magnet magnet(remanence, length_m, diameter_m, orientation);
    vec3             moment = magnet.magnetic_moment();

    // Since magnetic moment depends on remanence, volume, and orientation,
    // check that returned moment direction matches orientation (unit vector)
    double norm_orientation = orientation.norm();
    ASSERT_GT(norm_orientation, 0.0);
    vec3 normalized_orientation = orientation / norm_orientation;

    double dot_product = moment.normalized().dot(normalized_orientation);
    EXPECT_NEAR(dot_product, 1.0, 1e-6) << "Magnetic moment direction should match orientation";
}

// Test that magnetic moment magnitude is positive and reasonable
TEST_F(PermanentMagnetTest, MagneticMomentMagnitude) {
    permanent_magnet magnet(remanence, length_m, diameter_m, orientation);
    vec3             moment = magnet.magnetic_moment();

    // Volume of cylindrical magnet
    double volume             = std::numbers::pi * (diameter_m * diameter_m / 4.) * length_m;
    double expected_magnitude = (remanence / vacuum_permeability) * volume;
    EXPECT_NEAR(moment.norm(), expected_magnitude, 1e-6) << "Magnetic moment magnitude incorrect";
}
