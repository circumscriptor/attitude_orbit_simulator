#pragma once

#include <Eigen/Eigen>

namespace aos {

using MagneticDipoleMoment = Eigen::Matrix<double, 3, 1>;
using GeomagneticField     = Eigen::Matrix<double, 3, 1>;
using InertiaTensor        = Eigen::Matrix<double, 3, 3>;

constexpr InertiaTensor get_ineratia_tensor(double m, double a, double b, double c) {
    const double i_x = (1. / 12.) * m * (b * b + c * c);
    const double i_y = (1. / 12.) * m * (a * a + c * c);
    const double i_z = (1. / 12.) * m * (a * a + b * b);

    InertiaTensor tensor = InertiaTensor::Identity();
    tensor(0, 0)         = i_x;
    tensor(1, 1)         = i_y;
    tensor(2, 2)         = i_z;
    return tensor;
}

}  // namespace aos
