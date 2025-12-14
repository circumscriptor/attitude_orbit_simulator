#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace aos {

// NOLINTBEGIN
using mat3x3 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
using quat   = Eigen::Quaternion<double>;
using vec3   = Eigen::Matrix<double, 3, 1>;
using vecX   = Eigen::VectorX<double>;
using aaxis  = Eigen::AngleAxisd;
// NOLINTEND

}  // namespace aos
