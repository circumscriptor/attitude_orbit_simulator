#pragma once

#include <Eigen/Dense>

namespace aos {

// NOLINTBEGIN
using mat3x3 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
using quat   = Eigen::Quaternion<double>;
using vec3   = Eigen::Matrix<double, 3, 1>;
using vecX   = Eigen::VectorX<double>;
// NOLINTEND

using inertia_t            = mat3x3;
using velocity_t           = vec3;
using rod_magnetizations_t = vecX;

}  // namespace aos
