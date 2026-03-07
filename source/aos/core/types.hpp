#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace toml::inline v3 {
class table;
}  // namespace toml::inline v3

namespace aos {

using real = double;

// NOLINTBEGIN
using mat3x3 = Eigen::Matrix<real, 3, 3, Eigen::RowMajor>;
using quat   = Eigen::Quaternion<real>;
using vec3   = Eigen::Matrix<real, 3, 1>;
using vec4   = Eigen::Matrix<real, 4, 1>;
using vecX   = Eigen::VectorX<real>;
using aaxis  = Eigen::AngleAxis<real>;
// NOLINTEND

using toml_table = toml::table;

}  // namespace aos
