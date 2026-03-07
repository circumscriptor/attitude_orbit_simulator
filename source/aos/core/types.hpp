#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace toml::inline v3 {
class table;
}  // namespace toml::inline v3

namespace aos {

using real_t = double;

// NOLINTBEGIN
using mat3x3 = Eigen::Matrix<real_t, 3, 3, Eigen::RowMajor>;
using quat   = Eigen::Quaternion<real_t>;
using vec3   = Eigen::Matrix<real_t, 3, 1>;
using vec4   = Eigen::Matrix<real_t, 4, 1>;
using vecX   = Eigen::VectorX<real_t>;
using aaxis  = Eigen::AngleAxis<real_t>;
// NOLINTEND

using toml_table = toml::table;

}  // namespace aos
