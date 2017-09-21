// Copyright (C) 2017 The Regents of the University of California (Regents).
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#include "theia/math/rotation.h"

#include <ceres/rotation.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <limits>

namespace theia {
// Eigen::Vector3d MultiplyRotations(const Eigen::Vector3d& rotation1,
//                                   const Eigen::Vector3d& rotation2) {
//   const double theta1_sq = rotation1.squaredNorm();
//   const double theta2_sq = rotation2.squaredNorm();

//   // Compute the sin and cosine terms below. Take care to ensure that there will
//   // not be any divide-by-zeros when the rotations are small.
//   double cos_a;
//   double cos_b;
//   Eigen::Vector3d sin_a_times_v1;
//   Eigen::Vector3d sin_b_times_v2;

//   // We need to explicity handle the cases when the rotations are near zero.
//   // Near zero, the first order Taylor approximation of the rotation matrix R
//   // corresponding to a vector w and angle w is
//   //
//   //   R = I + hat(w) * sin(theta)
//   //
//   // But sintheta ~ theta and theta * w = angle_axis, which gives us
//   //
//   //  R = I + hat(angle_axis)
//   //
//   // We will use this in the special cases below to ensure stable computation
//   // when composing the rotations and avoid dividing by zero when theta is
//   // small.
//   if (theta1_sq < std::numeric_limits<double>::epsilon()) {
//     // Use the fact that theta ~ sin(theta) when theta is small. Since
//     // a = 0.5 * theta1, we end up with:
//     //   sin(a) * v1 = sin(theta / 2) * v1 = theta / 2 * v1 = rotation1 / 2.
//     sin_a_times_v1 = 0.5 * rotation1;
//     // When a is small, cos(a) ~ 1.0 by the first order taylor approximation.
//     cos_a = 1.0;
//   } else {
//     const double theta1 = std::sqrt(theta1_sq);
//     const double sin_a = std::sin(0.5 * theta1);
//     cos_a = std::cos(0.5 * theta1);
//     sin_a_times_v1 = sin_a * rotation1 / theta1;
//   }

//   // Same as above, but for theta2.
//   if (theta2_sq < std::numeric_limits<double>::epsilon()) {
//     sin_b_times_v2 = 0.5 * rotation2;
//     cos_b = 1.0;
//   } else {
//     const double theta2 = std::sqrt(theta2_sq);
//     const double sin_b = std::sin(0.5 * theta2);
//     cos_b = std::cos(0.5 * theta2);
//     sin_b_times_v2 = sin_b * rotation2 / theta2;
//   }

//   // Compute sin(c) * v3 using the formula above.
//   const Eigen::Vector3d sin_c_times_v3 = cos_b * sin_a_times_v1 +
//                                          cos_a * sin_b_times_v2 +
//                                          sin_a_times_v1.cross(sin_b_times_v2);

//   // If sin(c) is near zero then we again need to take care to avoid dividing by
//   // zero. We can use the first order Taylor approximation again, noting that
//   // sin(c) ~ c, which gives us:
//   //   rotation3 = theta * v3 = 2 * c * v3 ~ 2 * sin(c) * v3
//   const double sin_c_sq = sin_c_times_v3.squaredNorm();
//   if (sinc_c < std::numeric_limits<double>::epsilon()) {
//     const double diff = (2.0 * sin_c_times_v3 - rotation_aa).norm();
//     return 2.0 * sin_c_times_v3;
//   } else {
//     // Otherwise, we use the formula above. The angle axis rotation is the axis
//     // (v3) times the angle theta3.
//     const double sin_c = std::sqrt(sin_c_sq);
//     const double theta3 = 2.0 * std::asin(sin_c);
//     const Eigen::Vector3d v3 = sin_c_times_v3 / sin_c;

//     return theta3 * v3;
//   }
// }

// Use Ceres to perform a stable composition of rotations. This is not as
// efficient as directly composing angle axis vectors (see the old
// implementation commented above) but is more stable.
Eigen::Vector3d MultiplyRotations(const Eigen::Vector3d& rotation1,
                                  const Eigen::Vector3d& rotation2) {
  Eigen::Matrix3d rotation1_mat, rotation2_mat;
  ceres::AngleAxisToRotationMatrix(rotation1.data(), rotation1_mat.data());
  ceres::AngleAxisToRotationMatrix(rotation2.data(), rotation2_mat.data());

  const Eigen::Matrix3d rotation = rotation1_mat * rotation2_mat;
  Eigen::Vector3d rotation_aa;
  ceres::RotationMatrixToAngleAxis(rotation.data(), rotation_aa.data());
  return rotation_aa;
}

}  // namespace theia
