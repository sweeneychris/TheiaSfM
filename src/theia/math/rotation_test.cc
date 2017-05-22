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

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "gtest/gtest.h"

#include "theia/math/rotation.h"
#include "theia/math/util.h"

namespace theia {

void TestTwoRotations(const Eigen::Vector3d& rotation1,
                      const Eigen::Vector3d& rotation2,
                      const double tolerance_degrees) {
  // Convert to rotation matrices.
  Eigen::Matrix3d rotation1_matrix, rotation2_matrix;
  ceres::AngleAxisToRotationMatrix(
      rotation1.data(), ceres::ColumnMajorAdapter3x3(rotation1_matrix.data()));
  ceres::AngleAxisToRotationMatrix(
      rotation2.data(), ceres::ColumnMajorAdapter3x3(rotation2_matrix.data()));

  // Compose the rotations to obtain R3 = R1 * R2;
  const Eigen::Matrix3d rotation_matrix = rotation1_matrix * rotation2_matrix;

  // Get the estimated rotation.
  const Eigen::Vector3d estimated_rotation =
      MultiplyRotations(rotation1, rotation2);

  // Convert the estimated rotation to a rotation matrix.
  Eigen::Matrix3d estimated_rotation_matrix;
  ceres::AngleAxisToRotationMatrix(
      estimated_rotation.data(),
      ceres::ColumnMajorAdapter3x3(estimated_rotation_matrix.data()));

  // Determine the angle between the two rotations.
  const Eigen::AngleAxisd rotation_difference(rotation_matrix.transpose() *
                                              estimated_rotation_matrix);
  EXPECT_LT(RadToDeg(rotation_difference.angle()), tolerance_degrees);
}

TEST(MultiplyRotations, BasicTest) {
  static const double kToleranceDegrees = 1e-8;
  const Eigen::Vector3d rotation1(-0.3, -0.2, 0.1);
  const Eigen::Vector3d rotation2(0.13, 0.06, -0.4);
  TestTwoRotations(rotation1, rotation2, kToleranceDegrees);
}

TEST(MultiplyRotations, SmallFirstRotation) {
  static const double kToleranceDegrees = 1e-8;
  const Eigen::Vector3d rotation1(1e-12, 4e-10, 1e-16);
  const Eigen::Vector3d rotation2(0.13, 0.06, -0.4);
  TestTwoRotations(rotation1, rotation2, kToleranceDegrees);
}

TEST(MultiplyRotations, SmallSecondRotation) {
  static const double kToleranceDegrees = 1e-8;
  const Eigen::Vector3d rotation1(-0.3, -0.2, 0.1);
  const Eigen::Vector3d rotation2(1e-12, 4e-10, 1e-16);
  TestTwoRotations(rotation1, rotation2, kToleranceDegrees);
}

TEST(MultiplyRotations, BothSmallRotations) {
  static const double kToleranceDegrees = 1e-8;
  const Eigen::Vector3d rotation1(-4e-12, -2e-10, 1e-16);
  const Eigen::Vector3d rotation2(3e-12, 1e-10, 8e-16);
  TestTwoRotations(rotation1, rotation2, kToleranceDegrees);
}

TEST(MultiplyRotations, OppositeRotations) {
  static const double kToleranceDegrees = 1e-8;
  const Eigen::Vector3d rotation1(-0.3, -0.2, 0.1);
  const Eigen::Vector3d rotation2(0.3, 0.2, -0.1);
  TestTwoRotations(rotation1, rotation2, kToleranceDegrees);
}

}  // namespace theia
