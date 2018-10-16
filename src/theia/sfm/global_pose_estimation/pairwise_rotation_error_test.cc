// Copyright (C) 2014 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "gtest/gtest.h"
#include "theia/sfm/global_pose_estimation/pairwise_rotation_error.h"
#include "theia/math/util.h"

namespace theia {
namespace {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Vector3d;

static const double kRelativeRotationWeight = 1.0;

void PairwiseRotationErrorTest(const Matrix3d& relative_rot_mat,
                               const double weight,
                               const Matrix3d& global_rot_1_mat,
                               const Matrix3d& global_rot_2_mat) {
  Vector3d relative_rot, rotation1, rotation2;
  ceres::RotationMatrixToAngleAxis(relative_rot_mat.data(),
                                   relative_rot.data());
  ceres::RotationMatrixToAngleAxis(global_rot_1_mat.data(), rotation1.data());
  ceres::RotationMatrixToAngleAxis(global_rot_2_mat.data(), rotation2.data());

  // Compute ground truth angular error.
  const AngleAxisd loop_rotation(global_rot_2_mat *
                                 global_rot_1_mat.transpose());
  const Matrix3d gt_rotation_error_mat = loop_rotation * relative_rot_mat.transpose();

  // Convert matrices to angle-axis.
  Vector3d gt_rotation_error;
  ceres::RotationMatrixToAngleAxis(gt_rotation_error_mat.data(), gt_rotation_error.data());
  gt_rotation_error *= weight;

  // Initialize error function and compute rotation error.
  const PairwiseRotationError global_rotation_error(relative_rot, weight);
  Vector3d estimated_rotation_error;
  global_rotation_error(rotation1.data(), rotation2.data(),
                        estimated_rotation_error.data());

  static const double kTolerance = 1e-12;
  EXPECT_NEAR(estimated_rotation_error(0), gt_rotation_error(0), kTolerance);
  EXPECT_NEAR(estimated_rotation_error(1), gt_rotation_error(1), kTolerance);
  EXPECT_NEAR(estimated_rotation_error(2), gt_rotation_error(2), kTolerance);
}

}  // namespace

TEST(PairwiseRotationError, SmallRotation) {
  const Matrix3d global_rot_1_mat = Matrix3d::Identity();
  const Matrix3d global_rot_2_mat =
      AngleAxisd(DegToRad(2.0), Vector3d::UnitZ()).toRotationMatrix();
  const Matrix3d relative_rot =
      AngleAxisd(DegToRad(1.0), Vector3d::UnitZ()).toRotationMatrix();

  PairwiseRotationErrorTest(relative_rot, kRelativeRotationWeight,
                            global_rot_1_mat, global_rot_2_mat);
}

TEST(PairwiseRotationError, NontrivialRotation) {
  const Matrix3d global_rot_1_mat = Matrix3d::Identity();
  const Matrix3d global_rot_2_mat =
      (AngleAxisd(DegToRad(5.3), Vector3d::UnitX()) *
       AngleAxisd(DegToRad(1.2), Vector3d::UnitY()) *
       AngleAxisd(DegToRad(8.1), Vector3d::UnitZ())).toRotationMatrix();
  const Matrix3d relative_rot =
      (AngleAxisd(DegToRad(5.9), Vector3d::UnitX()) *
       AngleAxisd(DegToRad(1.8), Vector3d::UnitY()) *
       AngleAxisd(DegToRad(7.6), Vector3d::UnitZ())).toRotationMatrix();

  PairwiseRotationErrorTest(relative_rot, kRelativeRotationWeight,
                            global_rot_1_mat, global_rot_2_mat);
}

TEST(PairwiseRotationError, OneHundredEightyDegreeRotation) {
  const Matrix3d global_rot_1_mat = Matrix3d::Identity();
  const Matrix3d global_rot_2_mat =
      AngleAxisd(DegToRad(179.0), Vector3d::UnitZ()).toRotationMatrix();
  const Matrix3d relative_rot =
      AngleAxisd(DegToRad(-179.0), Vector3d::UnitZ()).toRotationMatrix();

  PairwiseRotationErrorTest(relative_rot, kRelativeRotationWeight,
                            global_rot_1_mat, global_rot_2_mat);
  PairwiseRotationErrorTest(relative_rot, kRelativeRotationWeight,
                            global_rot_2_mat, global_rot_1_mat);
}

TEST(PairwiseRotationError, Weight) {
  const Matrix3d global_rot_1_mat = Matrix3d::Identity();
  const Matrix3d global_rot_2_mat =
      (AngleAxisd(DegToRad(5.3), Vector3d::UnitX()) *
       AngleAxisd(DegToRad(1.2), Vector3d::UnitY()) *
       AngleAxisd(DegToRad(8.1), Vector3d::UnitZ())).toRotationMatrix();
  const Matrix3d relative_rot =
      (AngleAxisd(DegToRad(5.9), Vector3d::UnitX()) *
       AngleAxisd(DegToRad(1.8), Vector3d::UnitY()) *
       AngleAxisd(DegToRad(7.6), Vector3d::UnitZ())).toRotationMatrix();

  static const double kNontrivialRelativeRotationWeight = 2.0;
  PairwiseRotationErrorTest(relative_rot, kNontrivialRelativeRotationWeight,
                            global_rot_1_mat, global_rot_2_mat);
}

}  // namespace theia
