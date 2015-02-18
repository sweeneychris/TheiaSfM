// Copyright (C) 2015 The Regents of the University of California (Regents).
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
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/sfm/transformation/align_rotations.h"

namespace theia {

namespace {

void ApplyRotation(const Eigen::Matrix3d& rotation_transformation,
                   const double noise,
                   Eigen::Vector3d* rotation) {
  const Eigen::Matrix3d noisy_rotation =
      Eigen::AngleAxisd(DegToRad(noise), Eigen::Vector3d::Random().normalized())
          .toRotationMatrix();

  // Apply the transformation to the rotation.
  Eigen::Matrix3d rotation_mat;
  ceres::AngleAxisToRotationMatrix(
      rotation->data(), ceres::ColumnMajorAdapter3x3(rotation_mat.data()));
  const Eigen::Matrix3d transformed_rotation =
      rotation_mat * (noisy_rotation * rotation_transformation);

  // Convert back to angle axis.
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(transformed_rotation.data()),
      rotation->data());
}

void TestAlignRotations(const int num_views,
                        const double noise_degrees,
                        const double tolerance) {
  std::vector<Eigen::Vector3d> gt_rotations(num_views);
  std::vector<Eigen::Vector3d> rotations(num_views);

  Eigen::Matrix3d rotation_transformation = Eigen::AngleAxisd(
      15.0, Eigen::Vector3d::Random().normalized()).toRotationMatrix();
  for (int i = 0; i < num_views; i++) {
    gt_rotations[i] = Eigen::Vector3d::Random();
    rotations[i] = gt_rotations[i];
    ApplyRotation(rotation_transformation, noise_degrees, &rotations[i]);
  }

  AlignRotations(gt_rotations, &rotations);

  for (int i = 0; i < num_views; i++) {
    EXPECT_LT((gt_rotations[i] - rotations[i]).norm(), tolerance);
  }
}

}  // namespace

TEST(AlignRotations, NoNoise) {
  static const int kNumViews = 20;
  static const double kTolerance = 1e-8;
  static const double kNoise = 0.0;
  TestAlignRotations(kNumViews, kNoise, kTolerance);
}

TEST(AlignRotations, Noise) {
  static const int kNumViews = 20;
  static const double kTolerance = 5e-2;
  static const double kNoise = 1.0;
  TestAlignRotations(kNumViews, kNoise, kTolerance);
}

}  // namespace theia
