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

#include "theia/sfm/pose/fundamental_matrix_util.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

#include "gtest/gtest.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/types.h"
#include "theia/util/random.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

RandomNumberGenerator rng(51);

TEST(FundamentalMatrixUtil, FocalLengths) {
  static const double kTolerance = 1e-6;

  for (int i = 0; i < 100; i++) {
    const double gt_focal_length1 = 800.0;
    const double gt_focal_length2 = 1000.0;
    const Vector3d rotation_angle_axis = rng.RandVector3d();
    const Matrix3d rotation =
        Eigen::AngleAxisd(rotation_angle_axis.norm(),
                          rotation_angle_axis.normalized()).toRotationMatrix();
    const Vector3d translation = rng.RandVector3d().normalized();

    // Create calibration matrices.
    Matrix3d fundamental_matrix;
    ComposeFundamentalMatrix(gt_focal_length1, gt_focal_length2,
                             rotation.data(), translation.data(),
                             fundamental_matrix.data());

    double focal_length1, focal_length2;
    EXPECT_TRUE(FocalLengthsFromFundamentalMatrix(
        fundamental_matrix.data(), &focal_length1, &focal_length2));
    EXPECT_NEAR(focal_length1, gt_focal_length1, kTolerance);
    EXPECT_NEAR(focal_length2, gt_focal_length2, kTolerance);
  }
}

TEST(FundamentalMatrixUtil, SharedFocalLengthsZeroIntrinsics) {
  static const double kTolerance = 1e-4;

  for (int i = 0; i < 100; i++) {
    const double gt_focal_length = rng.RandDouble(800, 1600);
    const Vector3d rotation_angle_axis = rng.RandVector3d();
    const Matrix3d rotation =
        Eigen::AngleAxisd(rotation_angle_axis.norm(),
                          rotation_angle_axis.normalized())
            .toRotationMatrix();
    const Vector3d translation = rng.RandVector3d().normalized();

    // Create calibration matrices.
    Matrix3d fundamental_matrix;
    ComposeFundamentalMatrix(gt_focal_length,
                             gt_focal_length,
                             rotation.data(),
                             translation.data(),
                             fundamental_matrix.data());

    double focal_length;
    EXPECT_TRUE(SharedFocalLengthsFromFundamentalMatrix(
        fundamental_matrix.data(), &focal_length));
    EXPECT_NEAR(focal_length, gt_focal_length, kTolerance);
  }
}

TEST(FundamentalMatrixUtil, FundamentalMatrixFromProjectionMatrices) {
  static const double kTolerance = 1e-12;
  static const int kNumPoints = 10;

  for (int i = 0; i < 100; i++) {
    // Set up model points.
    Vector3d points_3d[kNumPoints];
    for (int j = 0; j < kNumPoints; j++) {
      points_3d[j] = rng.RandVector3d() + Vector3d(0, 0, 10);
    }

    // Set up projection matrices.
    const Vector3d rotation_angle_axis = rng.RandVector3d();
    const Matrix3d rotation =
        Eigen::AngleAxisd(rotation_angle_axis.norm(),
                          rotation_angle_axis.normalized()).toRotationMatrix();
    const Vector3d translation = rng.RandVector3d();
    Matrix3x4d pmatrix1, pmatrix2;
    pmatrix1 << Matrix3d::Identity(), Vector3d::Zero();
    pmatrix2 << rotation, translation;

    // Get the fundamental matrix.
    Matrix3d fmatrix;
    FundamentalMatrixFromProjectionMatrices(pmatrix1.data(), pmatrix2.data(),
                                            fmatrix.data());
    // Ensure the epipolar constraint holds.
    for (int j = 0; j < kNumPoints; j++) {
      const Vector3d image_point1 = pmatrix1 * points_3d[j].homogeneous();
      const Vector3d image_point2 = pmatrix2 * points_3d[j].homogeneous();

      EXPECT_LT(fabs(image_point1.transpose() * fmatrix * image_point2),
                kTolerance);
    }
  }
}

}  // namespace theia
