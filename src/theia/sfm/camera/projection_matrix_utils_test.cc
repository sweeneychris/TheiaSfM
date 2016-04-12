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

#include <Eigen/Core>

#include "gtest/gtest.h"
#include "theia/test/test_utils.h"
#include "theia/sfm/camera/projection_matrix_utils.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

TEST(CameraIntrinsics, IntrinsicsToCalibrationMatrix) {
  const double focal_length = 800.0;
  const double skew = 1.0;
  const double aspect_ratio = 1.25;
  const double principal_point_x = 320;
  const double principal_point_y = 240;

  Matrix3d gt_calibration_matrix;
  gt_calibration_matrix <<
      focal_length, skew, principal_point_x,
      0, focal_length * aspect_ratio, principal_point_y,
      0, 0, 1;

  Matrix3d calibration_matrix;
  IntrinsicsToCalibrationMatrix(focal_length,
                                skew,
                                aspect_ratio,
                                principal_point_x, principal_point_y,
                                &calibration_matrix);
  for (int i = 0; i < calibration_matrix.rows(); i++) {
    for (int j = 0; j < calibration_matrix.cols(); j++) {
      EXPECT_DOUBLE_EQ(gt_calibration_matrix(i, j), calibration_matrix(i, j));
    }
  }
}

TEST(CameraIntrinsics, CalibrationMatrixToIntrinsics) {
  Matrix3d calibration_matrix;
  calibration_matrix <<
      800, 1.0, 320,
      0, 1000, 240,
      0, 0, 1;
  double focal_length, skew, aspect_ratio, principal_point_x, principal_point_y;
  CalibrationMatrixToIntrinsics(calibration_matrix,
                                &focal_length,
                                &skew,
                                &aspect_ratio,
                                &principal_point_x,
                                &principal_point_y);
  EXPECT_DOUBLE_EQ(focal_length, 800);
  EXPECT_DOUBLE_EQ(skew, 1.0);
  EXPECT_DOUBLE_EQ(aspect_ratio, 1.25);
  EXPECT_DOUBLE_EQ(principal_point_x, 320);
  EXPECT_DOUBLE_EQ(principal_point_y, 240);
}


TEST(DecomposeProjectMatrix, Random) {
  for (int i = 0; i < 1000; i++) {
    const Matrix3x4d random_pmatrix = Matrix3x4d::Random();
    Matrix3d calibration_matrix;
    Matrix3d rotation_matrix;
    Vector3d rotation, position;
    EXPECT_TRUE(DecomposeProjectionMatrix(random_pmatrix,
                                          &calibration_matrix,
                                          &rotation_matrix,
                                          &position));
  }
}

TEST(ComposeProjectionMatrix, Identity) {
  const Matrix3d calibration_matrix = Matrix3d::Identity();
  const Matrix3d rotation_matrix = Matrix3d::Identity();
  const Vector3d position = Vector3d::Zero();
  Matrix3x4d projection_matrix;
  EXPECT_TRUE(ComposeProjectionMatrix(calibration_matrix,
                                      rotation_matrix,
                                      position,
                                      &projection_matrix));
  EXPECT_TRUE(projection_matrix == Matrix3x4d::Identity());
}

TEST(ComposeProjectionMatrix, Random) {
  for (int i = 0; i < 1000; i++) {
    // Generate a random rotation.
    const Vector3d rotation_angle_axis = Vector3d::Random();
    const Matrix3d rotation_matrix = Eigen::AngleAxisd(rotation_angle_axis.norm(), rotation_angle_axis.normalized()).toRotationMatrix();
    const Vector3d position = Vector3d::Random();

    Matrix3d calibration_matrix;
    const double focal_length = 800;
    const double principal_point[2] = { 400, 300 };
    IntrinsicsToCalibrationMatrix(focal_length,
                                  0.0,
                                  1.0,
                                  principal_point[0],
                                  principal_point[1],
                                  &calibration_matrix);

    Matrix3x4d projection_matrix;
    EXPECT_TRUE(ComposeProjectionMatrix(calibration_matrix,
                                        rotation_matrix,
                                        position,
                                        &projection_matrix));
  }
}

TEST(ProjectionMatrix, Consistency) {
  const double kTolerance = 1e-6;

  for (int i = 0; i < 1000; i++) {
    const Matrix3x4d random_pmatrix = Matrix3x4d::Random();
    Matrix3d calibration_matrix;
    Matrix3d rotation_matrix;
    Vector3d position;
    EXPECT_TRUE(DecomposeProjectionMatrix(random_pmatrix,
                                          &calibration_matrix,
                                          &rotation_matrix,
                                          &position));

    Matrix3x4d composed_pmatrix;
    EXPECT_TRUE(ComposeProjectionMatrix(calibration_matrix,
                                        rotation_matrix,
                                        position,
                                        &composed_pmatrix));

    EXPECT_TRUE(test::ArraysEqualUpToScale(12, random_pmatrix.data(),
                                           composed_pmatrix.data(), kTolerance))
        << "gt pmatrix: \n" << random_pmatrix << "\ncomposed pmatrix: \n"
        << composed_pmatrix;
  }
}

}  // namespace theia
