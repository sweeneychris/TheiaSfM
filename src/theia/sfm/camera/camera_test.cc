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
//         Torsten Sattler (sattlert@inf.ethz.ch)

#include <Eigen/Dense>

#include <ceres/rotation.h>
#include <math.h>
#include "gtest/gtest.h"

#include "theia/alignment/alignment.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/sfm/camera/camera.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

TEST(Camera, ProjectionMatrix) {
  const double kTolerance = 1e-12;

  Camera camera;
  const double image_size = 500;
  for (int i = 0; i < 100; i++) {
    const Matrix3x4d gt_projection_matrix = Matrix3x4d::Random();
    EXPECT_TRUE(camera.InitializeFromProjectionMatrix(image_size,
                                                      image_size,
                                                      gt_projection_matrix));
    Matrix3x4d projection_matrix;
    camera.GetProjectionMatrix(&projection_matrix);

    EXPECT_TRUE(test::ArraysEqualUpToScale(12,
                                           gt_projection_matrix.data(),
                                           projection_matrix.data(),
                                           kTolerance));
  }
}

TEST(Camera, InternalParameterGettersAndSetters) {
  Camera camera;

  // Check that default values are set
  EXPECT_EQ(camera.FocalLength(), 1.0);
  EXPECT_EQ(camera.AspectRatio(), 1.0);
  EXPECT_EQ(camera.Skew(), 0.0);
  EXPECT_EQ(camera.PrincipalPointX(), 0.0);
  EXPECT_EQ(camera.PrincipalPointY(), 0.0);
  EXPECT_EQ(camera.RadialDistortion1(), 0.0);
  EXPECT_EQ(camera.RadialDistortion2(), 0.0);

  // Set parameters to different values.
  camera.SetFocalLength(600.0);
  camera.SetAspectRatio(0.9);
  camera.SetSkew(0.01);
  camera.SetPrincipalPoint(300.0, 400.0);
  camera.SetRadialDistortion(0.01, 0.001);

  // Check that the values were updated.
  EXPECT_EQ(camera.FocalLength(), 600.0);
  EXPECT_EQ(camera.AspectRatio(), 0.9);
  EXPECT_EQ(camera.Skew(), 0.01);
  EXPECT_EQ(camera.PrincipalPointX(), 300.0);
  EXPECT_EQ(camera.PrincipalPointY(), 400.0);
  EXPECT_EQ(camera.RadialDistortion1(), 0.01);
  EXPECT_EQ(camera.RadialDistortion2(), 0.001);
}

TEST(Camera, ExternalParameterGettersAndSetters) {
  const double kTolerance = 1e-16;

  Camera camera;

  // Check that the default values are set.
  EXPECT_DOUBLE_EQ(camera.GetPosition().squaredNorm(), 0.0);
  EXPECT_DOUBLE_EQ(camera.GetOrientationAsAngleAxis().squaredNorm(), 0.0);
  EXPECT_DOUBLE_EQ((camera.GetOrientationAsRotationMatrix() -
                    Matrix3d::Identity()).squaredNorm(),
                   0.0);

  // Check that position getter/setters work.
  camera.SetPosition(Vector3d::Ones());
  EXPECT_DOUBLE_EQ((camera.GetPosition() - Vector3d::Ones()).squaredNorm(),
                   0.0);

  // Check that angle axis getter/setters work.
  Vector3d gt_angle_axis(1.0, 1.0, 1.0);
  gt_angle_axis = Vector3d(0.3, 0.7, 0.4);
  Matrix3d gt_rotation_matrix;
  ceres::AngleAxisToRotationMatrix(gt_angle_axis.data(),
                                   gt_rotation_matrix.data());
  camera.SetOrientationFromRotationMatrix(gt_rotation_matrix);
  EXPECT_LT(
      (camera.GetOrientationAsAngleAxis() - gt_angle_axis).squaredNorm(),
      kTolerance);
  EXPECT_LT((camera.GetOrientationAsRotationMatrix() - gt_rotation_matrix)
                .squaredNorm(),
              kTolerance);

  // Check that rotation matrix getter/setters work.
  gt_angle_axis = Vector3d(0.3, 0.7, 0.4);
  ceres::AngleAxisToRotationMatrix(gt_angle_axis.data(),
                                   gt_rotation_matrix.data());
  camera.SetOrientationFromRotationMatrix(gt_rotation_matrix);
  EXPECT_LT(
      (camera.GetOrientationAsAngleAxis() - gt_angle_axis).squaredNorm(),
      kTolerance);
  EXPECT_LT((camera.GetOrientationAsRotationMatrix() - gt_rotation_matrix)
                .squaredNorm(),
              kTolerance);
}


void ReprojectionTest(const Camera& camera) {
  const double kTolerance = 1e-5;

  for (int i = 0; i < 10; i++) {
    // Get a random pixel within the image.
    const Vector2d pixel =
        camera.ImageWidth() * (Vector2d::Random() + Vector2d::Ones()) / 2.0;

    // Get the normalized ray of that pixel.
    const Vector3d normalized_ray = camera.PixelToUnitDepthRay(pixel);

    const double random_depth = RandDouble(0.01, 100.0);
    const Vector4d random_point =
        (camera.GetPosition() + normalized_ray * random_depth)
            .homogeneous();

    Vector2d reprojected_pixel;
    const double depth =
        camera.ProjectPoint(random_point, &reprojected_pixel);

    // Expect the reconstructed 3d points to be close.
    EXPECT_LT(std::abs(random_depth - depth),
              kTolerance * random_depth)
        << "real depth = " << random_depth
        << " and reconstructed depth = " << depth;

    // Expect the reprojection to be close.
    EXPECT_LT((pixel - reprojected_pixel).norm(), kTolerance)
        << "gt pixel: " << pixel.transpose()
        << "\nreprojected pixel: " << reprojected_pixel.transpose();
  }
}

TEST(Camera, ReprojectionNoDistortion) {
  InitRandomGenerator();
  Camera camera;
  const double image_size = 600;
  for (int i = 0; i < 1; i++) {
    // Initialize a random camera.
    camera.InitializeFromProjectionMatrix(image_size, image_size,
                                          Matrix3x4d::Random());

    // Initialize random positive radial distortion parameters.
    camera.SetRadialDistortion(0, 0);
    ReprojectionTest(camera);
  }
}

TEST(Camera, ReprojectionOneDistortion) {
  InitRandomGenerator();
  Camera camera;
  const double image_size = 600;
  for (int i = 0; i < 1; i++) {
    // Initialize a random camera.
    camera.InitializeFromProjectionMatrix(image_size, image_size,
                                          Matrix3x4d::Random());

    camera.SetRadialDistortion(RandDouble(0.0, 0.2), 0.0);
    ReprojectionTest(camera);
  }
}

TEST(Camera, ReprojectionTwoDistortion) {
    InitRandomGenerator();
  Camera camera;
  const double image_size = 600;
  for (int i = 0; i < 1; i++) {
    // Initialize a random camera.
    camera.InitializeFromProjectionMatrix(image_size, image_size,
                                          Matrix3x4d::Random());

    // Initialize random positive radial distortion parameters.
    camera.SetRadialDistortion(RandDouble(0.0, 0.2),
                               RandDouble(0.0, 0.02));
    ReprojectionTest(camera);
  }
}

}  // namespace theia
