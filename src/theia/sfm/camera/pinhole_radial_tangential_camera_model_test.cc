// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include <Eigen/Dense>

#include <ceres/rotation.h>
#include <math.h>
#include "gtest/gtest.h"

#include "theia/alignment/alignment.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/sfm/camera/pinhole_radial_tangential_camera_model.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

static const int kNumTrials = 100;

TEST(PinholeRadialTangentialCameraModel, InternalParameterGettersAndSetters) {
  PinholeRadialTangentialCameraModel camera;

  EXPECT_EQ(camera.Type(),
            CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL);

  // Check that default values are set
  EXPECT_EQ(camera.FocalLength(), 1.0);
  EXPECT_EQ(camera.AspectRatio(), 1.0);
  EXPECT_EQ(camera.Skew(), 0.0);
  EXPECT_EQ(camera.PrincipalPointX(), 0.0);
  EXPECT_EQ(camera.PrincipalPointY(), 0.0);
  EXPECT_EQ(camera.RadialDistortion1(), 0.0);
  EXPECT_EQ(camera.RadialDistortion2(), 0.0);
  EXPECT_EQ(camera.RadialDistortion3(), 0.0);
  EXPECT_EQ(camera.TangentialDistortion1(), 0.0);
  EXPECT_EQ(camera.TangentialDistortion2(), 0.0);

  // Set parameters to different values.
  camera.SetFocalLength(600.0);
  camera.SetAspectRatio(0.9);
  camera.SetSkew(0.01);
  camera.SetPrincipalPoint(300.0, 400.0);
  camera.SetRadialDistortion(0.01, 0.001, 0.0001);
  camera.SetTangentialDistortion(0.02, 0.002);

  // Check that the values were updated.
  EXPECT_EQ(camera.FocalLength(), 600.0);
  EXPECT_EQ(camera.AspectRatio(), 0.9);
  EXPECT_EQ(camera.Skew(), 0.01);
  EXPECT_EQ(camera.PrincipalPointX(), 300.0);
  EXPECT_EQ(camera.PrincipalPointY(), 400.0);
  EXPECT_EQ(camera.RadialDistortion1(), 0.01);
  EXPECT_EQ(camera.RadialDistortion2(), 0.001);
  EXPECT_EQ(camera.RadialDistortion3(), 0.0001);
  EXPECT_EQ(camera.TangentialDistortion1(), 0.02);
  EXPECT_EQ(camera.TangentialDistortion2(), 0.002);
}
// Test to ensure that the camera intrinsics are being set appropriately.
void TestSetFromCameraintrinsicsPrior(const CameraIntrinsicsPrior& prior) {
  const PinholeRadialTangentialCameraModel default_camera;
  PinholeRadialTangentialCameraModel camera;
  camera.SetFromCameraIntrinsicsPriors(prior);

  if (prior.focal_length.is_set) {
    EXPECT_EQ(camera.FocalLength(), prior.focal_length.value[0]);
  } else {
    EXPECT_EQ(camera.FocalLength(), default_camera.FocalLength());
  }

  if (prior.principal_point.is_set) {
    EXPECT_EQ(camera.PrincipalPointX(), prior.principal_point.value[0]);
    EXPECT_EQ(camera.PrincipalPointY(), prior.principal_point.value[1]);
  } else {
    EXPECT_EQ(camera.PrincipalPointX(), default_camera.PrincipalPointX());
    EXPECT_EQ(camera.PrincipalPointY(), default_camera.PrincipalPointY());
  }

  if (prior.aspect_ratio.is_set) {
    EXPECT_EQ(camera.AspectRatio(), prior.aspect_ratio.value[0]);
  } else {
    EXPECT_EQ(camera.AspectRatio(), default_camera.AspectRatio());
  }

  if (prior.skew.is_set) {
    EXPECT_EQ(camera.Skew(), prior.skew.value[0]);
  } else {
    EXPECT_EQ(camera.Skew(), default_camera.Skew());
  }

  if (prior.radial_distortion.is_set) {
    EXPECT_EQ(camera.RadialDistortion1(), prior.radial_distortion.value[0]);
    EXPECT_EQ(camera.RadialDistortion2(), prior.radial_distortion.value[1]);
  } else {
    EXPECT_EQ(camera.RadialDistortion1(), default_camera.RadialDistortion1());
    EXPECT_EQ(camera.RadialDistortion2(), default_camera.RadialDistortion2());
  }
}

// Gradually add one prior at a time and ensure that the method still works. We
// test before and after setting the "is_set" member variable to true to ensure
// that setting the value of priors when is_set=false is handled properly.
TEST(PinholeRadialTangentialCameraModel, SetFromCameraIntrinsicsPriors) {
  CameraIntrinsicsPrior prior;
  prior.focal_length.value[0] = 1000.0;
  prior.principal_point.value[0] = 400.0;
  prior.principal_point.value[1] = 300.0;
  prior.aspect_ratio.value[0] = 1.01;
  prior.skew.value[0] = 0.01;
  prior.radial_distortion.value[0] = 0.1;
  prior.radial_distortion.value[1] = 0.01;
  prior.radial_distortion.value[2] = 0.001;
  prior.tangential_distortion.value[0] = 0.1;
  prior.tangential_distortion.value[1] = 0.01;

  TestSetFromCameraintrinsicsPrior(prior);

  prior.focal_length.is_set = true;
  TestSetFromCameraintrinsicsPrior(prior);

  prior.principal_point.is_set = true;
  TestSetFromCameraintrinsicsPrior(prior);

  prior.aspect_ratio.is_set = true;
  TestSetFromCameraintrinsicsPrior(prior);

  prior.skew.is_set = true;
  TestSetFromCameraintrinsicsPrior(prior);

  prior.radial_distortion.is_set = true;
  TestSetFromCameraintrinsicsPrior(prior);

  prior.tangential_distortion.is_set = true;
  TestSetFromCameraintrinsicsPrior(prior);
}

TEST(PinholeRadialTangentialCameraModel, GetSubsetFromOptimizeIntrinsicsType) {
  PinholeRadialTangentialCameraModel camera;
  std::vector<int> constant_subset;

  constant_subset =
      camera.GetSubsetFromOptimizeIntrinsicsType(OptimizeIntrinsicsType::NONE);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters());

  // Test that optimizing for focal length works correctly.
  constant_subset = camera.GetSubsetFromOptimizeIntrinsicsType(
      OptimizeIntrinsicsType::FOCAL_LENGTH);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters() - 1);
  for (int i = 0; i < constant_subset.size(); i++) {
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::FOCAL_LENGTH);
  }

  // Test that optimizing for principal_points works correctly.
  constant_subset = camera.GetSubsetFromOptimizeIntrinsicsType(
      OptimizeIntrinsicsType::PRINCIPAL_POINTS);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters() - 2);
  for (int i = 0; i < constant_subset.size(); i++) {
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::PRINCIPAL_POINT_X);
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::PRINCIPAL_POINT_Y);
  }

  // Test that optimizing for aspect ratio works correctly.
  constant_subset = camera.GetSubsetFromOptimizeIntrinsicsType(
      OptimizeIntrinsicsType::ASPECT_RATIO);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters() - 1);
  for (int i = 0; i < constant_subset.size(); i++) {
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::ASPECT_RATIO);
  }

  // Test that optimizing for skew works correctly.
  constant_subset = camera.GetSubsetFromOptimizeIntrinsicsType(
      OptimizeIntrinsicsType::SKEW);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters() - 1);
  for (int i = 0; i < constant_subset.size(); i++) {
    EXPECT_NE(constant_subset[i], PinholeRadialTangentialCameraModel::SKEW);
  }

  // Test that optimizing for radial distortion works correctly.
  constant_subset = camera.GetSubsetFromOptimizeIntrinsicsType(
      OptimizeIntrinsicsType::RADIAL_DISTORTION);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters() - 3);
  for (int i = 0; i < constant_subset.size(); i++) {
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_1);
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_2);
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_3);
  }

  // Test that optimizing for tangential distortion does not optimize any
  // parameters.
  constant_subset = camera.GetSubsetFromOptimizeIntrinsicsType(
      OptimizeIntrinsicsType::TANGENTIAL_DISTORTION);
  EXPECT_EQ(constant_subset.size(), camera.NumParameters() - 2);
  for (int i = 0; i < constant_subset.size(); i++) {
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_1);
    EXPECT_NE(constant_subset[i],
              PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_2);
  }
}

void ReprojectionTest(const PinholeRadialTangentialCameraModel& camera) {
  const double kTolerance = 1e-5;
  const double kImageWidth = 1200.0;

  for (int i = 0; i < 10; i++) {
    // Get a random pixel within the image.
    const Vector2d pixel =
        kImageWidth * (Vector2d::Random() + Vector2d::Ones()) / 2.0;

    // Get the normalized ray of that pixel.
    const Vector3d normalized_ray = camera.ImageToCameraCoordinates(pixel);

    const double random_depth = RandDouble(0.01, 100.0);
    const Vector3d random_point = normalized_ray * random_depth;
    const Vector2d reprojected_pixel =
        camera.CameraToImageCoordinates(random_point);

    // Expect the reprojection to be close.
    EXPECT_LT((pixel - reprojected_pixel).norm(), kTolerance)
        << "gt pixel: " << pixel.transpose()
        << "\nreprojected pixel: " << reprojected_pixel.transpose();
  }
}

TEST(PinholeRadialTangentialCameraModel, ReprojectionNoDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeRadialTangentialCameraModel camera;
  for (int i = 0; i < kNumTrials; i++) {
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0, 0, 0);
    ReprojectionTest(camera);
  }
}

TEST(PinholeRadialTangentialCameraModel, ReprojectionOneRadialDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeRadialTangentialCameraModel camera;
  for (int i = 0; i < kNumTrials; i++) {
    // Initialize a random camera.
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0.01, 0.0, 0.0);
    ReprojectionTest(camera);
  }
}

TEST(PinholeRadialTangentialCameraModel, ReprojectionTwoRadialDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeRadialTangentialCameraModel camera;
  for (int i = 0; i < kNumTrials; i++) {
    // Initialize a random camera.
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0.01, 0.001, 0.0);
    ReprojectionTest(camera);
  }
}

TEST(PinholeRadialTangentialCameraModel, ReprojectionOneTangentialDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeRadialTangentialCameraModel camera;
  for (int i = 0; i < kNumTrials; i++) {
    // Initialize a random camera.
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0.01, 0.001, 0.0001);
    camera.SetTangentialDistortion(0.01, 0.0);
    ReprojectionTest(camera);
  }
}

TEST(PinholeRadialTangentialCameraModel, ReprojectionTwoTangentialDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeRadialTangentialCameraModel camera;
  for (int i = 0; i < kNumTrials; i++) {
    // Initialize a random camera.
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0.01, 0.001, 0.0001);
    camera.SetTangentialDistortion(0.01, 0.001);
    ReprojectionTest(camera);
  }
}

}  // namespace theia
