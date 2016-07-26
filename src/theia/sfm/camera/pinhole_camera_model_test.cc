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
//         Torsten Sattler (sattlert@inf.ethz.ch)

#include <Eigen/Dense>

#include <ceres/rotation.h>
#include <math.h>
#include "gtest/gtest.h"

#include "theia/alignment/alignment.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/sfm/camera/pinhole_camera_model.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

TEST(PinholeCameraModel, InternalParameterGettersAndSetters) {
  PinholeCameraModel camera;

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

void ReprojectionTest(const PinholeCameraModel& camera) {
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

TEST(PinholeCameraModel, ReprojectionNoDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeCameraModel camera;
  for (int i = 0; i < 1; i++) {
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0, 0);
    ReprojectionTest(camera);
  }
}

TEST(PinholeCameraModel, ReprojectionOneDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeCameraModel camera;
  for (int i = 0; i < 1; i++) {
    // Initialize a random camera.
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0.01, 0.0);
    ReprojectionTest(camera);
  }
}

TEST(PinholeCameraModel, ReprojectionTwoDistortion) {
  static const double kPrincipalPoint[2] = {600.0, 400.0};
  static const double kFocalLength = 1200;
  InitRandomGenerator();
  PinholeCameraModel camera;
  for (int i = 0; i < 100; i++) {
    // Initialize a random camera.
    camera.SetFocalLength(kFocalLength);
    camera.SetPrincipalPoint(kPrincipalPoint[0], kPrincipalPoint[1]);
    camera.SetRadialDistortion(0.01, 0.001);
    ReprojectionTest(camera);
  }
}

}  // namespace theia
