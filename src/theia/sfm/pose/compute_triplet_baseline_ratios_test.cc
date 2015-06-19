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


#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "gtest/gtest.h"

#include "theia/sfm/pose/compute_triplet_baseline_ratios.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/view_triplet.h"
#include "theia/util/random.h"

namespace theia {

using Eigen::Vector3d;

void RelativeRotationFromTwoRotations(const Vector3d& rotation1,
                                      const Vector3d& rotation2,
                                      Vector3d* relative_rotation) {
  Eigen::Matrix3d rotation_matrix1, rotation_matrix2;
  ceres::AngleAxisToRotationMatrix(rotation1.data(), rotation_matrix1.data());
  ceres::AngleAxisToRotationMatrix(rotation2.data(), rotation_matrix2.data());

  const Eigen::Matrix3d relative_rotation_mat =
      rotation_matrix2 * rotation_matrix1.transpose();
  ceres::RotationMatrixToAngleAxis(relative_rotation_mat.data(),
                                   relative_rotation->data());
}

void RelativeTranslationFromTwoPositions(const Vector3d& position1,
                                         const Vector3d& position2,
                                         const Vector3d& rotation1,
                                         Vector3d* translation) {
  Eigen::Matrix3d rotation_matrix1;
  ceres::AngleAxisToRotationMatrix(rotation1.data(), rotation_matrix1.data());
  *translation = rotation_matrix1 * (position2 - position1).normalized();
}

void TwoViewInfoFromCameras(const Camera& camera1,
                            const Camera& camera2,
                            TwoViewInfo* info) {
  RelativeRotationFromTwoRotations(camera1.GetOrientationAsAngleAxis(),
                                   camera2.GetOrientationAsAngleAxis(),
                                   &info->rotation_2);
  RelativeTranslationFromTwoPositions(camera1.GetPosition(),
                                      camera2.GetPosition(),
                                      camera1.GetOrientationAsAngleAxis(),
                                      &info->position_2);
}

// Cameras are generated within a 2x2x2 box around the origin looking
// approximately toward the negative z-axis. The scene depth indicates roughly
// how many units away from the origin the 3D points will be when created. The
// points will always be in front of all cameras.
void TestTripletBaselineComputation(const double pixel_noise,
                                    const double scene_depth,
                                    const double tolerance) {
  InitRandomGenerator();

  static const int kNumTracks = 100;
  static const double kFocalLength = 1000;

  // Set up 3 views.
  Camera camera1, camera2, camera3;
  camera1.SetPosition(Vector3d::Random());
  camera1.SetOrientationFromAngleAxis(0.2 * Vector3d::Random());
  camera2.SetPosition(Vector3d::Random());
  camera2.SetOrientationFromAngleAxis(0.2 * Vector3d::Random());
  camera3.SetPosition(Vector3d::Random());
  camera3.SetOrientationFromAngleAxis(0.2 * Vector3d::Random());

  // Add tracks.
  std::vector<Feature> feature1(kNumTracks), feature2(kNumTracks),
      feature3(kNumTracks);
  for (int i = 0; i < kNumTracks; i++) {
    const Eigen::Vector4d point = Eigen::Vector3d::Random().homogeneous() +
                                  Eigen::Vector4d(0, 0, scene_depth, 0);

    // Add the observations (plus noise if applicable).
    camera1.ProjectPoint(point, &feature1[i]);
    camera2.ProjectPoint(point, &feature2[i]);
    camera3.ProjectPoint(point, &feature3[i]);

    if (pixel_noise > 0) {
      feature1[i].x() += RandGaussian(0.0, pixel_noise / kFocalLength);
      feature1[i].y() += RandGaussian(0.0, pixel_noise / kFocalLength);
      feature2[i].x() += RandGaussian(0.0, pixel_noise / kFocalLength);
      feature2[i].y() += RandGaussian(0.0, pixel_noise / kFocalLength);
      feature3[i].x() += RandGaussian(0.0, pixel_noise / kFocalLength);
      feature3[i].y() += RandGaussian(0.0, pixel_noise / kFocalLength);
    }
  }

  // Compute the baseline ratios.
  Eigen::Vector3d baseline;
  ViewTriplet triplet;
  TwoViewInfoFromCameras(camera1, camera2, &triplet.info_one_two);
  TwoViewInfoFromCameras(camera1, camera3, &triplet.info_one_three);
  TwoViewInfoFromCameras(camera2, camera3, &triplet.info_two_three);
  ComputeTripletBaselineRatios(triplet, feature1, feature2, feature3,
                               &baseline);

  // Measure the error.
  const double baseline_12 =
      (camera2.GetPosition() - camera1.GetPosition()).norm();
  const double baseline_13 =
      (camera3.GetPosition() - camera1.GetPosition()).norm();
  const double baseline_23 =
      (camera3.GetPosition() - camera2.GetPosition()).norm();
  const Eigen::Vector3d gt_baseline(1.0,
                                    baseline_13 / baseline_12,
                                    baseline_23 / baseline_12);

  EXPECT_NEAR(gt_baseline[0] / baseline[0], 1.0, tolerance);
  EXPECT_NEAR(gt_baseline[1] / baseline[1], 1.0, tolerance);
  EXPECT_NEAR(gt_baseline[2] / baseline[2], 1.0, tolerance);
}

TEST(ComputeTripletBaselineRatios, NoNoiseSmallScene) {
  TestTripletBaselineComputation(0, 5, 1e-8);
}

TEST(ComputeTripletBaselineRatios, NoiseSmallScene) {
  TestTripletBaselineComputation(1, 5, 0.05);
}

TEST(ComputeTripletBaselineRatios, NoNoiseLargeScene) {
  TestTripletBaselineComputation(0, 100, 1e-8);
}

TEST(ComputeTripletBaselineRatios, NoiseLargeScene) {
  TestTripletBaselineComputation(1, 100, 0.1);
}

}  // namespace theia
