// Copyright (C) 2013 The Regents of the University of California (Regents).
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
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <algorithm>
#include <limits>
#include <vector>
#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/util/random.h"
#include "theia/sfm/estimators/uncalibrated_relative_pose_estimator.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/test/test_utils.h"

namespace theia {
namespace {
using Eigen::AngleAxisd;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector2d;
using Eigen::Vector3d;

// Tests that the five point algorithm works correctly be checking the returned
// solutions and ensuring that at least one of the candidate solutions
// corresponds to the  we used to construct the problem.
void TestUncalibratedRelativePoseEstimator(
    const std::vector<Vector3d> points_3d,
    const Matrix3d& expected_rotation,
    const Vector3d& expected_translation,
    const double expected_focal_length1,
    const double expected_focal_length2,
    const double pose_tolerance,
    const double focal_length_tolerance) {
  const Vector3d gt_position =
      -expected_rotation.transpose() * expected_translation;

  // Calculates the image points in both views.
  std::vector<FeatureCorrespondence> features(8);
  for (int i = 0; i < points_3d.size(); ++i) {
    const Vector3d proj_3d =
        expected_rotation * points_3d[i] + expected_translation;
    features[i].feature1 = expected_focal_length1 * points_3d[i].hnormalized();
    features[i].feature2 = expected_focal_length2 * proj_3d.hnormalized();
  }

  // Compute the relative pose.
  UncalibratedRelativePoseEstimator relative_pose_estimator;
  std::vector<UncalibratedRelativePose> soln_relative_poses;
  EXPECT_TRUE(
      relative_pose_estimator.EstimateModel(features, &soln_relative_poses));
  CHECK_EQ(soln_relative_poses.size(), 1);

  // Expect the returned pose to match well.
  const double rotation_error =
      RadToDeg(AngleAxisd(soln_relative_poses[0].rotation.transpose() *
                          expected_rotation).angle());
  EXPECT_LT(rotation_error, pose_tolerance);

  const double position_error = RadToDeg(std::acos(
      std::min(1.0, std::max(-1.0, gt_position.normalized().dot(
                                       soln_relative_poses[0].position)))));
  EXPECT_LT(position_error, pose_tolerance);

  // Expect focal lengths to match well.
  EXPECT_NEAR(soln_relative_poses[0].focal_length1,
              expected_focal_length1,
              focal_length_tolerance);
  EXPECT_NEAR(soln_relative_poses[0].focal_length2,
              expected_focal_length2,
              focal_length_tolerance);
}


TEST(UncalibratedRelativePoseEstimator, BasicTest) {
  InitRandomGenerator();

  // Ground truth.
  const std::vector<Vector3d> points_3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0),
                                            Vector3d(-1.0, -3.0, 2.0),
                                            Vector3d(1.0, -2.0, 1.0),
                                            Vector3d(-1.0, 4.0, 2.0),
                                            Vector3d(-2.0, 2.0, 3.0)
  };

  for (int i = 0; i < 1000; i++) {
    const Matrix3d soln_rotation = ProjectToRotationMatrix(
        Matrix3d::Identity() + 0.3 * Matrix3d::Random());
    const Vector3d soln_translation = Vector3d::Random();
    const double soln_focal_length1 = RandDouble(800, 1600);
    const double soln_focal_length2 = RandDouble(800, 1600);

    const double kTolerance = 1e-4;
    TestUncalibratedRelativePoseEstimator(points_3d,
                                          soln_rotation,
                                          soln_translation.normalized(),
                                          soln_focal_length1,
                                          soln_focal_length2,
                                          kTolerance,
                                          kTolerance);
  }
}

}  // namespace
}  // namespace theia
