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
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/estimators/relative_pose_estimator.h"
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
void TestRelativePoseEstimator(const Vector3d points_3d[5],
                                  const Matrix3d& expected_rotation,
                                  const Vector3d& expected_translation,
                                  const double tolerance) {
  InitRandomGenerator();

  const Vector3d gt_position =
      -expected_rotation.transpose() * expected_translation;

  // Calculates the image points in both views.
  std::vector<FeatureCorrespondence> features(5);
  for (int i = 0; i < 5; ++i) {
    const Vector3d proj_3d =
        expected_rotation * points_3d[i] + expected_translation;
    features[i].feature1 = points_3d[i].hnormalized();
    features[i].feature2 = proj_3d.hnormalized();
  }

  RelativePoseEstimator relative_pose_estimator;
  std::vector<RelativePose> soln_relative_poses;
  EXPECT_TRUE(
      relative_pose_estimator.EstimateModel(features, &soln_relative_poses));
  CHECK_GT(soln_relative_poses.size(), 0);

  // Among the returned solutions verify that at least one is close to the
  // expected translation and rotation.
  bool matched_transform = false;
  for (int n = 0; n < soln_relative_poses.size(); ++n) {
    const double rotation_error =
        AngleAxisd(soln_relative_poses[n].rotation.transpose() *
                   expected_rotation).angle();
    const double position_error =
        (gt_position - soln_relative_poses[n].position).norm();

    if (rotation_error < tolerance && position_error < tolerance) {
      matched_transform = true;
    }
  }
  EXPECT_TRUE(matched_transform);
}


TEST(RelativePoseEstimator, BasicTest) {
  // Ground truth.
  const Vector3d points_3d[5] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(3.0, 1.0, 2.5),
                                  Vector3d(-1.0, 1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0) };
  const Matrix3d soln_rotation =
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0)).toRotationMatrix();
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double kTolerance = 1e-4;
  TestRelativePoseEstimator(points_3d,
                            soln_rotation,
                            soln_translation.normalized(),
                            kTolerance);
}

TEST(RelativePoseEstimator, BehindCamerasTest) {
  const Camera camera1;
  RelativePoseEstimator estimator;
  for (int i = 0; i < 1000; i++) {
    Camera camera2;
    camera2.SetOrientationFromAngleAxis(Vector3d::Random());
    camera2.SetPosition(Vector3d::Random());

    RelativePose pose;
    pose.rotation = camera2.GetOrientationAsRotationMatrix();
    pose.position = camera2.GetPosition();
    const Vector3d translation = -pose.rotation * pose.position;
    pose.essential_matrix = CrossProductMatrix(translation) * pose.rotation;

    for (int j = 0; j < 100; j++) {
      const Eigen::Vector4d point3d = Eigen::Vector4d::Random();
      FeatureCorrespondence feature;
      const double depth1 = camera1.ProjectPoint(point3d, &feature.feature1);
      const double depth2 = camera2.ProjectPoint(point3d, &feature.feature2);

      const double error = estimator.Error(feature, pose);
      if (depth1 < 0 || depth2 < 0) {
        EXPECT_EQ(error, std::numeric_limits<double>::max());
      } else {
        EXPECT_NE(error, std::numeric_limits<double>::max());
      }
    }
  }
}

}  // namespace
}  // namespace theia
