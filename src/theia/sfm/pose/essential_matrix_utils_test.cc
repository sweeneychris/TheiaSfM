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

#include "gtest/gtest.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/essential_matrix_utils.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/util/random.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

RandomNumberGenerator rng(51);

TEST(DecomposeEssentialMatrix, BasicTest) {
  const double kTranslationTolerance = 1e-6;
  const double kRotationTolerance = 1e-4;

  for (int i = 0; i < 100; i++) {
    const Matrix3d gt_rotation = RandomRotation(10.0, &rng);
    const Vector3d gt_translation = rng.RandVector3d().normalized();
    const Matrix3d essential_matrix =
        CrossProductMatrix(gt_translation) * gt_rotation;

    Matrix3d rotation1, rotation2;
    Vector3d translation;
    DecomposeEssentialMatrix(essential_matrix,
                             &rotation1,
                             &rotation2,
                             &translation);

    const double translation_dist =
        std::min((translation - gt_translation).norm(),
                 (translation + gt_translation).norm());

    const Eigen::AngleAxisd rotation1_aa(gt_rotation.transpose() * rotation1);
    const Eigen::AngleAxisd rotation2_aa(gt_rotation.transpose() * rotation2);
    const double rotation1_dist = rotation1_aa.angle();
    const double rotation2_dist = rotation2_aa.angle();

    EXPECT_TRUE(translation_dist < kTranslationTolerance &&
                (rotation1_dist < kRotationTolerance ||
                 rotation2_dist < kRotationTolerance));
  }
}

void TestGetBestPoseFromEssentialMatrix(const int num_inliers,
                                        const int num_outliers) {
  static const double kTolerance = 1e-12;

  for (int i = 0; i < 100; i++) {
    const Matrix3d gt_rotation = RandomRotation(15.0, &rng);

    const Vector3d gt_translation = rng.RandVector3d().normalized();
    const Vector3d gt_position = -gt_rotation.transpose() * gt_translation;
    const Matrix3d essential_matrix =
        CrossProductMatrix(gt_translation) * gt_rotation;

    // Create Correspondences.
    std::vector<FeatureCorrespondence> correspondences;
    for (int j = 0; j < num_inliers; j++) {
      // Make sure the point is in front of the camera.
      const Vector3d point_3d = rng.RandVector3d() + Vector3d(0, 0, 100);
      const Vector3d proj_3d = gt_rotation * point_3d + gt_translation;

      FeatureCorrespondence correspondence;
      correspondence.feature1 = point_3d.hnormalized();
      correspondence.feature2 = proj_3d.hnormalized();
      correspondences.emplace_back(correspondence);
    }

    // Add outliers
    for (int j = 0; j < num_outliers; j++) {
      // Make sure the point is in front of the camera.
      const Vector3d point_3d = rng.RandVector3d() + Vector3d(0, 0, -100);
      const Vector3d proj_3d = gt_rotation * point_3d + gt_translation;

      FeatureCorrespondence correspondence;
      correspondence.feature1 = point_3d.hnormalized();
      correspondence.feature2 = proj_3d.hnormalized();
      correspondences.emplace_back(correspondence);
    }

    Matrix3d estimated_rotation;
    Vector3d estimated_position;
    const int num_points_in_front = GetBestPoseFromEssentialMatrix(
        essential_matrix,
        correspondences,
        &estimated_rotation,
        &estimated_position);

    // Ensure that the results are correct. Sincer there is no noise we can
    // expect te number of point in front to be exact.
    EXPECT_EQ(num_points_in_front, num_inliers);
    EXPECT_LT((gt_rotation - estimated_rotation).norm(), kTolerance);
    EXPECT_LT((gt_position - estimated_position).norm(), kTolerance);
  }
}

TEST(GetBestPoseFromEssentialMatrix, AllInliers) {
  static const int kNumInliers = 100;
  static const int kNumOutliers = 0;
  TestGetBestPoseFromEssentialMatrix(kNumInliers, kNumOutliers);
}

TEST(GetBestPoseFromEssentialMatrix, MostlyInliers) {
  static const int kNumInliers = 100;
  static const int kNumOutliers = 50;
  TestGetBestPoseFromEssentialMatrix(kNumInliers, kNumOutliers);
}

}  // namespace theia
