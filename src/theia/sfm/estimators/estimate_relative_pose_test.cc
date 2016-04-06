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
#include <Eigen/Geometry>
#include <glog/logging.h>

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/test/test_utils.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/estimate_relative_pose.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

static const double kFocalLength = 1000.0;
static const double kSampsonError = 2.0;
static const double kErrorThreshold =
    (kSampsonError * kSampsonError) / (kFocalLength * kFocalLength);

// Generate points in a grid so that they are repeatable.
void GeneratePoints(std::vector<Vector3d>* points) {
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = 4; k <= 6; k++) {
        points->emplace_back(Vector3d(i, j, k));
      }
    }
  }
}

void ExecuteRandomTest(const RansacParameters& options,
                       const Matrix3d& rotation,
                       const Vector3d& position,
                       const double inlier_ratio,
                       const double noise,
                       const double tolerance_degrees) {
  InitRandomGenerator();

  // Create feature correspondences (inliers and outliers) and add noise if
  // appropriate.
  std::vector<Vector3d> points3d;
  GeneratePoints(&points3d);

  const Vector3d translation = (-rotation * position).normalized();
  std::vector<FeatureCorrespondence> correspondences;
  for (int i = 0; i < points3d.size(); i++) {
    FeatureCorrespondence correspondence;
    // Add an inlier or outlier.
    if (i < inlier_ratio * points3d.size()) {
      // Make sure the point is in front of the camera.
      correspondence.feature1 = points3d[i].hnormalized();
      correspondence.feature2 =
          (rotation * points3d[i] + translation).hnormalized();
    } else {
      correspondence.feature1 = Vector2d::Random();
      correspondence.feature2 = Vector2d::Random();
    }
    correspondences.emplace_back(correspondence);
  }

  if (noise) {
    for (int i = 0; i < points3d.size(); i++) {
      AddNoiseToProjection(noise / kFocalLength, &correspondences[i].feature1);
      AddNoiseToProjection(noise / kFocalLength, &correspondences[i].feature2);
    }
  }

  // Estimate the relative pose.
  RelativePose relative_pose;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateRelativePose(options,
                                   RansacType::RANSAC,
                                   correspondences,
                                   &relative_pose,
                                   &ransac_summary));

  // Expect that the inlier ratio is close to the ground truth.
  EXPECT_GT(static_cast<double>(ransac_summary.inliers.size()), 5);

  Eigen::AngleAxisd rotation_loop(rotation *
                                  relative_pose.rotation.transpose());
  EXPECT_LT(RadToDeg(rotation_loop.angle()), tolerance_degrees);

  const double translation_diff_rad = std::acos(
      Clamp(position.normalized().dot(relative_pose.position), -1.0, 1.0));

  EXPECT_LT(RadToDeg(translation_diff_rad), tolerance_degrees);
}

TEST(EstimateRelativePose, AllInliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;
  const double kPoseToleranceDegrees = 1e-4;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
    AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
        .toRotationMatrix()
  };
  const std::vector<Vector3d> positions = { Vector3d(-1.3, 0, 0),
                                            Vector3d(0, 0.1, 0.5) };
  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseToleranceDegrees);
    }
  }
}

TEST(EstimateRelativePose, AllInliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;
  const double kPoseToleranceDegrees = 5.0;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
    AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
        .toRotationMatrix()
  };
  const std::vector<Vector3d> positions = { Vector3d(-1.3, 0, 0),
                                            Vector3d(0, 0.1, 0.5) };

  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseToleranceDegrees);
    }
  }
}

TEST(EstimateRelativePose, OutliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.0001;
  const double kInlierRatio = 0.7;
  const double kNoise = 0.0;
  const double kPoseToleranceDegrees = 5.0;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    ProjectToRotationMatrix(Matrix3d::Identity() + 0.2 * Matrix3d::Random())
  };
  const std::vector<Vector3d> positions = { Vector3d(1, 0.2, 0),
                                            Vector3d(0, 1, 0.1) };

  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseToleranceDegrees);
    }
  }
}

TEST(EstimateRelativePose, OutliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 0.7;
  const double kNoise = 1.0;
  const double kPoseToleranceDegrees = 5.0;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    ProjectToRotationMatrix(Matrix3d::Identity() + 0.2 * Matrix3d::Random())
  };
  const std::vector<Vector3d> positions = { Vector3d(1, 0.2, 0),
                                            Vector3d(0, 1, 0.1) };

  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseToleranceDegrees);
    }
  }
}

}  // namespace theia
