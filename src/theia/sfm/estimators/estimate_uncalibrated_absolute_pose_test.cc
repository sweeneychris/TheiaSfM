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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/test/test_utils.h"
#include "theia/sfm/estimators/estimate_uncalibrated_absolute_pose.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/util/random.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

static const int kNumPoints = 100;
static const double kFocalLength = 1000.0;
static const double kReprojectionError = 4.0;
static const double kErrorThreshold = kReprojectionError * kReprojectionError;

void ExecuteRandomTest(const RansacParameters& options,
                       const Matrix3d& rotation,
                       const Vector3d& position,
                       const double inlier_ratio,
                       const double noise,
                       const double tolerance) {
  InitRandomGenerator();

  // Create feature correspondences (inliers and outliers) and add noise if
  // appropriate.
  std::vector<FeatureCorrespondence2D3D> correspondences;
  for (int i = 0; i < kNumPoints; i++) {
    FeatureCorrespondence2D3D correspondence;
    correspondence.world_point = Vector3d(
        RandDouble(-2.0, 2.0), RandDouble(-2.0, 2.0), RandDouble(6.0, 10.0));

    // Add an inlier or outlier.
    if (i < inlier_ratio * kNumPoints) {
      // Make sure the point is in front of the camera.
      correspondence.feature =
          kFocalLength *
          ((rotation * (correspondence.world_point - position)).hnormalized());
    } else {
      correspondence.feature = kFocalLength * Vector2d::Random();
    }
    correspondences.emplace_back(correspondence);
  }

  if (noise) {
    for (int i = 0; i < kNumPoints; i++) {
      AddNoiseToProjection(noise, &correspondences[i].feature);
    }
  }

  // Estimate the absolute pose.
  UncalibratedAbsolutePose pose;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateUncalibratedAbsolutePose(options,
                                             RansacType::RANSAC,
                                             correspondences,
                                             &pose,
                                             &ransac_summary));

  // Expect that the inlier ratio is close to the ground truth.
  EXPECT_GT(static_cast<double>(ransac_summary.inliers.size()), 3);

  // Expect poses are near.
  EXPECT_TRUE(test::ArraysEqualUpToScale(9,
                                         rotation.data(),
                                         pose.rotation.data(),
                                         tolerance));
  // The position is more noisy than the rotation usually.
  EXPECT_TRUE(test::ArraysEqualUpToScale(3,
                                         position.data(),
                                         pose.position.data(),
                                         2.0 * tolerance));

  // Expect focal length is near.
  static const double kFocalLengthTolerance = 0.05;
  EXPECT_NEAR(pose.focal_length,
              kFocalLength,
              kFocalLengthTolerance * kFocalLength);
}

TEST(EstimateUncalibratedAbsolutePose, AllInliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;
  const double kPoseTolerance = 1e-4;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
    AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
        .toRotationMatrix()
  };
  const std::vector<Vector3d> positions = { Vector3d(-1.3, 0, 0),
                                            Vector3d(0, 0, 0.5) };
  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseTolerance);
    }
  }
}

TEST(EstimateUncalibratedAbsolutePose, AllInliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;
  const double kPoseTolerance = 1e-2;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
    AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
        .toRotationMatrix()
  };
  const std::vector<Vector3d> positions = { Vector3d(-1.3, 0, 0),
                                            Vector3d(0, 0, 0.5) };

  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseTolerance);
    }
  }
}

TEST(EstimateUncalibratedAbsolutePose, OutliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.7;
  const double kNoise = 0.0;
  const double kPoseTolerance = 1e-2;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    ProjectToRotationMatrix(Matrix3d::Identity() + 0.3 * Matrix3d::Random())
  };
  const std::vector<Vector3d> positions = { Vector3d(1, 0, 0),
                                            Vector3d(0, 1, 0) };

  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseTolerance);
    }
  }
}

TEST(EstimateUncalibratedAbsolutePose, OutliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.7;
  const double kNoise = 1.0;
  const double kPoseTolerance = 0.1;

  const std::vector<Matrix3d> rotations = {
    Matrix3d::Identity(),
    ProjectToRotationMatrix(Matrix3d::Identity() + 0.3 * Matrix3d::Random())
  };
  const std::vector<Vector3d> positions = { Vector3d(1, 0, 0),
                                            Vector3d(0, 1, 0) };

  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options,
                        rotations[i],
                        positions[j],
                        kInlierRatio,
                        kNoise,
                        kPoseTolerance);
    }
  }
}

}  // namespace theia
