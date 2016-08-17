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
#include "theia/sfm/estimators/estimate_uncalibrated_relative_pose.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

static const int kNumTrials = 100;
RandomNumberGenerator rng(60);

void ExecuteRandomTest(const RansacParameters& options,
                       const Matrix3d& rotation,
                       const Vector3d& position,
                       const double focal_length1,
                       const double focal_length2,
                       const double inlier_ratio,
                       const double noise) {
  static const int kNumCorrespondences = 600;
  InitRandomGenerator();

  // Create feature correspondences (inliers and outliers) and add noise if
  // appropriate.
  const Vector3d translation = (-rotation * position).normalized();
  std::vector<FeatureCorrespondence> correspondences;
  for (int i = 0; i < kNumCorrespondences; i++) {
    FeatureCorrespondence correspondence;
    // Add an inlier or outlier.
    if (i < inlier_ratio * kNumCorrespondences) {
      // Make sure the point is in front of the camera.
      const Vector3d point_3d = rng.RandVector3d() + Vector3d(0, 0, 4);
      correspondence.feature1 = focal_length1 * point_3d.hnormalized();
      correspondence.feature2 =
          focal_length2 * (rotation * point_3d + translation).hnormalized();

      AddNoiseToProjection(noise, &rng, &correspondence.feature1);
      AddNoiseToProjection(noise, &rng, &correspondence.feature2);
    } else {
      correspondence.feature1 =
          focal_length1 *
          Vector2d(rng.RandDouble(-1.0, 1.0), rng.RandDouble(-1.0, 1.0));
      correspondence.feature2 =
          focal_length2 *
          Vector2d(rng.RandDouble(-1.0, 1.0), rng.RandDouble(-1.0, 1.0));
    }
    correspondences.emplace_back(correspondence);
  }

  // Estimate the relative pose.
  UncalibratedRelativePose relative_pose;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateUncalibratedRelativePose(options,
                                               RansacType::RANSAC,
                                               correspondences,
                                               &relative_pose,
                                               &ransac_summary));

  // Expect that the inlier ratio is close to the ground truth.
  const double observed_inlier_ratio =
      ransac_summary.inliers.size() /
      static_cast<double>(correspondences.size());
  EXPECT_GT(observed_inlier_ratio, 0.7 * inlier_ratio);
}

TEST(EstimateUncalibratedRelativePose, AllInliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = 2;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;

  for (int k = 0; k < kNumTrials; k++) {
    const Matrix3d rotation = RandomRotation(10.0, &rng);
    const Vector3d position = rng.RandVector3d();
    const double focal_length1 = rng.RandDouble(800, 1600);
    const double focal_length2 = rng.RandDouble(800, 1600);
    ExecuteRandomTest(options,
                      rotation,
                      position,
                      focal_length1,
                      focal_length2,
                      kInlierRatio,
                      kNoise);
  }
}

TEST(EstimateUncalibratedRelativePose, AllInliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = 2;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;

  for (int k = 0; k < kNumTrials; k++) {
const Matrix3d rotation = RandomRotation(10.0, &rng);
const Vector3d position = rng.RandVector3d();
    const double focal_length1 = rng.RandDouble(800, 1600);
    const double focal_length2 = rng.RandDouble(800, 1600);
    ExecuteRandomTest(options,
                      rotation,
                      position,
                      focal_length1,
                      focal_length2,
                      kInlierRatio,
                      kNoise);
  }
}

TEST(EstimateUncalibratedRelativePose, OutliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = 2;
  options.failure_probability = 0.001;
  const double kInlierRatio = 0.7;
  const double kNoise = 0.0;

  for (int k = 0; k < kNumTrials; k++) {
    const Matrix3d rotation = RandomRotation(10.0, &rng);
    const Vector3d position = rng.RandVector3d();
    const double focal_length1 = rng.RandDouble(800, 1600);
    const double focal_length2 = rng.RandDouble(800, 1600);
    ExecuteRandomTest(options,
                      rotation,
                      position,
                      focal_length1,
                      focal_length2,
                      kInlierRatio,
                      kNoise);
  }
}

TEST(EstimateUncalibratedRelativePose, OutliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.failure_probability = 0.001;
  options.error_thresh = 4.0 * 4.0;
  options.max_iterations = 1000;
  const double kInlierRatio = 0.7;
  const double kNoise = 1.0;

  for (int k = 0; k < kNumTrials; k++) {
    const Matrix3d rotation = RandomRotation(10.0, &rng);
    const Vector3d position = rng.RandVector3d();
    const double focal_length1 = rng.RandDouble(800, 1600);
    const double focal_length2 = rng.RandDouble(800, 1600);

    ExecuteRandomTest(options,
                      rotation,
                      position,
                      focal_length1,
                      focal_length2,
                      kInlierRatio,
                      kNoise);
  }
}

}  // namespace theia
