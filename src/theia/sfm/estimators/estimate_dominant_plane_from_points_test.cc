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
// Author: Benjamin Nuernberger (bnuernberger@cs.ucsb.edu)

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <algorithm>
#include <vector>
#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/util/random.h"
#include "theia/sfm/estimators/estimate_dominant_plane_from_points.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/test/test_utils.h"

namespace theia {
namespace {
using Eigen::Vector3d;

static const double kErrorThreshold = 1.0;

// Generate points on a plane
void GeneratePoints(std::vector<Vector3d>* points) {
  static const int kDepth = 5.0;
  for (int i = -4; i <= 4; i++) {
    for (int j = -4; j <= 4; j++) {
      // includes collinear points
      points->emplace_back(Vector3d(i, j, kDepth));
    }
  }
}

void ExecuteRandomTest(const RansacParameters& options,
                       const double inlier_ratio,
                       const double noise) {
  InitRandomGenerator();

  // Create 3D points (inliers and outliers) and add noise if appropriate.
  std::vector<Vector3d> points3d;
  GeneratePoints(&points3d);

  for (int i = 0; i < points3d.size(); i++) {
    // Add an inlier or outlier.
    if (i >= inlier_ratio * points3d.size()) {
      points3d[i] = Vector3d::Random();
    }
  }

  if (noise) {
    for (int i = 0; i < points3d.size(); i++) {
      AddNoiseToPoint(noise, &points3d[i]);
    }
  }

  // Estimate the dominant plane.
  Plane plane;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateDominantPlaneFromPoints(options,
                                 RansacType::RANSAC,
                                 points3d,
                                 &plane,
                                 &ransac_summary));

  // Expect that the inlier ratio is close to the ground truth.
  EXPECT_GT(static_cast<double>(ransac_summary.inliers.size()), 3);
}

TEST(EstimateDominantPlane, MinimalCase) {
  std::vector<Vector3d> points3d;
  points3d.emplace_back(Vector3d(1, 1, 1));
  points3d.emplace_back(Vector3d(-1, 1, 0));
  points3d.emplace_back(Vector3d(2, 0, 3));

  // Estimate the dominant plane.
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;

  Plane plane;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateDominantPlaneFromPoints(options,
                                 RansacType::RANSAC,
                                 points3d,
                                 &plane,
                                 &ransac_summary));

  EXPECT_EQ(static_cast<double>(ransac_summary.inliers.size()), 3);
}

TEST(EstimateDominantPlane, AllInliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;

  ExecuteRandomTest(options,
                    kInlierRatio,
                    kNoise);
}

TEST(EstimateDominantPlane, AllInliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;

  ExecuteRandomTest(options,
                    kInlierRatio,
                    kNoise);
}

TEST(EstimateDominantPlane, OutliersNoNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 0.7;
  const double kNoise = 0.0;

  ExecuteRandomTest(options,
                    kInlierRatio,
                    kNoise);
}

TEST(EstimateDominantPlane, OutliersWithNoise) {
  RansacParameters options;
  options.use_mle = true;
  options.error_thresh = kErrorThreshold;
  options.failure_probability = 0.001;
  const double kInlierRatio = 0.7;
  const double kNoise = 1.0;

  ExecuteRandomTest(options,
                    kInlierRatio,
                    kNoise);
}

}  // namespace
}  // namespace theia
