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
#include <vector>

#include "gtest/gtest.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/sfm/estimators/estimate_triangulation.h"
#include "theia/util/random.h"

namespace theia {

namespace {

using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

RandomNumberGenerator rng(151);

void CreateObservations(
    const int num_observations,
    const int num_outliers,
    std::vector<Matrix3x4d>* projection_matrices,
    std::vector<Vector2d>* features) {
  const Vector4d point3d(0, 0, 8, 1.0);
  projection_matrices->resize(num_observations + num_outliers);
  features->resize(num_observations + num_outliers);

  // Initialize a random rotation and position within a 2x2x2 box around the
  // origin.
  for (int i = 0; i < num_observations; i++) {
    const Matrix3d rotation = RandomRotation(5.0, &rng);
    const Vector3d position =
        rng.RandVector3d() + Vector3d(i / 2.0, i / 2.0, i / 2.0);
    const Vector3d translation = -rotation * position;
    projection_matrices->at(i) << rotation, translation;
    features->at(i) = (projection_matrices->at(i) * point3d).hnormalized();
  }

  for (int i = 0; i < num_outliers; i++) {
    rng.SetRandom(&(*projection_matrices)[i]);
    features->at(i) = (projection_matrices->at(i) * point3d).hnormalized();
  }
}

}  // namespace

TEST(EstimateTriangulation, InsufficientObservations) {
  std::vector<Matrix3x4d> projection_matrices;
  std::vector<Vector2d> features;
  CreateObservations(1, 0, &projection_matrices, &features);

  RansacParameters params;
  params.error_thresh = 1.0 * 1.0;
  RansacSummary summary;
  Vector4d triangulated_point;
  EXPECT_FALSE(EstimateTriangulation(params,
                                     projection_matrices,
                                     features,
                                     &triangulated_point,
                                     &summary));
}

TEST(EstimateTriangulation, TwoViews) {
  std::vector<Matrix3x4d> projection_matrices;
  std::vector<Vector2d> features;
  CreateObservations(2, 0, &projection_matrices, &features);

  RansacParameters params;
  params.error_thresh = 1.0 * 1.0;
  RansacSummary summary;
  Vector4d triangulated_point;
  EXPECT_TRUE(EstimateTriangulation(params,
                                    projection_matrices,
                                    features,
                                    &triangulated_point,
                                    &summary));

  const Vector4d gt_point(0, 0, 8, 1);
  static const double kTolerance = 1e-6;
  EXPECT_LT((triangulated_point.hnormalized() - gt_point.hnormalized()).norm(),
            kTolerance);
  EXPECT_EQ(summary.inliers.size(), 2);
}

TEST(EstimateTriangulation, WithOutliers) {
  std::vector<Matrix3x4d> projection_matrices;
  std::vector<Vector2d> features;
  CreateObservations(10, 2, &projection_matrices, &features);

  RansacParameters params;
  params.error_thresh = 1.0 * 1.0;
  RansacSummary summary;
  Vector4d triangulated_point;
  EXPECT_TRUE(EstimateTriangulation(params,
                                    projection_matrices,
                                    features,
                                    &triangulated_point,
                                    &summary));

  const Vector4d gt_point(0, 0, 8, 1);
  static const double kTolerance = 1e-6;
  EXPECT_LT((triangulated_point.hnormalized() - gt_point.hnormalized()).norm(),
            kTolerance);
  EXPECT_GE(summary.inliers.size(), 6);
}

}  // namespace theia
