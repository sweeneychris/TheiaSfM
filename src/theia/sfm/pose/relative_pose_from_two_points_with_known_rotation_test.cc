// Copyright (C) 2013 The Regents of the University of California (Regents)
// and Google, Inc. All rights reserved.
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
//     * Neither the name of The Regents or University of California, Google,
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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
#include <cmath>

#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/sfm/pose/relative_pose_from_two_points_with_known_rotation.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/util/util.h"

namespace theia {

using Eigen::Map;
using Eigen::Quaterniond;
using Eigen::Vector2d;
using Eigen::Vector3d;

namespace {

RandomNumberGenerator rng(58);

void TestRelativePositionFromTwoRaysWithNoise(
    const Vector3d points_3d[2],
    const double model_noise,
    const double projection_noise,
    const Vector3d& test_position,
    const double max_difference_between_position) {
  static const double kTolerance = 1e-4;

  // Sets up the model points by transforming the points_3d by the inverse
  // of the the test position.
  Vector2d rotated_rays1[2];
  Vector2d rotated_rays2[2];
  for (int i = 0; i < 2; ++i) {
    rotated_rays1[i] = points_3d[i].hnormalized();
    rotated_rays2[i] = (points_3d[i] - test_position).hnormalized();
  }

  // Adds projection noise if required.
  if (projection_noise) {
    for (int i = 0; i < 2; ++i) {
      AddNoiseToProjection(projection_noise, &rng, &rotated_rays1[i]);
      AddNoiseToProjection(projection_noise, &rng, &rotated_rays2[i]);
    }
  }

  // Computes the pose.
  Vector3d soln_position;
  EXPECT_TRUE(RelativePoseFromTwoPointsWithKnownRotation(
      rotated_rays1, rotated_rays2, &soln_position));
  const Eigen::Vector3d relative_position = test_position.normalized();
  EXPECT_LT((soln_position - relative_position).norm(), kTolerance)
      << "Expected position: " << relative_position.transpose()
      << " vs estimated: " << soln_position.transpose();
}

// Checks with several points, optionally adds model and projection noise.
void ManyPointsTest(const double model_noise,
                    const double projection_noise,
                    const double max_difference_between_position) {
  static const Vector3d kPositions[8] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 91.0),
      Vector3d(1.0, 7.0, 11.0),
      Vector3d(1.0, 0.0, 0.0),
      Vector3d(0.0, 0.0, 1.0)
  };

  // Sets up some test points.
  static const double kTestPoints[][6] = {
      {-1.62, -2.99, 6.12, 4.42, -1.53, 9.83},
      {1.45, -0.59, 5.29, 1.89, -1.10, 8.22},
      {-0.21, 2.38, 5.63, 0.61, -0.97, 7.49},
      {0.48, 0.70, 8.94, 1.65, -2.56, 8.63},
      {2.44, -0.20, 7.78, 2.84, -2.58, 7.35},
      {-1.35, -2.84, 7.33, -0.42, 1.54, 8.86},
      {2.56, 1.72, 7.86, 1.75, -1.39, 5.73},
      {2.08, -3.91, 8.37, -0.91, 1.36, 9.16},
      {2.84, 1.54, 8.74, -1.01, 3.02, 8.18},
      {-3.73, -0.62, 7.81, -2.98, -1.88, 6.23},
      {2.39, -0.19, 6.47, -0.63, -1.05, 7.11},
      {-1.76, -0.55, 5.18, -3.19, 3.27, 8.18},
      {0.31, -2.77, 7.54, 0.54, -3.77, 9.77},
  };

  for (int i = 0; i < THEIA_ARRAYSIZE(kTestPoints); ++i) {
    const Vector3d points_3d[2] = {
      Vector3d(kTestPoints[i][0], kTestPoints[i][1], kTestPoints[i][2]),
      Vector3d(kTestPoints[i][3], kTestPoints[i][4], kTestPoints[i][5]),
    };

    for (int transform_index = 0;
         transform_index < THEIA_ARRAYSIZE(kPositions);
         ++transform_index) {
      TestRelativePositionFromTwoRaysWithNoise(points_3d,
                                               model_noise,
                                               projection_noise,
                                               kPositions[transform_index],
                                               max_difference_between_position);
    }
  }
}

// Tests a single set of model to image correspondences and a single
// transformation consisting of a position and a rotation around kAxis.
void BasicTest() {
  // Sets up some points in the 3D scene
  const Vector3d points_3d[2] = { Vector3d(5.0, 20.0, 23.0),
                                  Vector3d(-6.0, 16.0, 33.0) };

  const Vector3d kAxis(0.0, 1.0, 0.0);
  const Vector3d kExpectedPosition(-3.0, 1.5, 11.0);

  static const double kModelNoise = 0.0;
  static const double kProjectionNoise = 0.0 / 512;
  static const double kMaxAllowedPositionDifference = 1.0e-5;
  TestRelativePositionFromTwoRaysWithNoise(points_3d,
                                           kModelNoise,
                                           kProjectionNoise,
                                           kExpectedPosition,
                                           kMaxAllowedPositionDifference);
}

TEST(PositionFromTwoRaysTest, BasicTest) {
  BasicTest();
}

// Tests a single set of model to image correspondences and multiple
// transformations consisting of position and rotations around different
// axes.
TEST(PositionFromTwoRaysTest, DifferentAxesTest) {
  const Vector3d points_3d[2] = {
      Vector3d(5.0, 20.0, 23.0),
      Vector3d(-6.0, 16.0, 33.0)
  };

  static const Vector3d kPositions[8] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 21.0),
      Vector3d(1.0, 7.0, 11.0),
      Vector3d(1.0, 0.0, 0.0),
      Vector3d(0.0, 0.0, 1.0)
  };

  for (int i = 0; i < THEIA_ARRAYSIZE(kPositions); ++i) {
    static const double kModelNoise = 0.0;
    static const double kProjectionNoise = 0.0 / 512;

    static const double kMaxAllowedPositionDifference = 1.0e-5;
    TestRelativePositionFromTwoRaysWithNoise(points_3d,
                                             kModelNoise,
                                             kProjectionNoise,
                                             kPositions[i],
                                             kMaxAllowedPositionDifference);
  }
}

// Tests many sets of model to image correspondences and multiple
// transformations consisting of position and rotations around different
// axes, does not add any noise.
TEST(PositionFromTwoRaysTest, ManyPointsTest) {
  static const double kModelNoise = 0.0;
  static const double kProjectionNoise = 0.0;
  static const double kMaxAllowedPositionDifference = 1.0e-5;

  ManyPointsTest(kModelNoise, kProjectionNoise, kMaxAllowedPositionDifference);
}

}  // namespace
}  // namespace theia
