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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu), John Flynn (jflynn@google.com)

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <cmath>

#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/two_point_pose_partial_rotation.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/util/util.h"

namespace theia {

using Eigen::Map;
using Eigen::Quaterniond;
using Eigen::Vector3d;

namespace {

RandomNumberGenerator rng(58);

void TestTwoPointPoseWithNoise(
    const Vector3d points_3d[2],
    const Vector3d& axis,
    const double model_noise,
    const double projection_noise,
    const Quaterniond& test_rotation,
    const Vector3d& test_translation,
    const double max_rotation_difference,
    const double max_difference_between_translation) {
  // Sets up the image rays, since the points_3d are in the image's coordinate
  // space they are just these points after normalization.
  Vector3d image_rays[2] = { points_3d[0].normalized(),
                             points_3d[1].normalized() };

  // Sets up the model points by transforming the points_3d by the inverse
  // of the the test translation and rotation.
  Vector3d model_points[2];
  const Quaterniond inverse_test_rotation = test_rotation.inverse();
  const Vector3d inverse_test_translation =
      -1.0 * (inverse_test_rotation * test_translation);
  for (int i = 0; i < 2; ++i) {
    model_points[i] =
        inverse_test_rotation * points_3d[i] + inverse_test_translation;
  }

  // Adds projection noise if required.
  if (projection_noise) {
    for (int i = 0; i < 2; ++i) {
      AddNoiseToRay(projection_noise, &rng, &image_rays[i]);
    }
  }

  // Adds model noise if required.
  if (model_noise) {
    for (int i = 0; i < 2; ++i) {
      AddNoiseToPoint(model_noise, &rng, &model_points[i]);
    }
  }

  Quaterniond soln_rotations[2];
  Vector3d soln_translations[2];

  // Computes the pose.
  int num_solutions = TwoPointPosePartialRotation(axis,
                                                  model_points[0],
                                                  model_points[1],
                                                  image_rays[0],
                                                  image_rays[1],
                                                  soln_rotations,
                                                  soln_translations);

  // Checks that one of the solutions matches the test rotation and
  // translations.
  bool matched_transform = false;
  for (int n = 0; n < num_solutions; ++n) {
    const Quaterniond& solved_rotation = soln_rotations[n];
    const Vector3d& solved_translation = soln_translations[n];

    // Checks that the calculated quaternions are unit length.
    ASSERT_NEAR(solved_rotation.norm(), 1.0, 1e-5);

    const double angle_between_rotations =
        test_rotation.angularDistance(solved_rotation);

    bool matched_rotation = (angle_between_rotations < max_rotation_difference);

    const double translation_difference =
        (test_translation - solved_translation).norm();

    bool matched_translation =
        (translation_difference < max_difference_between_translation);

    if (matched_translation && matched_rotation) {
      matched_transform = true;
    }

    // Checks that the solution is valid, i.e. that it gives a very
    // small reprojection error.
    for (int i = 0; i < 2; ++i) {
      const Vector3d proj_3d =
          solved_rotation * model_points[i] + solved_translation;
      static const double kMaxProjectionError = 1e-2;
      for (int d = 0; d < 2; ++d) {
        EXPECT_NEAR(image_rays[i][d] / image_rays[i].z(),
                    proj_3d[d] / proj_3d.z(), kMaxProjectionError);
      }
    }
  }
  EXPECT_TRUE(matched_transform);
}

// Checks with several points, optionally adds model and projection noise.
void ManyPointsTest(const double model_noise,
                    const double projection_noise,
                    const double max_rotation_difference,
                    const double max_difference_between_translation) {
  // Sets some test rotations and translations.
  static const Vector3d kAxes[] = {
      Vector3d(0.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized()
  };

  static const double kRotationAngles[THEIA_ARRAYSIZE(kAxes)] = {
      DegToRad(7.0),
      DegToRad(12.0),
      DegToRad(15.0),
      DegToRad(20.0),
      DegToRad(11.0),
      DegToRad(0.0),  // Tests no rotation.
      DegToRad(5.0),
      DegToRad(0.0)  // Tests no rotation and no translation.
  };

  static const Vector3d kTranslations[THEIA_ARRAYSIZE(kAxes)] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 91.0),
      Vector3d(1.0, 7.0, 11.0),
      Vector3d(0.0, 0.0, 0.0),  // Tests no translation.
      Vector3d(0.0, 0.0, 0.0)  // Tests no translation and no rotation.
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
         transform_index < THEIA_ARRAYSIZE(kAxes);
         ++transform_index) {
      const Quaterniond test_rotation(Eigen::AngleAxisd(
          kRotationAngles[transform_index], kAxes[transform_index]));

      TestTwoPointPoseWithNoise(points_3d,
                                kAxes[transform_index],
                                model_noise,
                                projection_noise,
                                test_rotation,
                                kTranslations[transform_index],
                                max_rotation_difference,
                                max_difference_between_translation);
    }
  }
}

// Tests a single set of model to image correspondences and a single
// transformation consisting of a translation and a rotation around kAxis.
void BasicTest() {
  // Sets up some points in the 3D scene
  const Vector3d points_3d[2] = { Vector3d(5.0, 20.0, 23.0),
                                  Vector3d(-6.0, 16.0, 33.0) };

  const Vector3d kAxis(0.0, 1.0, 0.0);
  const Quaterniond kExpectedRotation(Eigen::AngleAxisd(0.15, kAxis));
  const Vector3d kExpectedTranslation(-3.0, 1.5, 11.0);

  static const double kModelNoise = 0.0;
  static const double kProjectionNoise = 0.0 / 512;
  static const double kMaxAllowedRotationDifference = DegToRad(1.0e-5);
  static const double kMaxAllowedTranslationDifference = 1.0e-5;
  TestTwoPointPoseWithNoise(points_3d,
                            kAxis,
                            kModelNoise,
                            kProjectionNoise,
                            kExpectedRotation,
                            kExpectedTranslation,
                            kMaxAllowedRotationDifference,
                            kMaxAllowedTranslationDifference);
}

TEST(TwoPointPoseTest, BasicTest) {
  BasicTest();
}


// Tests a single set of model to image correspondences and multiple
// transformations consisting of translation and rotations around different
// axes.
TEST(TwoPointPoseTest, DifferentAxesTest) {
  const Vector3d points_3d[2] = {
      Vector3d(5.0, 20.0, 23.0),
      Vector3d(-6.0, 16.0, 33.0)
  };

  static const Vector3d kAxes[] = {
      Vector3d(0.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized()
  };

  static const double kRotationAngles[THEIA_ARRAYSIZE(kAxes)] = {
      DegToRad(7.0),
      DegToRad(12.0),
      DegToRad(15.0),
      DegToRad(20.0),
      DegToRad(11.0),
      DegToRad(0.0),  // Tests no rotation.
      DegToRad(5.0),
      DegToRad(0.0)  // Tests no rotation and no translation.
  };

  static const Vector3d kTranslations[THEIA_ARRAYSIZE(kAxes)] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 21.0),
      Vector3d(1.0, 7.0, 11.0),
      Vector3d(0.0, 0.0, 0.0),  // Tests no translation.
      Vector3d(0.0, 0.0, 0.0)  // Tests no translation and no rotation.
  };

  for (int i = 0; i < THEIA_ARRAYSIZE(kAxes); ++i) {
    Quaterniond rotation(Eigen::AngleAxisd(kRotationAngles[i], kAxes[i]));
    static const double kModelNoise = 0.0;
    static const double kProjectionNoise = 0.0 / 512;
    static const double kMaxAllowedRotationDifference = DegToRad(1.0e-5);
    static const double kMaxAllowedTranslationDifference = 1.0e-5;
    TestTwoPointPoseWithNoise(points_3d,
                              kAxes[i],
                              kModelNoise,
                              kProjectionNoise,
                              rotation,
                              kTranslations[i],
                              kMaxAllowedRotationDifference,
                              kMaxAllowedTranslationDifference);
  }
}

// Tests many sets of model to image correspondences and multiple
// transformations consisting of translation and rotations around different
// axes, does not add any noise.
TEST(TwoPointPoseTest, ManyPointsTest) {
  static const double kModelNoise = 0.0;
  static const double kProjectionNoise = 0.0;
  static const double kMaxAllowedRotationDifference = DegToRad(1.0e-5);
  static const double kMaxAllowedTranslationDifference = 1.0e-5;

  ManyPointsTest(kModelNoise,
                 kProjectionNoise,
                 kMaxAllowedRotationDifference,
                 kMaxAllowedTranslationDifference);
}

// Tests many sets of model to image correspondences and multiple
// transformations consisting of translation and rotations around different
// axes, with noise on the image correspondeces.
TEST(TwoPointPoseTest, ManyPointsAndCorrespondenceNoiseTest) {
  static const double kModelNoise = 0.0;
  static const double kProjectionNoise = 0.5 / 512;
  static const double kMaxAllowedRotationDifference = DegToRad(10.0);
  static const double kMaxAllowedTranslationDifference = 5.0;

  ManyPointsTest(kModelNoise,
                 kProjectionNoise,
                 kMaxAllowedRotationDifference,
                 kMaxAllowedTranslationDifference);
}

// Tests many sets of model to image correspondences and multiple
// transformations consisting of translation and rotations around different
// axes, with noise on the model points.
TEST(TwoPointPoseTest, ManyPointsAndModelNoiseTest) {
  static const double kModelNoise = 0.01;
  static const double kProjectionNoise = 0.0;
  static const double kMaxAllowedRotationDifference = DegToRad(15.0);
  static const double kMaxAllowedTranslationDifference = 10.0;

  ManyPointsTest(kModelNoise,
                 kProjectionNoise,
                 kMaxAllowedRotationDifference,
                 kMaxAllowedTranslationDifference);
}

}  // namespace
}  // namespace theia
