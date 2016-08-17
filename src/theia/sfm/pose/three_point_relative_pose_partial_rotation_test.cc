// Copyright (C) 2014 The Regents of the University of California (Regents)
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
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include <limits>
#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/test/test_utils.h"
#include "theia/util/util.h"
#include "theia/util/random.h"
#include "theia/sfm/pose/three_point_relative_pose_partial_rotation.h"
#include "theia/sfm/pose/test_util.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Quaterniond;

RandomNumberGenerator rng(60);

// Tests that the three point pose works correctly by taking the passed
// points_3d, projecting them to get image one rays, transforming by
// (expected_rotation, expected_translation) to get image two rays and then
// verifying that the ThreePointRelativePosePartialRotation function returns
// (test_rotation, test_translation_direction) among its solutions.
// Noise can be added to the image projections by setting projection_noise.
// The thresholds for rotation and translation similarity can be controlled
// by max_rotation_difference and max_difference_between_translation.
// If projection_noise_std_dev is non-zero then a random noise generator must
// be passed in rng.
void TestThreePointResultWithNoise(const Vector3d& axis,
                                   const Vector3d points_3d[3],
                                   const double projection_noise_std_dev,
                                   const Quaterniond& expected_rotation,
                                   const Vector3d& expected_translation,
                                   const double max_rotation_difference,
                                   const double max_angle_between_translation) {
  // Calculates the image rays in both views.
  Vector3d view_one_rays[3];
  Vector3d view_two_rays[3];
  for (int i = 0; i < 3; ++i) {
    const Vector3d proj_3d =
        expected_rotation * points_3d[i] + expected_translation;
    view_one_rays[i] = points_3d[i].normalized();
    view_two_rays[i] = proj_3d.normalized();
  }
  if (projection_noise_std_dev) {
    // Adds noise to both of the rays.
    for (int i = 0; i < 3; ++i) {
      Eigen::Vector2d view_one_ray = view_one_rays[i].hnormalized();
      AddNoiseToProjection(projection_noise_std_dev, &rng, &view_one_ray);
      view_one_rays[i] = view_one_ray.homogeneous().normalized();

      Eigen::Vector2d view_two_ray = view_two_rays[i].hnormalized();
      AddNoiseToProjection(projection_noise_std_dev, &rng, &view_two_ray);
      view_two_rays[i] = view_two_ray.homogeneous().normalized();
    }
  }

  // Calculates the essential matrix, this may return multiple solutions.
  std::vector<Quaterniond> soln_rotations;
  std::vector<Vector3d> soln_translations;
  ThreePointRelativePosePartialRotation(axis,
                                        view_one_rays,
                                        view_two_rays,
                                        &soln_rotations,
                                        &soln_translations);
  EXPECT_GT(soln_rotations.size(), 0);

  // Among the returned solutions verify that at least one is close to the
  // expected translation and rotation.
  bool matched_transform = false;
  for (int n = 0; n < soln_rotations.size(); ++n) {
    double rotation_difference =
        expected_rotation.angularDistance(soln_rotations[n]);

    bool matched_rotation = (rotation_difference < max_rotation_difference);

    // The translation is only known up to scale so this verifies that the
    // translations have matching directions.
    const double translation_angle_difference =
        AngleAxisd(Quaterniond::FromTwoVectors(soln_translations[n],
                                               expected_translation)).angle();

    bool matched_translation =
        (translation_angle_difference < max_angle_between_translation);

    if (matched_translation && matched_rotation) {
      matched_transform = true;
    }
  }
  CHECK(matched_transform);
}

// A basic test with a few points and no noise.
void BasicTest() {
  const Vector3d points_3d[3] = {Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0)};
  const Vector3d axis = Vector3d(0.0, 0.0, 1.0).normalized();
  Quaterniond kExpectedRotation(AngleAxisd(DegToRad(51.0), axis));

  const Vector3d kExpectedTranslation(1.0, 1.0, 1.0);

  const double kProjectionNoise = 0.0;
  const double kMaxAllowedRotationDifference = 1e-5;
  const double kMaxAllowedAngleBetweenTranslations = 1e-8;

  TestThreePointResultWithNoise(axis,
                                points_3d,
                                kProjectionNoise,
                                kExpectedRotation,
                                kExpectedTranslation,
                                kMaxAllowedRotationDifference,
                                kMaxAllowedAngleBetweenTranslations);
}

TEST(ThreePointRelativePosePartialRotationTest, Basic) {
  BasicTest();
}

// A basic test with a few points, no noise and no rotation.
TEST(ThreePointRelativePosePartialRotationTest, NoRotationTest) {
  const Vector3d points_3d[3] = {Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0)};
  // This test fails with these points.
  // const Vector3d points_3d[3] = {Vector3d(-1.0, 3.0, -3.0),
  //                                 Vector3d(1.0, -1.0, 2.0),
  //                                 Vector3d(2.0, 1.0, 9.0)};
  const Vector3d axis = Vector3d(0.0, 0.0, 1.0).normalized();
  Quaterniond kExpectedRotation(AngleAxisd(0.0, axis));

  const Vector3d kExpectedTranslation(1.0, 1.0, 1.0);

  double kProjectionNoise = 0.0;
  double kMaxAllowedRotationDifference = 1e-5;
  double kMaxAllowedAngleBetweenTranslations = 1e-4;
  TestThreePointResultWithNoise(axis,
                                points_3d,
                                kProjectionNoise,
                                kExpectedRotation,
                                kExpectedTranslation,
                                kMaxAllowedRotationDifference,
                                kMaxAllowedAngleBetweenTranslations);
}

// Tests a variety of axes, angles and translations with added projection noise.
TEST(ThreePointRelativePosePartialRotationTest, NoiseTest) {
  const Vector3d kPoints3D[3] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0) };

  const Vector3d kAxes[7] = {
      Vector3d(0.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized()
  };

  const double kAngles[7] = { 11.0, 7.0, 5.0, 2.0, 3.0, 13.0, 12.0 };

  const Vector3d kTranslations[7] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 91.0),
      Vector3d(6.0, 3.0, 2.0),
      Vector3d(13.0, 1.0, 15.0)
  };

  for (int transform_index = 0; transform_index < THEIA_ARRAYSIZE(kAxes);
       ++transform_index) {
    Quaterniond kExpectedRotation(
        AngleAxisd(DegToRad(kAngles[transform_index]), kAxes[transform_index]));
    kExpectedRotation.normalize();

    const double kProjectionNoise = 1.0 / 512;
    const double kMaxAllowedRotationDifference = DegToRad(8.0);
    const double kMaxAllowedAngleBetweenTranslations = DegToRad(2.0);
    TestThreePointResultWithNoise(kAxes[transform_index],
                                  kPoints3D,
                                  kProjectionNoise,
                                  kExpectedRotation,
                                  kTranslations[transform_index],
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedAngleBetweenTranslations);
  }
}

// Tests that the solver degrade gracefully when the passed axis is not exactly
// correct.
TEST(ThreePointRelativePosePartialRotationTest, IncorrectAxisTest) {
  const Vector3d kPoints3D[3] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0) };

  const Vector3d kAxes[7] = {
      Vector3d(0.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized()
  };

  const double kAngles[7] = { 11.0, 7.0, 5.0, 2.0, 3.0, 13.0, 12.0 };

  const Vector3d kTranslations[7] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 91.0),
      Vector3d(6.0, 3.0, 2.0),
      Vector3d(13.0, 1.0, 15.0)
  };

  // The axes are perturbed by these rotation to simulate the axis not being
  // perfectly known.
  const Quaterniond kAxisPerturbations[6] = {
    Quaterniond(AngleAxisd(DegToRad(1.0), Vector3d(1.0, 0, 0))),
    Quaterniond(AngleAxisd(DegToRad(-1.0), Vector3d(1.0, 0, 0))),
    Quaterniond(AngleAxisd(DegToRad(1.0), Vector3d(0.0, 1.0, 0))),
    Quaterniond(AngleAxisd(DegToRad(-1.0), Vector3d(0.0, 1.0, 0))),
    Quaterniond(AngleAxisd(DegToRad(1.0), Vector3d(0.0, 0.0, 1.0))),
    Quaterniond(AngleAxisd(DegToRad(-1.0), Vector3d(0.0, 0.0, 1.0)))
  };

  for (int transform_index = 0;
       transform_index < THEIA_ARRAYSIZE(kAxes);
       ++transform_index) {
    for (int axis_rotation_index = 0;
         axis_rotation_index < THEIA_ARRAYSIZE(kAxisPerturbations);
         ++axis_rotation_index) {
      // Perturbs the axis by the axis perturbation.
      const Vector3d perturbed_axis =
          kAxisPerturbations[axis_rotation_index] * kAxes[transform_index];
      // Uses the perturbed rotation as the expected rotation, but the original
      // axis is passed as the axis.
      Quaterniond kExpectedRotation(
          AngleAxisd(DegToRad(kAngles[transform_index]), perturbed_axis));

      // Tests the ThreePointRelativePosePartialRotation function with fairly
      // large thresholds on the allowed translation and rotation differences.
      const double kProjectionNoise = 0.0;
      const double kMaxAllowedRotationDifference =
          DegToRad(2.0);
      const double kMaxAllowedAngleBetweenTranslations =
          DegToRad(2.0);

      TestThreePointResultWithNoise(kAxes[transform_index],
                                    kPoints3D,
                                    kProjectionNoise,
                                    kExpectedRotation,
                                    kTranslations[transform_index],
                                    kMaxAllowedRotationDifference,
                                    kMaxAllowedAngleBetweenTranslations);
    }
  }
}

// Comprehensive test with many points and different axes and translations.
TEST(ThreePointRelativePosePartialRotationTest, ManyPoints) {
  const Vector3d kAxes[7] = {
      Vector3d(0.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 1.0).normalized(),
      Vector3d(0.0, 1.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized()
  };
  const double kAngles[7] = {11.0, 7.0, 5.0, 2.0, 3.0, 13.0, 12.0};

  const Vector3d kTranslations[7] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(3.0, 2.0, 13.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
      Vector3d(3.0, 1.5, 91.0),
      Vector3d(6.0, 3.0, 2.0),
      Vector3d(13.0, 1.0, 15.0)
  };

  for (int i = 0; i < 1000; ++i) {
    for (int transform_index = 0; transform_index < THEIA_ARRAYSIZE(kAxes);
         ++transform_index) {
      std::vector<Vector3d> random_points;
      CreateRandomPointsInFrustum(1.0,
                                  1.0,
                                  2.0,
                                  10.0,
                                  3,
                                  &rng,
                                  &random_points);

      const Vector3d points_3d[3] = { random_points[0],
                                      random_points[1],
                                      random_points[2] };
      Quaterniond kExpectedRotation(AngleAxisd(
          DegToRad(kAngles[transform_index]), kAxes[transform_index]));

      const double kProjectionNoise = 0.0;
      const double kMaxAllowedRotationDifference = 1e-5;
      const double kMaxAllowedAngleBetweenTranslations = 1e-4;
      TestThreePointResultWithNoise(kAxes[transform_index],
                                    points_3d,
                                    kProjectionNoise,
                                    kExpectedRotation,
                                    kTranslations[transform_index],
                                    kMaxAllowedRotationDifference,
                                    kMaxAllowedAngleBetweenTranslations);
    }
  }
}

}  // namespace theia
