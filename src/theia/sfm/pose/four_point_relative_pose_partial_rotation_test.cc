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
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/test/test_utils.h"
#include "theia/util/util.h"
#include "theia/util/random.h"
#include "theia/sfm/pose/four_point_relative_pose_partial_rotation.h"
#include "theia/sfm/pose/test_util.h"

namespace theia {
namespace {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

// Tests that the four point pose works correctly by taking the passed
// points_3d, projecting them to view_1_origins to get image one rays,
// transforming by (expected_rotation, expected_translation) to get
// image two rays and then verifying that the FourPointEssentialMatrix
// function returns (expected_rotation, expected_translation) among its
// solutions.
// Noise can be added to the image projections by setting projection_noise.
// The thresholds for rotation and translation similarity can be controlled
// by max_rotation_difference and max_difference_between_translation.
// If projection_noise_std_dev is non-zero then a random noise generator must
// be passed in rng.
void TestFourPointResultWithNoise(const Vector3d& axis,
                                  const Vector3d points_3d[4],
                                  const Vector3d view_1_origins[4],
                                  const Vector3d view_2_origins[4],
                                  double projection_noise_std_dev,
                                  const Quaterniond& expected_rotation,
                                  const Vector3d& expected_translation,
                                  double max_rotation_difference,
                                  double max_translation_difference) {
  Vector3d image_one_rays[4];
  Vector3d image_two_rays[4];

  for (int i = 0; i < 4; ++i) {
    image_one_rays[i] = (points_3d[i] - view_1_origins[i]).normalized();
    const Vector3d transformed_point =
        expected_rotation * points_3d[i] + expected_translation;
    image_two_rays[i]= (transformed_point - view_2_origins[i]).normalized();
  }
  if (projection_noise_std_dev) {
    for (int i = 0; i < 4; ++i) {
      AddNoiseToRay(projection_noise_std_dev, &image_one_rays[i]);
      AddNoiseToRay(projection_noise_std_dev, &image_two_rays[i]);
    }
  }
  std::vector<Quaterniond> soln_rotations;
  std::vector<Vector3d> soln_translations;

  FourPointRelativePosePartialRotation(axis,
                                       image_one_rays,
                                       view_1_origins,
                                       image_two_rays,
                                       view_2_origins,
                                       &soln_rotations,
                                       &soln_translations);

  bool matched_transform = false;
  for (int n = 0; n < soln_rotations.size(); ++n) {
    const double rotation_difference =
        expected_rotation.angularDistance(soln_rotations[n]);

    bool matched_rotation = (rotation_difference < max_rotation_difference);

    const double translation_difference =
        (expected_translation - soln_translations[n]).norm();

    bool matched_translation =
        (translation_difference < max_translation_difference);

    if (matched_translation && matched_rotation) {
      matched_transform = true;
      break;
    }
  }
  EXPECT_TRUE(matched_transform);
}

void BasicTest() {
  // Sets up some points in the 3D scene
  const Vector3d kPoints3D[4] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0),
                                  Vector3d(4.0, 3.0, 5.0) };

  const Vector3d kImageOneOrigins[4] = { Vector3d(-1.0, 0.0, 0.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(2.0, 0.0, 0.0),
                                         Vector3d(3.0, 0.0, 0.0) };

  const Vector3d kImageTwoOrigins[4] = { Vector3d(0.0, 1.0, 0.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(0.0, 2.0, 0.0),
                                         Vector3d(0.0, 3.0, 0.0) };

  const Vector3d axis = Vector3d(1.0, 1.0, 1.0).normalized();
  Quaterniond kExpectedRotation(AngleAxisd(DegToRad(51.0), axis));

  const Vector3d kExpectedTranslation(-2.0, 3.0, -5.0);
  double kProjectionNoise = 0.0;
  double kMaxAllowedRotationDifference = 1e-5;
  double kMaxAllowedTranslationDifference = 1e-4;
  TestFourPointResultWithNoise(axis,
                               kPoints3D,
                               kImageOneOrigins,
                               kImageTwoOrigins,
                               kProjectionNoise,
                               kExpectedRotation,
                               kExpectedTranslation,
                               kMaxAllowedRotationDifference,
                               kMaxAllowedTranslationDifference);
}

TEST(FourPointEssentialMatrixTest, Basic) {
  BasicTest();
}

TEST(FourPointEssentialMatrixTest, NoRotation) {
  // Sets up some points in the 3D scene
  const Vector3d kPoints3D[4] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 2.0, 5.0),
                                  Vector3d(4.0, 3.0, 5.0) };

  const Vector3d kImageOneOrigins[4] = { Vector3d(-1.0, 0.0, -1.5),
                                         Vector3d(0.0, 0.0, -1.0),
                                         Vector3d(2.0, 0.0, 0.0),
                                         Vector3d(3.0, 0.0, 0.0) };

  const Vector3d kImageTwoOrigins[4] = { Vector3d(0.0, 1.0, 2.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(0.0, 2.0, 1.0),
                                         Vector3d(0.0, 3.0, 0.0) };

  const Vector3d axis = Vector3d(1.0, 1.0, 1.0).normalized();
  Quaterniond kExpectedRotation(AngleAxisd(DegToRad(0.0), axis));

  const Vector3d kExpectedTranslation(1.0, 1.0, 1.0);
  double kProjectionNoise = 0.0;
  double kMaxAllowedRotationDifference = 1e-5;
  double kMaxAllowedTranslationDifference = 1e-4;
  TestFourPointResultWithNoise(axis,
                               kPoints3D,
                               kImageOneOrigins,
                               kImageTwoOrigins,
                               kProjectionNoise,
                               kExpectedRotation,
                               kExpectedTranslation,
                               kMaxAllowedRotationDifference,
                               kMaxAllowedTranslationDifference);
}

// Tests a variety of axes, angles and translations with added projection noise.
TEST(FourPointEssentialMatrixTest, NoiseTest) {
  const Vector3d kPoints3D[4] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0),
                                  Vector3d(4.0, 3.0, 5.0) };

  const Vector3d kImageOneOrigins[4] = { Vector3d(-1.0, 0.0, 0.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(2.0, 0.0, 0.0),
                                         Vector3d(3.0, 0.0, 0.0) };

  const Vector3d kImageTwoOrigins[4] = { Vector3d(0.0, 1.0, 0.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(0.0, 2.0, 0.0),
                                         Vector3d(0.0, 3.0, 0.0) };

  const Vector3d kAxes[5] = { Vector3d(0.0, 0.0, 1.0).normalized(),
                              Vector3d(1.0, 0.0, 0.0).normalized(),
                              Vector3d(1.0, 0.0, 1.0).normalized(),
                              Vector3d(1.0, 1.0, 0.0).normalized(),
                              Vector3d(1.0, 1.0, 1.0).normalized() };

  const double kAngles[5] = {11.0, 5.0, 2.0, 13.0, 12.0};

  const Vector3d kTranslations[5] = {
    Vector3d(1.0, 1.0, 1.0), Vector3d(4.0, 5.0, 11.0), Vector3d(1.0, 2.0, 15.0),
    Vector3d(6.0, 3.0, 2.0), Vector3d(13.0, 1.0, 15.0)
  };

  InitRandomGenerator();

  for (int transform_index = 0;
       transform_index < THEIA_ARRAYSIZE(kAxes);
       ++transform_index) {
    Quaterniond kExpectedRotation(
        AngleAxisd(DegToRad(kAngles[transform_index]), kAxes[transform_index]));

    const double kProjectionNoise = 0.5 / 512;
    const double kMaxAllowedRotationDifference = DegToRad(4.0);
    const double kMaxAllowedTranslationDifference = 2.0;

    TestFourPointResultWithNoise(kAxes[transform_index],
                                 kPoints3D,
                                 kImageOneOrigins,
                                 kImageTwoOrigins,
                                 kProjectionNoise,
                                 kExpectedRotation,
                                 kTranslations[transform_index],
                                 kMaxAllowedRotationDifference,
                                 kMaxAllowedTranslationDifference);
  }
}

// Tests that the solver degrades gracefully when the passed axis is not exactly
// correct.
TEST(FourPointEssentialMatrixTest, IncorrectAxisTest) {
  const Vector3d kPoints3D[4] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0),
                                  Vector3d(4.0, 3.0, 5.0) };

  const Vector3d kImageOneOrigins[4] = { Vector3d(-1.0, 0.0, 0.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(2.0, 0.0, 0.0),
                                         Vector3d(3.0, 0.0, 0.0) };

  const Vector3d kImageTwoOrigins[4] = { Vector3d(0.0, 1.0, 0.0),
                                         Vector3d(0.0, 0.0, 0.0),
                                         Vector3d(0.0, 2.0, 0.0),
                                         Vector3d(0.0, 3.0, 0.0) };

  const Vector3d kAxes[5] = {
      Vector3d(0.0, 0.0, 1.0).normalized(),
      Vector3d(1.0, 0.0, 0.0).normalized(),
      Vector3d(1.0, 0.0, 1.0).normalized(),
      Vector3d(1.0, 1.0, 0.0).normalized(),
      Vector3d(1.0, 1.0, 1.0).normalized()
  };

  const double kAngles[5] = {11.0, 5.0, 2.0, 13.0, 12.0};

  const Vector3d kTranslations[5] = {
      Vector3d(1.0, 1.0, 1.0),
      Vector3d(4.0, 5.0, 11.0),
      Vector3d(1.0, 2.0, 15.0),
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

      // Tests the ThreePointEssentialMatrix function with fairly large
      // thresholds on the allowed translation and rotation differences.
      const double kProjectionNoise = 0.0;
      const double kMaxAllowedRotationDifference = DegToRad(2.0);
      const double kMaxAllowedTranslationDifference = 2.0;

      TestFourPointResultWithNoise(kAxes[transform_index],
                                   kPoints3D,
                                   kImageOneOrigins,
                                   kImageTwoOrigins,
                                   kProjectionNoise,
                                   kExpectedRotation,
                                   kTranslations[transform_index],
                                   kMaxAllowedRotationDifference,
                                   kMaxAllowedTranslationDifference);
    }
  }
}

}  // namespace
}  // namespace theia
