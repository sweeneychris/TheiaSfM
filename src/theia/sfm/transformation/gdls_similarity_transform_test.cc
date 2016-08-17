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

#include "theia/alignment/alignment.h"
#include "theia/math/util.h"
#include "theia/util/random.h"
#include "theia/util/util.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/transformation/gdls_similarity_transform.h"

namespace theia {
namespace {
using Eigen::AngleAxisd;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

RandomNumberGenerator rng(57);

void TestGdlsSimilarityTransformWithNoise(
    const std::vector<Vector3d>& camera_centers,
    const std::vector<Vector3d>& world_points,
    const double projection_noise_std_dev,
    const Quaterniond& expected_rotation,
    const Vector3d& expected_translation,
    const double expected_scale,
    const double max_reprojection_error,
    const double max_rotation_difference,
    const double max_translation_difference,
    const double max_scale_difference) {
  const int num_points = world_points.size();
  const int num_cameras = camera_centers.size();

  std::vector<Vector3d> camera_rays;
  std::vector<Vector3d> ray_origins;
  camera_rays.reserve(num_points);
  ray_origins.reserve(num_points);
  for (int i = 0; i < num_points; i++) {
    Vector3d ray_origin = (expected_rotation * camera_centers[i % num_cameras] +
                           expected_translation) / expected_scale;

    ray_origins.push_back(ray_origin);

    // Reproject 3D points into camera frame.
    camera_rays.push_back(
        (expected_rotation * world_points[i] + expected_translation -
         expected_scale * ray_origins[i]).normalized());
  }

  if (projection_noise_std_dev) {
    // Adds noise to both of the rays.
    for (int i = 0; i < num_points; i++) {
      AddNoiseToRay(projection_noise_std_dev, &rng, &camera_rays[i]);
    }
  }

  // Run DLS Similarity Transform.
  std::vector<Quaterniond> soln_rotation;
  std::vector<Vector3d> soln_translation;
  std::vector<double> soln_scale;
  GdlsSimilarityTransform(ray_origins, camera_rays, world_points,
                          &soln_rotation, &soln_translation, &soln_scale);

  // Check solutions and verify at least one is close to the actual solution.

  const int num_solutions = soln_rotation.size();
  EXPECT_GT(num_solutions, 0);
  bool matched_transform = false;
  for (int i = 0; i < num_solutions; i++) {
    // Check that reprojection errors are small.
    double max_reproj_err = 0.0;
    for (int j = 0; j < num_points; j++) {
      const Quaterniond unrot =
          Quaterniond::FromTwoVectors(camera_rays[j], Vector3d(0, 0, 1));
      const Vector3d reprojected_point =
          (soln_rotation[i] * world_points[j] + soln_translation[i]) /
              soln_scale[i] - ray_origins[j];

      const Vector3d unrot_cam_ray = unrot * camera_rays[j];
      const Vector3d unrot_reproj_pt = unrot * reprojected_point;

      const double reprojection_error =
          (unrot_cam_ray.hnormalized() - unrot_reproj_pt.hnormalized()).norm();
      if (reprojection_error > max_reproj_err)
        max_reproj_err = reprojection_error;
      EXPECT_LE(reprojection_error, max_reprojection_error)
          << "Reproj error is " << reprojection_error * 512.0;
    }

    // Check that the solution is accurate.
    const double rotation_difference =
        expected_rotation.angularDistance(soln_rotation[i]);
    const bool matched_rotation =
        (rotation_difference < max_rotation_difference);
    const double translation_difference =
        (expected_translation - soln_translation[i]).squaredNorm();
    const bool matched_translation =
        (translation_difference < max_translation_difference);
    const double scale_difference =
        fabs(expected_scale - soln_scale[i]) / expected_scale;
    const bool matched_scale = (scale_difference < max_scale_difference);

    if (matched_translation && matched_rotation && matched_scale) {
      matched_transform = true;
    }
  }
  EXPECT_TRUE(matched_transform);
}

void BasicTest() {
  const std::vector<Vector3d> points_3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double soln_scale = 2.5;
  const double kNoise = 0.0;
  const double kMaxReprojectionError = 1.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1e-4);
  const double kMaxAllowedTranslationDifference = 1e-6;
  const double kMaxAllowedScaleDifference = 1e-6;

  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };

  TestGdlsSimilarityTransformWithNoise(
      kImageOrigins,
      points_3d,
      kNoise,
      soln_rotation,
      soln_translation,
      soln_scale,
      kMaxReprojectionError,
      kMaxAllowedRotationDifference,
      kMaxAllowedTranslationDifference,
      kMaxAllowedScaleDifference);
}

TEST(GdlsSimilarityTransform, Basic) {
  BasicTest();
}

TEST(GdlsSimilarityTransform, NoiseTest) {
  const std::vector<Vector3d> points_3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0),
                                            Vector3d(-1.0, -3.0, 2.0),
                                            Vector3d(1.0, -2.0, 1.0),
                                            Vector3d(-1.0, 4.0, 2.0),
                                            Vector3d(-2.0, 2.0, 3.0)
  };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(0.0, 1.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(0.0, 2.0, 0.0),
                                                Vector3d(0.0, 3.0, 0.0) };

  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double soln_scale = 2.5;
  const double kNoise = 1.0 / 512.0;
  const double kMaxReprojectionError = 5.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1.0);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const double kMaxAllowedScaleDifference = 1e-2;

  TestGdlsSimilarityTransformWithNoise(
      kImageOrigins,
      points_3d,
      kNoise,
      soln_rotation,
      soln_translation,
      soln_scale,
      kMaxReprojectionError,
      kMaxAllowedRotationDifference,
      kMaxAllowedTranslationDifference,
      kMaxAllowedScaleDifference);
}

TEST(GdlsSimilarityTransform, ManyPoints) {
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
    Vector3d(3.0, 1.5, 18.0),
    Vector3d(1.0, 7.0, 11.0),
    Vector3d(0.0, 0.0, 0.0),  // Tests no translation.
    Vector3d(0.0, 0.0, 0.0)  // Tests no translation and no rotation.
  };

  static const double kScales[THEIA_ARRAYSIZE(kAxes)] = {
    0.33,
    1.0,
    4.2,
    10.13,
    3.14,
    7.22,
    0.1,
    10.0
  };

  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };

  static const int num_points[3] = { 100, 500, 1000 };
  const double kNoise = 0.5 / 512.0;
  const double kMaxReprojectionError = 10.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(0.3);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const double kMaxAllowedScaleDifference = 1e-2;

  for (int i = 0; i < THEIA_ARRAYSIZE(kAxes); i++) {
    const Quaterniond soln_rotation(AngleAxisd(kRotationAngles[i], kAxes[i]));
    for (int j = 0; j < THEIA_ARRAYSIZE(num_points); j++) {
      std::vector<Vector3d> points_3d;
      points_3d.reserve(num_points[j]);
      for (int k = 0; k < num_points[j]; k++) {
        points_3d.push_back(Vector3d(RandDouble(-5.0, 5.0),
                                     RandDouble(-5.0, 5.0),
                                     RandDouble(2.0, 10.0)));
      }

      TestGdlsSimilarityTransformWithNoise(
          kImageOrigins,
          points_3d,
          kNoise,
          soln_rotation,
          kTranslations[i],
          kScales[i],
          kMaxReprojectionError,
          kMaxAllowedRotationDifference,
          kMaxAllowedTranslationDifference,
          kMaxAllowedScaleDifference);
    }
  }
}

TEST(GdlsSimilarityTransform, NoRotation) {
  const std::vector<Vector3d> points_3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0),
                                            Vector3d(-1.0, -3.0, 2.0),
                                            Vector3d(1.0, -2.0, 1.0),
                                            Vector3d(-1.0, 4.0, 2.0),
                                            Vector3d(-2.0, 2.0, 3.0)
  };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };

  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(0.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double soln_scale = 0.77;
  const double kNoise = 0.5 / 512.0;
  const double kMaxReprojectionError = 2.0 / 5.0;
  const double kMaxAllowedRotationDifference = DegToRad(0.2);
  const double kMaxAllowedTranslationDifference = 5e-4;
  const double kMaxAllowedScaleDifference = 1e-2;

  TestGdlsSimilarityTransformWithNoise(
      kImageOrigins,
      points_3d,
      kNoise,
      soln_rotation,
      soln_translation,
      soln_scale,
      kMaxReprojectionError,
      kMaxAllowedRotationDifference,
      kMaxAllowedTranslationDifference,
      kMaxAllowedScaleDifference);
}

TEST(GdlsSimilarityTransform, NoTranslation) {
  const std::vector<Vector3d> points_3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0),
                                            Vector3d(-1.0, -3.0, 2.0),
                                            Vector3d(1.0, -2.0, 1.0),
                                            Vector3d(-1.0, 4.0, 2.0),
                                            Vector3d(-2.0, 2.0, 3.0)
  };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };

  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(0.0, 0.0, 0.0);
  const double soln_scale = 2.5;
  const double kNoise = 1.0 / 512.0;
  const double kMaxReprojectionError = 4.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(0.5);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const double kMaxAllowedScaleDifference = 1e-2;

  TestGdlsSimilarityTransformWithNoise(
      kImageOrigins,
      points_3d,
      kNoise,
      soln_rotation,
      soln_translation,
      soln_scale,
      kMaxReprojectionError,
      kMaxAllowedRotationDifference,
      kMaxAllowedTranslationDifference,
      kMaxAllowedScaleDifference);
}

}  // namespace
}  // namespace theia
