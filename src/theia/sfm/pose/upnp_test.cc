// Copyright (C) 2018 The Regents of the University of California (Regents).
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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include "gtest/gtest.h"

#include <vector>

#include "theia/math/util.h"
#include "theia/sfm/pose/upnp.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/util/random.h"

namespace theia {
namespace {

using Eigen::Quaterniond;
using Eigen::AngleAxisd;
using Eigen::Vector3d;

RandomNumberGenerator rng(57);

struct InputDatum {
  std::vector<Eigen::Vector3d> ray_origins;
  std::vector<Eigen::Vector3d> ray_directions;
  std::vector<Eigen::Vector3d> world_points;
};

InputDatum ComputeInputDatum(
    const std::vector<Vector3d>& world_points,
    const std::vector<Vector3d>& camera_centers,
    const Quaterniond& expected_rotation,
    const Vector3d& expected_translation) {
  const int num_points = world_points.size();
  const int num_cameras = camera_centers.size();
  
  InputDatum input_datum;
  std::vector<Vector3d>& ray_directions = input_datum.ray_directions;
  std::vector<Vector3d>& ray_origins = input_datum.ray_origins;
  ray_directions.reserve(num_points);
  ray_origins.reserve(num_points);

  for (int i = 0; i < num_points; ++i) {
    // Ray origin wrt to coordinate system of the camera rig.
    const Vector3d ray_origin =
        expected_rotation * camera_centers[i % num_cameras] +
        expected_translation;

    ray_origins.emplace_back(std::move(ray_origin));

    // Reproject 3D points into camera frame.
    ray_directions.emplace_back(
        (expected_rotation * world_points[i] + expected_translation -
         ray_origins[i]).normalized());
  }

  input_datum.world_points = world_points;
  return input_datum;
}

bool CheckReprojectionErrors(const InputDatum& input_datum,
                             const Eigen::Quaterniond& soln_rotation,
                             const Eigen::Vector3d& soln_translation,
                             const double kMaxReprojectionError) {
  const int num_points = input_datum.world_points.size();
  double good_reprojection_errors = true;
  for (int i = 0; i < num_points; ++i) {
    const Quaterniond unrot =
        Quaterniond::FromTwoVectors(input_datum.ray_directions[i],
                                    Vector3d(0, 0, 1));
    const Vector3d reprojected_point =
        soln_rotation * input_datum.world_points[i] + soln_translation -
        input_datum.ray_origins[i];
  
    const Vector3d unrot_cam_ray = unrot * input_datum.ray_directions[i];
    const Vector3d unrot_reproj_pt = unrot * reprojected_point;

    const double reprojection_error =
        (unrot_cam_ray.hnormalized() - unrot_reproj_pt.hnormalized()).norm();

    good_reprojection_errors = (good_reprojection_errors &&
                                (reprojection_error < kMaxReprojectionError));
  }
  return good_reprojection_errors;
}

void TestUpnpPoseEstimationWithNoise(
    const Quaterniond& expected_rotation,
    const Vector3d& expected_translation,
    const double projection_noise_std_dev,
    const double max_reprojection_error,
    const double max_rotation_difference,
    const double max_translation_difference,
    InputDatum* input_datum) {
  const int num_points = input_datum->world_points.size();
  // Add noise to ray.
  if (projection_noise_std_dev > 0.0) {
    for (int i = 0; i < input_datum->world_points.size(); ++i) {
      AddNoiseToRay(projection_noise_std_dev,
                    &rng,
                    &input_datum->ray_directions[i]);
    }
  }

  // Estimate pose.
  std::vector<Quaterniond> solution_rotations;
  std::vector<Vector3d> solution_translations;
  const UpnpCostParameters upnp_params =
      Upnp(input_datum->ray_origins,
           input_datum->ray_directions,
           input_datum->world_points,
           &solution_rotations,
           &solution_translations);
  const double upnp_cost = EvaluateUpnpCost(upnp_params, expected_rotation);
  VLOG(2) << "Upnp cost with expected rotation: " << upnp_cost;

  // Check solutions and verify at least one is close to the actual solution.
  const int num_solutions = solution_rotations.size();
  EXPECT_GT(num_solutions, 0);
  bool matched_transform = false;
  for (int i = 0; i < num_solutions; ++i) {
    const bool good_reprojection_errors =
        CheckReprojectionErrors(*input_datum,
                                solution_rotations[i],
                                solution_translations[i],
                                max_reprojection_error);
    const double rotation_difference =
        expected_rotation.angularDistance(solution_rotations[i]);
    const bool matched_rotation = rotation_difference < max_rotation_difference;
    const double translation_difference =
        (expected_translation - solution_translations[i]).squaredNorm();
    const bool matched_translation =
        translation_difference < max_translation_difference;
    VLOG(2) << "Matched rotation: " << matched_rotation
            << " rotation error [deg]=" << RadToDeg(rotation_difference);
    VLOG(2) << "Matched translation: " << matched_translation
            << " translation error=" << translation_difference;
    VLOG(2) << "Good reprojection errors: " << good_reprojection_errors;
    if (matched_rotation && matched_translation && good_reprojection_errors) {
      matched_transform = true;
    }
  }
  EXPECT_TRUE(matched_transform);
}

// Verifies that the cost-function parameters are correct for central cameras.
TEST(UpnpTests, ComputeCostParametersForCentralCameraPoseEstimation) {
  const std::vector<Vector3d> kPoints3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const std::vector<Vector3d> kImageOrigin = { Vector3d(2.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double kNoise = 0.0;
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigin,
                                             soln_rotation,
                                             soln_translation);

  std::vector<Eigen::Quaterniond> solution_rotations;
  std::vector<Eigen::Vector3d> solution_translations;
  const UpnpCostParameters upnp_params = Upnp(input_datum.ray_origins,
                                              input_datum.ray_directions,
                                              input_datum.world_points,
                                              &solution_rotations,
                                              &solution_translations);

  const double upnp_cost = EvaluateUpnpCost(upnp_params, soln_rotation);
  EXPECT_NEAR(upnp_cost, 0.0, 1e-6);
}

// Verifies that the cost-function parameters are correct for non-central
// cameras.
TEST(UpnpTests, ComputeCostParametersForNonCentralCameraPoseEstimation) {
  const std::vector<Vector3d> kPoints3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double kNoise = 0.0;
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);

  std::vector<Eigen::Quaterniond> solution_rotations;
  std::vector<Eigen::Vector3d> solution_translations;
  const UpnpCostParameters upnp_params = Upnp(input_datum.ray_origins,
                                              input_datum.ray_directions,
                                              input_datum.world_points,
                                              &solution_rotations,
                                              &solution_translations);

  const double upnp_cost = EvaluateUpnpCost(upnp_params, soln_rotation);
  EXPECT_NEAR(upnp_cost, 0.0, 1e-6);
}

// Checks the case of a minimal sample and central camera pose estimation.
TEST(UpnpTests, MinimalSampleCentralCameraPoseEstimation) {
  const double kNoiseStdDev = 0.0;
  const double kMaxReprojectionError = 1.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1e-4);
  const double kMaxAllowedTranslationDifference = 1e-6;
  const std::vector<Vector3d> kPoints3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const std::vector<Vector3d> kImageOrigin = { Vector3d(2.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigin,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks the case of a minimal sample and central camera pose estimation.
TEST(UpnpTests, MinimalSampleCentralCameraPoseEstimationWithNoise) {
  const double kNoiseStdDev = 0.6 / 512.0;
  const double kMaxReprojectionError = 1.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1.0);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const std::vector<Vector3d> kPoints3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const std::vector<Vector3d> kImageOrigin = { Vector3d(2.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigin,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks the case of a minimal sample and a non-central camera pose estimation.
TEST(UpnpTests, MinimalSampleNonCentralCameraPoseEstimation) {
  const double kNoiseStdDev = 0.0;
  const double kMaxReprojectionError = 1.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1e-4);
  const double kMaxAllowedTranslationDifference = 1e-6;
  const std::vector<Vector3d> kPoints3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks the case of a minimal sample and a non-central camera pose estimation
// with noise.
TEST(UpnpTests, MinimalSampleNonCentralCameraPoseEstimationWithNoise) {
  const double kNoiseStdDev = 1e-3;
  const double kMaxReprojectionError = 3.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1.0);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const std::vector<Vector3d> kPoints3d = { Vector3d(-1.0, 3.0, 3.0),
                                            Vector3d(1.0, -1.0, 2.0),
                                            Vector3d(-1.0, 1.0, 2.0),
                                            Vector3d(2.0, 1.0, 3.0) };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(-1.0, 0.0, 0.0),
                                                Vector3d(0.0, 0.0, 0.0),
                                                Vector3d(2.0, 0.0, 0.0),
                                                Vector3d(3.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks that the estimator works well using a non-minimal sample of data
// points.
TEST(UpnpTests, NonMinimalSampleCentralCameraPoseEstimation) {
  const double kNoiseStdDev = 0.0;
  const double kMaxReprojectionError = 1.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1e-4);
  const double kMaxAllowedTranslationDifference = 1e-6;
  const std::vector<Vector3d> kPoints3d = {
    Vector3d(-1.0, 3.0, 3.0),
    Vector3d(1.0, -1.0, 2.0),
    Vector3d(-1.0, 1.0, 2.0),
    Vector3d(2.0, 1.0, 3.0),
    Vector3d(-1.0, -3.0, 2.0),
    Vector3d(1.0, -2.0, 1.0),
    Vector3d(-1.0, 4.0, 2.0),
    Vector3d(-2.0, 2.0, 3.0)
  };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(0.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks that the estimator works well using a non-minimal sample of data
// points.
TEST(UpnpTests, NonMinimalSampleCentralCameraPoseEstimationWithNoise) {
  const double kNoiseStdDev = 1e-3;
  const double kMaxReprojectionError = 3.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1.0);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const std::vector<Vector3d> kPoints3d = {
    Vector3d(-1.0, 3.0, 3.0),
    Vector3d(1.0, -1.0, 2.0),
    Vector3d(-1.0, 1.0, 2.0),
    Vector3d(2.0, 1.0, 3.0),
    Vector3d(-1.0, -3.0, 2.0),
    Vector3d(1.0, -2.0, 1.0),
    Vector3d(-1.0, 4.0, 2.0),
    Vector3d(-2.0, 2.0, 3.0)
  };
  const std::vector<Vector3d> kImageOrigins = { Vector3d(0.0, 0.0, 0.0) };
  const Quaterniond soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks that the estimator works well using a non-minimal sample of data.
TEST(UpnpTests, NonMinimalSampleNonCentralCameraPoseEstimation) {
  const double kNoiseStdDev = 0.0;
  const double kMaxReprojectionError = 1.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1e-4);
  const double kMaxAllowedTranslationDifference = 1e-6;
  const std::vector<Vector3d> kPoints3d = {
    Vector3d(-1.0, 3.0, 3.0),
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
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

// Checks that the estimator works well using a non-minimal sample of data
// points with noise.
TEST(UpnpTests, NonMinimalSampleNonCentralCameraPoseEstimationWithNoise) {
  const double kNoiseStdDev = 1e-3;
  const double kMaxReprojectionError = 3.0 / 512.0;
  const double kMaxAllowedRotationDifference = DegToRad(1.0);
  const double kMaxAllowedTranslationDifference = 1e-3;
  const std::vector<Vector3d> kPoints3d = {
    Vector3d(-1.0, 3.0, 3.0),
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
      AngleAxisd(DegToRad(13.0), Vector3d(1.0, 0.0, 1.0).normalized()));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  // Compute input datum.
  InputDatum input_datum = ComputeInputDatum(kPoints3d,
                                             kImageOrigins,
                                             soln_rotation,
                                             soln_translation);
  // Execute test.
  TestUpnpPoseEstimationWithNoise(soln_rotation,
                                  soln_translation,
                                  kNoiseStdDev,
                                  kMaxReprojectionError,
                                  kMaxAllowedRotationDifference,
                                  kMaxAllowedTranslationDifference,
                                  &input_datum);
}

}  // namespace
}  // namespace theia
