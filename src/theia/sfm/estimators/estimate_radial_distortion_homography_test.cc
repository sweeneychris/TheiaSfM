// Copyright (C) 2019 The Regents of the University of California (Regents).
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

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// May 2019

#include <glog/logging.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

#include "theia/matching/feature_correspondence.h"
#include "theia/math/util.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/estimate_radial_distortion_homography.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/test/test_utils.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

RandomNumberGenerator rng(53);

static constexpr double kFocalLength1 = 1000.0;
static constexpr double kFocalLength2 = 1500.0;
static constexpr double kRadialDistortion1 = -1e-7;
static constexpr double kRadialDistortion2 = -5e-7;
static constexpr double kReprojectionError = 3.0;

// 10 percent error
static constexpr double kRadDistThreshold1 =
    0.1 * std::abs(kRadialDistortion1 * std::pow(kFocalLength1, 2));
static constexpr double kRadDistThreshold2 =
    0.1 * std::abs(kRadialDistortion2 * std::pow(kFocalLength2, 2));

// Generate points on a plane so that a homography can accurately estimate the
// motion.
void GeneratePoints(std::vector<Vector3d>* points) {
  static const int kDepth = 3.0;
  for (int i = -4; i <= 4; i++) {
    for (int j = -4; j <= 4; j++) {
      points->emplace_back(Vector3d(i, j, kDepth));
    }
  }
}

void GenerateDistortedImagePoint(
    const Vector3d& point_3d, const double projection_noise_std_dev,
    const Matrix3d& expected_rotation, const Vector3d& expected_translation,
    const double focal_length1, const double focal_length2,
    const double radial_distortion1, const double radial_distortion2,
    const bool outlier, Vector3d& image_1_point_normalized,
    Vector3d& image_2_point_normalized, Vector2d& image_1_point,
    Vector2d& image_2_point) {
  Vector3d point3_cam2 = expected_rotation * point_3d + expected_translation;

  DistortPoint(point_3d, focal_length1, radial_distortion1, image_1_point);
  DistortPoint(point3_cam2, focal_length2, radial_distortion2, image_2_point);

  if (outlier) {
    image_1_point[0] += rng.RandDouble(-100.0, 100.0);
    image_2_point[0] += rng.RandDouble(-100.0, 100.0);
    image_1_point[1] += rng.RandDouble(-100.0, 100.0);
    image_2_point[1] += rng.RandDouble(-100.0, 100.0);
  }

  if (projection_noise_std_dev) {
    AddNoiseToProjection(projection_noise_std_dev, &rng, &image_1_point);
    AddNoiseToProjection(projection_noise_std_dev, &rng, &image_2_point);
  }

  // normalize points with focal length (and principal point) estimate for
  // estimation
  UndistortPoint(image_1_point, focal_length1, 0.0, image_1_point_normalized);
  UndistortPoint(image_2_point, focal_length2, 0.0, image_2_point_normalized);
}

// Creates a test scenario from ground truth 3D points and ground truth rotation
// and translation. Projection noise is optional (set to 0 for no
// noise).
void GenerateDistortedImagePoints(
    const std::vector<Vector3d>& points_3d,
    const double projection_noise_std_dev, const Matrix3d& expected_rotation,
    const Vector3d& expected_translation, const double focal_length1,
    const double focal_length2, const double radial_distortion1,
    const double radial_distortion2, const double inlier_ratio,
    std::vector<Vector2d>* image_1_points_normalized,
    std::vector<Vector2d>* image_2_points_normalized,
    std::vector<Vector2d>* image_1_points,
    std::vector<Vector2d>* image_2_points) {
  for (int i = 0; i < points_3d.size(); i++) {
    Vector2d distorted_point1, distorted_point2;
    Vector3d normalized_points2d_1, normalized_points2d_2;
    bool outlier = false;
    if (i > inlier_ratio * points_3d.size()) {
      outlier = true;
    }
    GenerateDistortedImagePoint(
        points_3d[i], projection_noise_std_dev, expected_rotation,
        expected_translation, focal_length1, focal_length2, radial_distortion1,
        radial_distortion2, outlier, normalized_points2d_1,
        normalized_points2d_2, distorted_point1, distorted_point2);

    image_1_points->push_back(distorted_point1);
    image_2_points->push_back(distorted_point2);

    image_1_points_normalized->push_back(normalized_points2d_1.hnormalized());
    image_2_points_normalized->push_back(normalized_points2d_2.hnormalized());
  }
}

void ExecuteRandomTest(const RansacParameters& options,
                       const Matrix3d& rotation, const Vector3d& position,
                       const double radial_distortion_1,
                       const double radial_distortion_2,
                       const double focal_length_1, const double focal_length_2,
                       const double inlier_ratio, const double noise) {
  // Create feature correspondences (inliers and outliers) and add noise if
  // appropriate.
  std::vector<Vector3d> points3d;
  GeneratePoints(&points3d);

  std::vector<Vector2d> image_1_points_normalized, image_2_points_normalized;
  std::vector<Vector2d> image_1_points, image_2_points;

  // Generate distorted points
  GenerateDistortedImagePoints(
      points3d, noise, rotation, position, focal_length_1, focal_length_2,
      radial_distortion_1, radial_distortion_2, inlier_ratio,
      &image_1_points_normalized, &image_2_points_normalized, &image_1_points,
      &image_2_points);
  // Get correspondence
  std::vector<RadialDistortionFeatureCorrespondence> correspondences;
  for (size_t i = 0; i < image_1_points.size(); ++i) {
    RadialDistortionFeatureCorrespondence correspondence;
    correspondence.feature_left = image_1_points[i];
    correspondence.feature_right = image_2_points[i];
    correspondence.normalized_feature_left = image_1_points_normalized[i];
    correspondence.normalized_feature_right = image_2_points_normalized[i];
    correspondence.focal_length_estimate_left = focal_length_1;
    correspondence.focal_length_estimate_right = focal_length_2;
    correspondences.push_back(correspondence);
  }
  // Estimate the radial distortion homography.
  RadialHomographyResult radial_homography_result;
  RansacSummary ransac_summary;
  EXPECT_TRUE(EstimateRadialHomographyMatrix(
      options, RansacType::RANSAC, correspondences, &radial_homography_result,
      &ransac_summary));
  // Expect that the radial distortion values are close to the ground truth
  EXPECT_NEAR(radial_homography_result.l1,
              radial_distortion_1 * (focal_length_1 * focal_length_1),
              kRadDistThreshold1);
  EXPECT_NEAR(radial_homography_result.l2,
              radial_distortion_2 * (focal_length_2 * focal_length_2),
              kRadDistThreshold2);
}

TEST(EstimateRadialHomographyMatrix, AllInliersNoNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(rng);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.01;
  options.max_iterations = 1000;
  options.min_iterations = 1;
  const double kInlierRatio = 1.0;
  const double kNoise = 0.0;

  const std::vector<Matrix3d> rotations = {
      Matrix3d::Identity(),
      AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
      AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
          .toRotationMatrix()};
  const std::vector<Vector3d> positions = {Vector3d(-1.3, 0, 0),
                                           Vector3d(0, 1.3, 0.0)};
  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options, rotations[i], positions[j], kRadialDistortion1,
                        kRadialDistortion2, kFocalLength1, kFocalLength2,
                        kInlierRatio, kNoise);
    }
  }
}

TEST(EstimateRadialHomographyMatrix, AllInliersNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(rng);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.01;
  options.max_iterations = 1000;
  options.min_iterations = 1;
  const double kInlierRatio = 1.0;
  const double kNoise = 1.0;

  const std::vector<Matrix3d> rotations = {
      Matrix3d::Identity(),
      AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
      AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
          .toRotationMatrix()};
  const std::vector<Vector3d> positions = {Vector3d(-1.3, 0, 0),
                                           Vector3d(0, 1.3, 0.0)};
  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options, rotations[i], positions[j], kRadialDistortion1,
                        kRadialDistortion2, kFocalLength1, kFocalLength2,
                        kInlierRatio, kNoise);
    }
  }
}

TEST(EstimateRadialHomographyMatrix, OutliersNoNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(rng);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.01;
  options.max_iterations = 1000;
  options.min_iterations = 1;
  const double kInlierRatio = 0.6;
  const double kNoise = 0.0;

  const std::vector<Matrix3d> rotations = {
      Matrix3d::Identity(),
      AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
      AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
          .toRotationMatrix()};
  const std::vector<Vector3d> positions = {Vector3d(-1.3, 0, 0),
                                           Vector3d(0, 1.3, 0.0)};
  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options, rotations[i], positions[j], kRadialDistortion1,
                        kRadialDistortion2, kFocalLength1, kFocalLength2,
                        kInlierRatio, kNoise);
    }
  }
}

TEST(EstimateRadialHomographyMatrix, OutliersNoise) {
  RansacParameters options;
  options.rng = std::make_shared<RandomNumberGenerator>(rng);
  options.use_mle = true;
  options.error_thresh = kReprojectionError;
  options.failure_probability = 0.01;
  options.max_iterations = 1000;
  options.min_iterations = 1;
  const double kInlierRatio = 0.6;
  const double kNoise = 1.0;

  const std::vector<Matrix3d> rotations = {
      Matrix3d::Identity(),
      AngleAxisd(DegToRad(12.0), Vector3d::UnitY()).toRotationMatrix(),
      AngleAxisd(DegToRad(-9.0), Vector3d(1.0, 0.2, -0.8).normalized())
          .toRotationMatrix()};
  const std::vector<Vector3d> positions = {Vector3d(-1.3, 0, 0),
                                           Vector3d(0, 1.3, 0.0)};
  for (int i = 0; i < rotations.size(); i++) {
    for (int j = 0; j < positions.size(); j++) {
      ExecuteRandomTest(options, rotations[i], positions[j], kRadialDistortion1,
                        kRadialDistortion2, kFocalLength1, kFocalLength2,
                        kInlierRatio, kNoise);
    }
  }
}

}  // theia
