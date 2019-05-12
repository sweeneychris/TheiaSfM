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
using Eigen::Quaterniond;
using Eigen::Vector2d;
using Eigen::Vector3d;

RandomNumberGenerator rng(53);

static const int kNumTrials = 100;

void GenerateDistortedImagePoint(
    const Vector3d& point_3d, const double projection_noise_std_dev,
    const Quaterniond& expected_rotation, const Vector3d& expected_translation,
    const double focal_length1, const double focal_length2,
    const double radial_distortion1, const double radial_distortion2,
    const bool outlier, Vector3d& image_1_point_normalized,
    Vector3d& image_2_point_normalized, Vector2d& image_1_point,
    Vector2d& image_2_point) {

  Vector3d point3_cam2 =
      (expected_rotation * point_3d + expected_translation);

  Vector2d distorted_point1, distorted_point2;
  DistortPoint(point_3d, focal_length1, radial_distortion1, image_1_point);
  DistortPoint(point3_cam2, focal_length2, radial_distortion2, image_2_point);

  if (outlier) {
    image_1_point[0] += rng.RandDouble(-100.0, 100.0);
    image_2_point[0] += rng.RandDouble(-100.0, 100.0);
    image_1_point[1] += rng.RandDouble(-100.0, 100.0);
    image_2_point[1] += rng.RandDouble(-100.0, 100.0);
  }

  if (projection_noise_std_dev) {
      AddNoiseToProjection(projection_noise_std_dev, &rng,
                           &image_1_point);
      AddNoiseToProjection(projection_noise_std_dev, &rng,
                           &image_2_point);
  }

  // normalize points with focal length (and principal point) estimate for
  // estimation
  UndistortPoint(distorted_point1, focal_length1, 0.0,
                 image_1_point_normalized);
  UndistortPoint(distorted_point2, focal_length2, 0.0,
                 image_2_point_normalized);
}

// Creates a test scenario from ground truth 3D points and ground truth rotation
// and translation. Projection noise is optional (set to 0 for no
// noise).
void GenerateDistortedImagePoints(
    const std::vector<Vector3d>& points_3d,
    const double projection_noise_std_dev, const Quaterniond& expected_rotation,
    const Vector3d& expected_translation, const double focal_length1,
    const double focal_length2, const double radial_distortion1,
    const double radial_distortion2, const double outlier_ratio,
    std::vector<Vector2d>* image_1_points_normalized,
    std::vector<Vector2d>* image_2_points_normalized,
    std::vector<Vector2d>* image_1_points,
    std::vector<Vector2d>* image_2_points) {

  for (int i = 0; i < points_3d.size(); i++) {
    Vector2d distorted_point1, distorted_point2;
    Vector3d normalized_points2d_1, normalized_points2d_2;
    bool outlier = false;
    if (rng.Rand(0., 1.) < outlier_ratio) {
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

}  // theia
