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
// January 2019

#include <glog/logging.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/sfm/pose/six_point_radial_distortion_homography.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/util/random.h"

namespace theia {
namespace {
using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector2d;
using Eigen::Vector3d;

RandomNumberGenerator rng(53);

// Creates a test scenario from ground truth 3D points and ground truth rotation
// and translation. Projection noise is optional (set to 0 for no
// noise).
void GenerateDistortedImagePoints(
    const std::vector<Vector3d>& points_3d,
    const double projection_noise_std_dev, const Quaterniond& expected_rotation,
    const Vector3d& expected_translation, const double focal_length1,
    const double focal_length2, const double radial_distortion1,
    const double radial_distortion2,
    std::vector<Vector2d>* image_1_points_normalized,
    std::vector<Vector2d>* image_2_points_normalized,
    std::vector<Vector2d>* image_1_points,
    std::vector<Vector2d>* image_2_points) {
  image_1_points->reserve(points_3d.size());
  image_2_points->reserve(points_3d.size());

  for (int i = 0; i < points_3d.size(); i++) {
    Vector3d point3_cam1 = points_3d[i];
    Vector3d point3_cam2 =
        (expected_rotation * points_3d[i] + expected_translation);

    Vector2d distorted_point1, distorted_point2;
    DistortPoint(point3_cam1, focal_length1, radial_distortion1,
                 distorted_point1);
    DistortPoint(point3_cam2, focal_length2, radial_distortion2,
                 distorted_point2);

    image_1_points->push_back(distorted_point1);
    image_2_points->push_back(distorted_point2);

    if (projection_noise_std_dev) {
      for (int i = 0; i < points_3d.size(); i++) {
        AddNoiseToProjection(projection_noise_std_dev, &rng,
                             &((*image_1_points)[i]));
        AddNoiseToProjection(projection_noise_std_dev, &rng,
                             &((*image_2_points)[i]));
      }
    }

    Vector3d points2d_1, points2d_2;
    // normalize points with focal length (and principal point) estimate for
    // estimation
    UndistortPoint(distorted_point1, focal_length1, 0.0, points2d_1);
    UndistortPoint(distorted_point2, focal_length2, 0.0, points2d_2);

    image_1_points_normalized->push_back(points2d_1.hnormalized());
    image_2_points_normalized->push_back(points2d_2.hnormalized());
  }
}

// Run a test for the homography with at least 4 points.
void SixPointHomographyWithNoiseTest(
    const std::vector<Vector3d>& points_3d,
    const double projection_noise_std_dev, const Quaterniond& expected_rotation,
    const Vector3d& expected_translation, const double kMaxSymmetricError,
    const double focal_length1, const double focal_length2,
    const double radial_distortion1, const double radial_distortion2) {
  std::vector<Vector2d> image_1_points, image_points_1_normalized;
  std::vector<Vector2d> image_2_points, image_points_2_normalized;
  // generate distorted points for both cameras
  GenerateDistortedImagePoints(
      points_3d, projection_noise_std_dev, expected_rotation,
      expected_translation, focal_length1, focal_length2, radial_distortion1,
      radial_distortion2, &image_points_1_normalized,
      &image_points_2_normalized, &image_1_points, &image_2_points);
  // Compute two-sided radial distortion homography matrix.
  std::vector<RadialHomographyResult> radial_homography_result;
  EXPECT_TRUE(SixPointRadialDistortionHomography(image_points_1_normalized,
                                                 image_points_2_normalized,
                                                 &radial_homography_result));
  bool one_correct_solution = false;
  for (int i = 0; i < radial_homography_result.size(); ++i) {
    // we need to scale the radial distortion values,
    // since we used normalized image points for estimation
    //radial_homography_result[i].l1 /= (focal_length1 * focal_length1);
    //radial_homography_result[i].l2 /= (focal_length2 * focal_length2);

    double sym_error = 0.0;
    for (int p = 0; p < image_1_points.size(); ++p) {
      sym_error += CheckRadialSymmetricError(
          radial_homography_result[i], image_1_points[p], image_2_points[p],
          focal_length1, focal_length2);
    }
    sym_error /= (double)image_1_points.size();
    if (sym_error < kMaxSymmetricError) {
      one_correct_solution = true;
    }
  }
  EXPECT_TRUE(one_correct_solution);
}

void BasicTest() {
  const std::vector<Vector3d> points_3d = {
      Vector3d(-1.0, 3.0, 1.0), Vector3d(1.0, -1.0, 1.0),
      Vector3d(-1.0, 1.0, 1.0), Vector3d(2.0, 1.0, 1.0),
      Vector3d(3.0, 1.0, 1.0),  Vector3d(2.0, 2.0, 1.0)};

  const Quaterniond soln_rotation(
      AngleAxisd(DegToRad(10.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double kNoise = 0.0 / 512.0;
  const double kMaxSymmetricError = 1e-4;

  const double focal_length1 = 1500.0;
  const double focal_length2 = 1600.0;
  const double radial_distortion1 = -1e-7;
  const double radial_distortion2 = -2e-7;

  SixPointHomographyWithNoiseTest(
      points_3d, kNoise, soln_rotation, soln_translation, kMaxSymmetricError,
      focal_length1, focal_length2, radial_distortion1, radial_distortion2);
}

TEST(SixPointRadialHomography, BasicTest) { BasicTest(); }

TEST(SixPointRadialHomography, NoiseTest) {
  const std::vector<Vector3d> points_3d = {
      Vector3d(-1.0, 3.0, 1.0), Vector3d(1.0, -1.0, 1.0),
      Vector3d(-1.0, 1.0, 1.0), Vector3d(2.0, 1.0, 1.0),
      Vector3d(3.0, 1.0, 1.0),  Vector3d(2.0, 2.0, 1.0)};

  const Quaterniond soln_rotation(
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0)));
  const Vector3d soln_translation(0.0, 1.0, 1.0);
  const double kNoise = 0.5;
  const double kMaxSymmetricError = 2.0;

  const double focal_length1 = 1500.0;
  const double focal_length2 = 1600.0;
  const double radial_distortion1 = -1e-7;
  const double radial_distortion2 = -2e-7;

  SixPointHomographyWithNoiseTest(
      points_3d, kNoise, soln_rotation, soln_translation, kMaxSymmetricError,
      focal_length1, focal_length2, radial_distortion1, radial_distortion2);
}

}  // namespace
}  // namespace theia
