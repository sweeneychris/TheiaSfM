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

void DistortPoint(const Eigen::Vector3d& point_in_camera,
                  const double focal_length, const double radial_distortion,
                  Eigen::Vector2d& distorted_point) {
  Eigen::Vector2d point_in_cam_persp_div;
  point_in_cam_persp_div[0] =
      focal_length * point_in_camera[0] / point_in_camera[2];
  point_in_cam_persp_div[1] =
      focal_length * point_in_camera[1] / point_in_camera[2];
  // see also division_undistortion_camera_model.h
  const double r_u_sq = point_in_cam_persp_div[0] * point_in_cam_persp_div[0] +
                        point_in_cam_persp_div[1] * point_in_cam_persp_div[1];

  const double denom = 2.0 * radial_distortion * r_u_sq;
  const double inner_sqrt = 1.0 - 4.0 * radial_distortion * r_u_sq;

  // If the denominator is nearly zero then we can evaluate the distorted
  // coordinates as k or r_u^2 goes to zero. Both evaluate to the identity.
  if (std::abs(denom) < std::numeric_limits<double>::epsilon() ||
      inner_sqrt < 0.0) {
    distorted_point = point_in_cam_persp_div;
  } else {
    const double scale = (1.0 - std::sqrt(inner_sqrt)) / denom;
    distorted_point = point_in_cam_persp_div * scale;
  }
}

void UndistortPoint(const Eigen::Vector2d& distorted_point,
                    const double focal_length, const double radial_distortion,
                    Eigen::Vector3d& undistorted_point) {
  double px = distorted_point[0];
  double py = distorted_point[1];
  // The squared radius of the distorted image point.
  const double r_d_sq = px * px + py * py;

  const double undistortion = 1.0 / (1.0 + radial_distortion * r_d_sq);
  undistorted_point[0] = px * undistortion / focal_length;
  undistorted_point[1] = py * undistortion / focal_length;
  undistorted_point[2] = 1.0;
}

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

Vector3d ProjectCameraToCamera(const Matrix3d& H, const Vector3d& X, Vector3d* Y) {
  (*Y) = H * X;
  (*Y) /= (*Y)(2);
}

double CheckRadialSymmetricError(
    const RadialHomographyResult& radial_homography, const Vector2d& pt_left,
    const Vector2d& pt_right, const double focal_length1,
    const double focal_length2) {
  Vector3d bearing_vector_left, bearing_vector_right;
  UndistortPoint(pt_left, focal_length1, radial_homography.l1,
                 bearing_vector_left);
  UndistortPoint(pt_right, focal_length2, radial_homography.l2,
                 bearing_vector_right);

  Eigen::Vector3d ray2_in_1, ray1_in_2;
  ProjectCameraToCamera(radial_homography.H, bearing_vector_right, &ray2_in_1);
  ProjectCameraToCamera(radial_homography.H.inverse(), bearing_vector_left, &ray1_in_2);

  Vector2d pt_left_from_right, pt_right_from_left;
  DistortPoint(ray2_in_1, focal_length1, radial_homography.l1,
               pt_left_from_right);
  DistortPoint(ray1_in_2, focal_length2, radial_homography.l2,
               pt_right_from_left);

  double dleft_x = pt_left(0) - pt_left_from_right(0);
  double dleft_y = pt_left(1) - pt_left_from_right(1);

  double dright_x = pt_right(0) - pt_right_from_left(0);
  double dright_y = pt_right(1) - pt_right_from_left(1);

  double sym_error = 0.5 * (dleft_x * dleft_x + dleft_y * dleft_y +
                            dright_x * dright_x + dright_y * dright_y);

  return sym_error;
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
    radial_homography_result[i].l1 /= (focal_length1 * focal_length1);
    radial_homography_result[i].l2 /= (focal_length2 * focal_length2);

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
