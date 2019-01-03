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
// Please contact the author of this file if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// December 2018

#include <math.h>
#include <Eigen/Core>
#include <random>
#include "gtest/gtest.h"

#include "theia/sfm/pose/four_point_focal_length_radial_distortion.h"
#include "theia/test/test_utils.h"

namespace theia {

namespace {

const int min_nr_points = 4;

using Eigen::Array;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix4d;

void P4pfrTestWithNoise(const Matrix3d &gt_rotation,
                        const Vector3d &gt_translation,
                        const double focal_length,
                        const double radial_distortion,
                        const std::vector<Vector3d> &world_points_vector,
                        const double noise, const double reproj_tolerance) {
  // Camera intrinsics matrix.
  const Matrix3d camera_matrix =
      Eigen::DiagonalMatrix<double, 3>(focal_length, focal_length, 1.0);
  Matrix<double, 3, 4> gt_transformation, gt_projection;
  gt_transformation << gt_rotation, gt_translation;
  gt_projection = camera_matrix * gt_transformation;

  // Reproject 3D points to get undistorted image points.

  Matrix<double, 2, min_nr_points> undistorted_image_point;
  for (int i = 0; i < min_nr_points; ++i) {
    undistorted_image_point.col(i) =
        (gt_projection * world_points_vector[i].homogeneous()).hnormalized();
  }
  // Determine radius of undistorted points and use that to compute the radius
  // of the distorted points.
  Array<double, 1, min_nr_points> radius_undistorted =
      undistorted_image_point.colwise().norm();
  Array<double, 1, min_nr_points> radius_distorted =
      (1.0 -
       (1.0 - 4.0 * radial_distortion * radius_undistorted.square()).sqrt()) /
      (2.0 * radial_distortion * radius_undistorted);
  Array<double, 1, min_nr_points> distortion_vec =
      radius_distorted / radius_undistorted;

  // Apply radial distortion.
  std::vector<Vector2d> distorted_image_points_vector(min_nr_points);
  Map<Matrix<double, 2, min_nr_points>> distorted_image_point(
      distorted_image_points_vector[0].data());
  distorted_image_point = undistorted_image_point.cwiseProduct(
      distortion_vec.matrix().replicate<2, 1>());

  // Add noise to distorted image points.
  if (noise) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, noise);
    for (int i = 0; i < min_nr_points; i++) {
      distorted_image_point.col(i).x() += distribution(generator);
      distorted_image_point.col(i).y() += distribution(generator);
    }
  }

  // Run the P4Pfr algorithm.
  std::vector<Matrix3d> soln_rotations;
  std::vector<Vector3d> soln_translations;
  std::vector<double> soln_focal_lenghts;
  std::vector<double> soln_radial_distortions;

  CHECK(FourPointsPoseFocalLengthRadialDistortion(
      distorted_image_points_vector, world_points_vector, &soln_rotations,
      &soln_translations, &soln_radial_distortions, &soln_focal_lenghts));
  bool matched_transform = false;
  for (int i = 0; i < soln_radial_distortions.size(); ++i) {
    matched_transform = true;
    // Check the reprojection error.
    for (int n = 0; n < min_nr_points; n++) {
      const double distortion_w =
          1.0 +
          soln_radial_distortions[i] *
              distorted_image_point.col(n).squaredNorm();
      Vector2d undist_pt = distorted_image_point.col(n) / distortion_w;
      Vector3d reproj_pt =
          soln_rotations[i] * world_points_vector[n] + soln_translations[i];
      Matrix3d K = Vector3d(soln_focal_lenghts[i], soln_focal_lenghts[i], 1.0)
                       .asDiagonal();
      reproj_pt = K * reproj_pt;
      const double reproj_error =
          (undist_pt - reproj_pt.hnormalized()).squaredNorm();
      if (reproj_error > reproj_tolerance) {
        matched_transform = false;
        break;
      }
    }

    if (matched_transform) {
      break;
    }
  }
  // One of the solutions must have been a valid solution.
  EXPECT_TRUE(matched_transform);
}

void BasicTest(const double noise, const double reproj_tolerance) {
  // focal length
  const double focal_length = 1000;
  // radial distortion
  const double radial_distortion = -1e-6;

  const double x = -0.10;  // rotation of the view around x axis
  const double y = -0.20;  // rotation of the view around y axis
  const double z = 0.30;   // rotation of the view around z axis

  // Create a ground truth pose.
  Matrix3d Rz, Ry, Rx;
  Rz << cos(z), sin(z), 0, -sin(z), cos(z), 0, 0, 0, 1;
  Ry << cos(y), 0, -sin(y), 0, 1, 0, sin(y), 0, cos(y);
  Rx << 1, 0, 0, 0, cos(x), sin(x), 0, -sin(x), cos(x);
  const Matrix3d gt_rotation = Rz * Ry * Rx;
  const Vector3d gt_translation =
      Vector3d(-0.00950692, 000.0171496, 000.0508743);

  // Create 3D world points that are viable based on the camera intrinsics and
  // extrinsics.
  std::vector<Vector3d> world_points_vector(min_nr_points);
  Map<Matrix<double, 3, min_nr_points>> world_points(
      world_points_vector[0].data());
  world_points << -0.42941, 0.000621211, -0.350949, -1.45205, 0.415794,
      -0.556605, -1.92898, -1.89976, 1.4949, 0.838307, 1.41972, 1.25756;
  P4pfrTestWithNoise(gt_rotation, gt_translation, focal_length,
                     radial_distortion, world_points_vector, noise,
                     reproj_tolerance);
}

void PlanarTestWithNoise(const double noise, const double reproj_tolerance) {
  // focal length
  const double focal_length = 1000;
  // radial distortion
  const double radial_distortion = -1e-6;

  const double size = 100;
  const double depth = 150;

  const double x = -0.10;  // rotation of the view around x axis
  const double y = -0.20;  // rotation of the view around y axis
  const double z = 0.30;   // rotation of the view around z axis

  // Create a ground truth pose.
  Matrix3d Rz, Ry, Rx;
  Rz << cos(z), sin(z), 0, -sin(z), cos(z), 0, 0, 0, 1;
  Ry << cos(y), 0, -sin(y), 0, 1, 0, sin(y), 0, cos(y);
  Rx << 1, 0, 0, 0, cos(x), sin(x), 0, -sin(x), cos(x);
  const Matrix3d gt_rotation = Rz * Ry * Rx;
  const Vector3d gt_translation =
      Vector3d(-0.00950692, 000.0171496, 000.0508743);

  // Create 3D world points that are viable based on the camera intrinsics and
  // extrinsics.
  std::vector<Vector3d> world_points_vector(min_nr_points);
  world_points_vector[0] = Vector3d(-size / 2, -size / 2, depth);
  world_points_vector[1] = Vector3d(size / 2, -size / 2, depth);
  world_points_vector[2] = Vector3d(size / 2, size / 2, depth);
  world_points_vector[3] = Vector3d(-size / 2, size / 2, depth);

  P4pfrTestWithNoise(gt_rotation, gt_translation, focal_length,
                     radial_distortion, world_points_vector, noise,
                     reproj_tolerance);
}

TEST(P4Pfr, BasicTest) { BasicTest(0.0, 1e-12); }

TEST(P4Pfr, BasicNoiseTest) { BasicTest(0.5, 5); }

TEST(P4Pfr, PlanarTestNoNoise) { PlanarTestWithNoise(0.0, 1e-12); }

TEST(P4Pfr, PlanarTestWithNoise) { PlanarTestWithNoise(0.5, 5); }
}
}  // namespace theia
