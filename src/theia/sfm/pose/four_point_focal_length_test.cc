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
// Please contact the author of this file if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include <math.h>
#include <glog/logging.h>
#include <Eigen/Core>
#include <ctime>
#include <random>
#include "gtest/gtest.h"

#include "theia/sfm/pose/four_point_focal_length.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"

namespace {
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector3d;

theia::RandomNumberGenerator rng(151);

void P4pfTestWithNoise(const Matrix3d& gt_rotation,
                       const Vector3d& gt_translation,
                       const double focal_length,
                       const std::vector<Vector3d>& world_points_vector,
                       const double noise,
                       const double reproj_tolerance) {
  Map<const Matrix<double, 3, 4> > world_points(world_points_vector[0].data());

  // Camera intrinsics matrix.
  const Matrix3d camera_matrix =
      Eigen::DiagonalMatrix<double, 3>(focal_length, focal_length, 1.0);
  // Create the projection matrix P = K * [R t].
  Matrix<double, 3, 4> gt_projection;
  gt_projection << gt_rotation, gt_translation;
  gt_projection = camera_matrix * gt_projection;

  // Reproject 3D points to get undistorted image points.
  std::vector<Eigen::Vector2d> image_points_vector(5);
  Map<Matrix<double, 2, 4> > image_point(image_points_vector[0].data());
  image_point = (gt_projection * world_points.colwise().homogeneous()).colwise()
      .hnormalized();

  // Add noise to distorted image points.
  if (noise) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, noise);
    for (int i = 0; i < 4; i++) {
      image_point.col(i).x() += distribution(generator);
      image_point.col(i).y() += distribution(generator);
    }
  }

  // Run P5pf algorithm.
  std::vector<Matrix<double, 3, 4> > soln_projection;
  int num_solns = theia::FourPointPoseAndFocalLength(
      image_points_vector, world_points_vector, &soln_projection);
  ASSERT_GT(num_solns, 0);

  bool matched_transform = false;
  for (int i = 0; i < num_solns; ++i) {
    matched_transform = true;
    // Check that the reprojection error is very small.
    for (int n = 0; n < 4; n++) {
      Vector3d reproj_point =
          soln_projection[i] * world_points.col(n).homogeneous();
      const double reproj_error =
          (reproj_point.hnormalized() - image_point.col(n)).norm();
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
  const double focal_length = 800;

  const double x = -0.10;  // rotation of the view around x axis
  const double y = -0.20;  // rotation of the view around y axis
  const double z = 0.30;   // rotation of the view around z axis

  // Create a ground truth pose.
  Matrix3d Rz, Ry, Rx;
  Rz << cos(z), sin(z), 0,
        -sin(z), cos(z), 0,
        0, 0, 1;
  Ry << cos(y), 0, -sin(y),
        0, 1, 0,
        sin(y), 0, cos(y);
  Rx << 1, 0, 0,
        0, cos(x), sin(x),
        0, -sin(x), cos(x);
  const Matrix3d gt_rotation = Rz * Ry * Rx;
  const Vector3d gt_translation = Vector3d(-0.00950692, 0.0171496, 0.0508743);

  // Create 3D world points that are viable based on the camera intrinsics and
  // extrinsics.
  std::vector<Vector3d> world_points_vector = { Vector3d(-1.0, 0.5, 1.2),
                                                Vector3d(-0.79, -0.68, 1.9),
                                                Vector3d(1.42, 1.01, 2.19),
                                                Vector3d(0.87, -0.49, 0.89) };

  P4pfTestWithNoise(gt_rotation, gt_translation, focal_length,
                    world_points_vector, noise, reproj_tolerance);
}

void RandomTestWithNoise(const double noise, const double reproj_tolerance) {
  const double kBaseline = 0.25;

  // focal length (values used in the ICCV paper)
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double focal_length = distribution(generator) * 50.0 + 600;

  // Rotation areound x, y, z axis.
  const double x = distribution(generator) * 0.5 - 0.25;
  const double y = distribution(generator) * 0.5 - 0.25;
  const double z = distribution(generator) * 0.5 - 0.25;

  // Create a ground truth pose.
  Matrix3d Rz, Ry, Rx;
  Rz << cos(z), sin(z), 0,
        -sin(z), cos(z), 0,
        0, 0, 1;
  Ry << cos(y), 0, -sin(y),
        0, 1, 0,
        sin(y), 0, cos(y);
  Rx << 1, 0, 0,
        0, cos(x), sin(x),
        0, -sin(x), cos(x);
  const Matrix3d gt_rotation = Rz * Ry * Rx;
  const Vector3d gt_translation = rng.RandVector3d() * kBaseline;

  // Create 3D world points that are viable based on the camera intrinsics and
  // extrinsics.
  std::vector<Vector3d> world_points_vector(4);
  Map<Matrix<double, 3, 4> > world_points(world_points_vector[0].data());
  world_points.row(2) = 2.0 * rng.RandVector4d().transpose().array() + 2.0;
  world_points.row(1) = 2.0 * rng.RandVector4d().transpose();
  world_points.row(0) = 2.0 * rng.RandVector4d().transpose();

  P4pfTestWithNoise(gt_rotation, gt_translation, focal_length,
                    world_points_vector, noise, reproj_tolerance);
}

TEST(P4pf, BasicTest) {
  BasicTest(0.0, 1e-4);
}

TEST(P4pf, BasicNoiseTest) {
  BasicTest(0.5, 10.0);
}

TEST(P4pf, RandomTest) {
  RandomTestWithNoise(0.0, 0.1);
}

}  // namespace
