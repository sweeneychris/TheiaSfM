// Copyright (C) 2013 The Regents of the University of California (Regents).
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
#include <math.h>
#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/test/test_utils.h"
#include "theia/util/random.h"
#include "theia/sfm/pose/perspective_three_point.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/types.h"

namespace theia {

using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

void PoseFromThreeCalibratedTest(const double noise) {
  // Projection matrix.
  const Matrix3d gt_rotation =
      (Eigen::AngleAxisd(15.0, Vector3d(1.0, 0.0, 0.0)) *
       Eigen::AngleAxisd(-10.0, Vector3d(0.0, 1.0, 0.0))).toRotationMatrix();
  const Vector3d gt_translation(0.3, -1.7, 1.15);
  Matrix3x4d projection_mat;
  projection_mat << gt_rotation, gt_translation;

  // Points in the 3D scene.
  const Vector3d kPoints3d[3] = { Vector3d(-0.3001, -0.5840, 1.2271),
                                  Vector3d(-1.4487, 0.6965, 0.3889),
                                  Vector3d(-0.7815, 0.7642, 0.1257)};

  // Points in the camera view.
  Vector2d kPoints2d[3];
  for (int i = 0; i < 3; i++) {
    kPoints2d[i] =
        (projection_mat * kPoints3d[i].homogeneous()).eval().hnormalized();
    if (noise) {
      AddNoiseToProjection(noise, &kPoints2d[i]);
    }
  }

  std::vector<Matrix3d> rotations;
  std::vector<Vector3d> translations;
  CHECK(
      PoseFromThreePoints(kPoints2d, kPoints3d, &rotations, &translations));

  bool matched_transform = false;
  for (int i = 0; i < rotations.size(); ++i) {
    // Check that the rotation and translation are close.
    double angular_diff = RadToDeg(Eigen::Quaterniond(
        rotations[i]).angularDistance(Eigen::Quaterniond(gt_rotation)));
    double trans_diff = ((-gt_rotation * gt_translation) -
                         (-rotations[i] * translations[i])).norm();
    bool rot_match = angular_diff < 1.0;
    bool trans_match = trans_diff < 0.1;
    if (rot_match && trans_match) {
      matched_transform = true;

      Matrix3x4d soln_proj;
      soln_proj << rotations[i], translations[i];
      // Check the reprojection error.
      for (int j = 0; j < 3; j++) {
        const Vector3d projected_pt = soln_proj * kPoints3d[j].homogeneous();
        EXPECT_LT((kPoints2d[j] - projected_pt.hnormalized()).norm() * 800.0,
                  2.0);
      }
    }
  }
  EXPECT_TRUE(matched_transform);
}

TEST(P3P, PoseFromThreeCalibrated) {
  PoseFromThreeCalibratedTest(0.0 / 800.0);
}

TEST(P3P, PoseFromThreeCalibratedNoise) {
  PoseFromThreeCalibratedTest(1.0 / 800.0);
}

}  // namespace theia
