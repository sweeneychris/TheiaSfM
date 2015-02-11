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

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/sfm/transformation/align_point_clouds.h"

namespace theia {
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::RowMajor;
using Eigen::Vector3d;

namespace {
double kEpsilon = 1e-6;

void UmeyamaSimpleTest() {
  std::vector<Vector3d> left = {
      Vector3d(0.4, -3.105, 2.147),
      Vector3d(1.293, 7.1982, -.068),
      Vector3d(-5.34, 0.708, -3.69),
      Vector3d(-.345, 1.987, 0.936),
      Vector3d(0.93, 1.45, 1.079),
      Vector3d(-3.15, -4.73, 2.49),
      Vector3d(2.401, -2.03, -1.87),
      Vector3d(3.192, -.573, 0.1),
      Vector3d(-2.53, 3.07, -5.19)};

  const Matrix3d rotation_mat =
      Eigen::AngleAxisd(DegToRad(15.0), Vector3d(1.0, -2.7, 1.9).normalized())
          .toRotationMatrix();
  const Vector3d translation_vec(0, 2, 2);
  const double expected_scale = 1.5;

  // Transform the points.
  std::vector<Vector3d> right;
  for (int i = 0; i < left.size(); i++) {
    Vector3d transformed_point =
        expected_scale * rotation_mat * left[i] + translation_vec;
    right.emplace_back(transformed_point);
  }

  // Compute the similarity transformation.
  Matrix3d rotation;
  Vector3d translation;
  double scale;
  AlignPointCloudsUmeyama(left, right, &rotation, &translation, &scale);

  // Ensure the calculated transformation is the same as the one we set.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_LT(std::abs(rotation(i, j) - rotation_mat(i, j)), kEpsilon);
    }
    ASSERT_LT(std::abs(translation(i) - translation_vec(i)), kEpsilon);
  }
  ASSERT_LT(fabs(expected_scale - scale), kEpsilon);
}

}  // namespace

TEST(AlignPointCloudsUmeyama, SimpleTest) {
  UmeyamaSimpleTest();
}

}  // namespace theia
