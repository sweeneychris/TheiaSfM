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
#include <algorithm>
#include <vector>
#include "gtest/gtest.h"

#include "theia/math/util.h"
#include "theia/util/random.h"
#include "theia/sfm/estimators/essential_matrix_estimator.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/test_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/test/test_utils.h"

namespace theia {
namespace {
using Eigen::AngleAxisd;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector2d;
using Eigen::Vector3d;

// Tests that the five point algorithm works correctly be checking the epipolar
// errors of the returned solutions and ensuring that at least one of the
// candidate solutions corresponds to the essential matrix we used to construct
// the problem.
void TestEssentialMatrixEstimator(const Vector3d points_3d[5],
                                  const Matrix3d& expected_rotation,
                                  const Vector3d& expected_translation,
                                  const double ematrix_tolerance) {
  InitRandomGenerator();

  // Calculates the image points in both views.
  std::vector<FeatureCorrespondence> features(5);
  for (int i = 0; i < 5; ++i) {
    const Vector3d proj_3d =
        expected_rotation * points_3d[i] + expected_translation;
    features[i].feature1 = points_3d[i].hnormalized();
    features[i].feature2 = proj_3d.hnormalized();
  }

  // Calculates the essential matrix, this may return multiple solutions.
  const Matrix3d gt_ematrix =
      CrossProductMatrix(expected_translation) * expected_rotation;

  EssentialMatrixEstimator ematrix_estimator;
  std::vector<Matrix3d> soln_ematrices;
  EXPECT_TRUE(ematrix_estimator.EstimateModel(features, &soln_ematrices));
  CHECK_GT(soln_ematrices.size(), 0);

  // Among the returned solutions verify that at least one is close to the
  // expected translation and rotation.
  bool matched_transform = false;
  const double kEpipolarTolerance = 1e-8;
  for (int n = 0; n < soln_ematrices.size(); ++n) {
    // All solutions should have valid epipolar constraints.
    for (int i = 0; i < 5; i++) {
      const double gt_error = SquaredSampsonDistance(
          soln_ematrices[n],
          features[i].feature1,
          features[i].feature2);
      const double estimated_error = ematrix_estimator.Error(features[i],
                                                             soln_ematrices[n]);
      EXPECT_DOUBLE_EQ(estimated_error, gt_error);
      EXPECT_LT(estimated_error, kEpipolarTolerance);
    }

    if (test::ArraysEqualUpToScale(9, soln_ematrices[n].data(),
                                   gt_ematrix.data(), ematrix_tolerance)) {
      matched_transform = true;
    }
  }
  EXPECT_TRUE(matched_transform);
}


TEST(EssentialMatrixEstimator, BasicTest) {
  // Ground truth essential matrix.
  const Vector3d points_3d[5] = { Vector3d(-1.0, 3.0, 3.0),
                                  Vector3d(1.0, -1.0, 2.0),
                                  Vector3d(3.0, 1.0, 2.5),
                                  Vector3d(-1.0, 1.0, 2.0),
                                  Vector3d(2.0, 1.0, 3.0) };
  const Matrix3d soln_rotation = Quaterniond(
      AngleAxisd(DegToRad(13.0), Vector3d(0.0, 0.0, 1.0))).toRotationMatrix();
  const Vector3d soln_translation(1.0, 1.0, 1.0);
  const double kEMatrixTolerance = 1e-4;
  TestEssentialMatrixEstimator(points_3d,
                               soln_rotation,
                               soln_translation,
                               kEMatrixTolerance);
}

}  // namespace
}  // namespace theia
