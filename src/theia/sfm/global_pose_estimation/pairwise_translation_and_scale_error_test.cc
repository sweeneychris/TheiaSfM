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
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include <ceres/rotation.h>
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/sfm/global_pose_estimation/pairwise_translation_and_scale_error.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace {

static const double kScale = 1.0;

void PairwiseTranslationAndScaleErrorTest(const Vector3d& orientation,
                                          const double scale,
                                          const Vector3d& position_1,
                                          const Vector3d& position_2) {
  static const double kTolerance = 1e-8;

  Matrix3d rotation;
  ceres::AngleAxisToRotationMatrix(orientation.data(), rotation.data());

  const Vector3d local_position1 = (rotation * position_1) / scale;
  const Vector3d local_position2 = (rotation * position_2) / scale;
  const Vector3d expected_error =
      (position_2 - position_1) -
    scale * rotation.transpose() * (local_position2 - local_position1);

  // Initialize error function and compute rotation error.
  const PairwiseTranslationAndScaleError translation_error(
      orientation, local_position1, local_position2);
  Vector3d error = Vector3d::Random();
  translation_error(position_1.data(),
                    position_2.data(),
                    &scale,
                    error.data());

  EXPECT_NEAR(error(0), expected_error(0), kTolerance);
  EXPECT_NEAR(error(1), expected_error(1), kTolerance);
  EXPECT_NEAR(error(2), expected_error(2), kTolerance);
}

}  // namespace

TEST(PairwiseTranslationAndScaleError, TranslationNoNoise) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  const Vector3d orientation(0.1, -0.3, 0.2);

  const double scale = (position_2 - position_1).norm();
  PairwiseTranslationAndScaleErrorTest(orientation,
                                       scale,
                                       position_1,
                                       position_2);
}

}  // namespace theia
