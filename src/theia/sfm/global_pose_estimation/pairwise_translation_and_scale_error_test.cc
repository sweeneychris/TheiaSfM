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

#include <ceres/ceres.h>
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/sfm/global_pose_estimation/pairwise_translation_and_scale_error.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;

namespace {

static const double kScale = 1.0;

void PairwiseTranslationAndScaleErrorTest(const Vector3d& known_translation,
                                          const double scale,
                                          const Vector3d& position_1,
                                          const Vector3d& position_2) {
  // Compute ground truth angular error.
  Vector3d translation = position_2 - position_1;
  const Vector3d expected_error = translation - scale * known_translation;

  // Initialize error function and compute rotation error.
  const PairwiseTranslationAndScaleError translation_error(known_translation);
  Vector4d error = Vector4d::Zero();
  double l1_weight = 1.0;
  translation_error(position_1.data(),
                    position_2.data(),
                    &scale,
                    &l1_weight,
                    error.data());

  EXPECT_DOUBLE_EQ(error(0), expected_error(0));
  EXPECT_DOUBLE_EQ(error(1), expected_error(1));
  EXPECT_DOUBLE_EQ(error(2), expected_error(2));
}

}  // namespace

TEST(PairwiseTranslationAndScaleError, TranslationNoNoise) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  const Vector3d relative_translation = (position_2 - position_1).normalized();

  const double scale = (position_2 - position_1).norm();
  PairwiseTranslationAndScaleErrorTest(relative_translation,
                                       scale,
                                       position_1,
                                       position_2);
}

TEST(PairwiseTranslationAndScaleError, TranslationWithNoise) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  Vector3d relative_translation = (position_2 - position_1).normalized();

  // Add noise.
  relative_translation =
      (relative_translation + Vector3d(0.01, 0.01, 0.01)).normalized();
  const double scale = (position_2 - position_1).norm();
  PairwiseTranslationAndScaleErrorTest(relative_translation,
                                       scale,
                                       position_1,
                                       position_2);
}

TEST(PairwiseTranslationAndScaleError, TranslationNoNoiseAndBadScale) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(5.0, 0.0, 0.0);
  const Vector3d relative_translation = (position_2 - position_1).normalized();

  PairwiseTranslationAndScaleErrorTest(relative_translation,
                                       kScale,
                                       position_1,
                                       position_2);
}

TEST(PairwiseTranslationAndScaleError, TranslationWithNoiseAndBadScale) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(5.0, 0.0, 0.0);
  Vector3d relative_translation = (position_2 - position_1).normalized();

  // Add noise.
  relative_translation =
      (relative_translation + Vector3d(0.01, 0.01, 0.01)).normalized();

  PairwiseTranslationAndScaleErrorTest(relative_translation,
                               kScale,
                               position_1,
                               position_2);
}

}  // namespace theia
