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
#include "theia/sfm/pose/pairwise_translation_error.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace {

static const double kRelativeTranslationWeight = 1.0;

void PairwiseTranslationErrorTest(const Vector3d& known_translation,
                                  const double weight,
                                  const Vector3d& position_1,
                                  const Vector3d& position_2) {
  // Compute ground truth angular error.
  Vector3d translation = position_2 - position_1;
  static const double kEpsilon = 1e-8;
  if (translation.norm() > kEpsilon) {
    translation.normalize();
  }

  const Vector3d expected_error = weight * (translation - known_translation);

  // Initialize error function and compute rotation error.
  const PairwiseTranslationError translation_error(known_translation, weight);
  Vector3d error = Vector3d::Zero();
  translation_error(position_1.data(),
                    position_2.data(),
                    error.data());

  EXPECT_DOUBLE_EQ(error(0), expected_error(0));
  EXPECT_DOUBLE_EQ(error(1), expected_error(1));
  EXPECT_DOUBLE_EQ(error(2), expected_error(2));
}

}  // namespace

TEST(PairwiseTranslationError, NoTranslation) {
  const Vector3d position_1(1.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  const Vector3d relative_translation(0.0, 0.0, 0.0);

  // Initialize error function and compute rotation error.
  const PairwiseTranslationError translation_error(relative_translation, 1.0);
  Vector3d error = Vector3d::Zero();
  EXPECT_TRUE(translation_error(position_1.data(),
                                position_2.data(),
                                error.data()));
}

TEST(PairwiseTranslationError, TranslationNoNoise) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  const Vector3d relative_translation = (position_2 - position_1).normalized();

  PairwiseTranslationErrorTest(relative_translation,
                               kRelativeTranslationWeight,
                               position_1,
                               position_2);
}

TEST(PairwiseTranslationError, TranslationWithNoise) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  Vector3d relative_translation = (position_2 - position_1).normalized();

  // Add noise.
  relative_translation =
      (relative_translation + Vector3d(0.01, 0.01, 0.01)).normalized();

  PairwiseTranslationErrorTest(relative_translation,
                               kRelativeTranslationWeight,
                               position_1,
                               position_2);
}

TEST(PairwiseTranslationError, NontrivialWeight) {
  const Vector3d position_1(0.0, 0.0, 0.0);
  const Vector3d position_2(1.0, 0.0, 0.0);
  Vector3d relative_translation = (position_2 - position_1).normalized();

  // Add noise.
  relative_translation =
      (relative_translation + Vector3d(0.01, 0.01, 0.01)).normalized();

  static const double kNontrivialWeight = 1.1;
  PairwiseTranslationErrorTest(relative_translation,
                               kNontrivialWeight,
                               position_1,
                               position_2);
}

}  // namespace theia
