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

#include <Eigen/Core>
#include <vector>

#include "theia/matching/brute_force_feature_matcher.h"
#include "theia/matching/distance.h"
#include "theia/matching/feature_matcher.h"
#include "theia/matching/image_pair_match.h"
#include "theia/matching/in_memory_features_and_matches_database.h"
#include "theia/matching/keypoints_and_descriptors.h"

#include "gtest/gtest.h"

namespace theia {

using Eigen::VectorXf;

static const int kNumDescriptors = 10;
static const int kNumDescriptorDimensions = 10;

TEST(BruteForceFeatureMatcherTest, NoOptions) {
  // Set up descriptors.
  KeypointsAndDescriptors features1, features2;
  features1.descriptors.resize(kNumDescriptors);
  features2.descriptors.resize(kNumDescriptors);
  for (int i = 0; i < kNumDescriptors; i++) {
    // Avoid a zero vector.
    features1.descriptors[i] = VectorXf::Constant(kNumDescriptorDimensions, 1);
    features2.descriptors[i] = VectorXf::Constant(kNumDescriptorDimensions, 1);
    features1.descriptors[i].normalize();
    features2.descriptors[i].normalize();
  }

  // Set options.
  FeatureMatcherOptions options;
  options.min_num_feature_matches = 0;
  options.keep_only_symmetric_matches = false;
  options.use_lowes_ratio = false;
  options.perform_geometric_verification = false;

  // Add features.
  features1.keypoints.resize(features1.descriptors.size());
  features2.keypoints.resize(features2.descriptors.size());
  InMemoryFeaturesAndMatchesDatabase database;
  database.PutFeatures("1", features1);
  database.PutFeatures("2", features2);

  BruteForceFeatureMatcher matcher(options, &database);
  matcher.AddImage("1");
  matcher.AddImage("2");

  // Match features
  matcher.MatchImages();

  // Check that the results are valid.
  EXPECT_GT(database.NumMatches(), 0);
}

TEST(BruteForceFeatureMatcherTest, RatioTest) {
  // Set up descriptors.
  KeypointsAndDescriptors features1, features2;
  features1.descriptors.resize(1);
  features2.descriptors.resize(2);

  features1.descriptors[0] =
      VectorXf::Constant(kNumDescriptorDimensions, 1).normalized();

  // Set the two descriptors to be very close to each other so that they do not
  // pass the ratio test.
  features2.descriptors[0] = VectorXf::Constant(kNumDescriptorDimensions, 1);
  features2.descriptors[0](0) = 0.9;
  features2.descriptors[0].normalize();
  features2.descriptors[1] = VectorXf::Constant(kNumDescriptorDimensions, 1);
  features2.descriptors[1](0) = 0.89;
  features2.descriptors[1].normalize();

  // Set options.
  FeatureMatcherOptions options;
  options.min_num_feature_matches = 0;
  options.keep_only_symmetric_matches = false;
  options.use_lowes_ratio = true;
  options.perform_geometric_verification = false;

  // Add features.
  features1.keypoints.resize(features1.descriptors.size());
  features2.keypoints.resize(features2.descriptors.size());

  InMemoryFeaturesAndMatchesDatabase database;
  database.PutFeatures("1", features1);
  database.PutFeatures("2", features2);

  BruteForceFeatureMatcher matcher(options, &database);
  matcher.AddImage("1");
  matcher.AddImage("2");

  // Match features.
  matcher.MatchImages();

  // Check that the results are valid.
  EXPECT_GT(database.NumMatches(), 0);
}

TEST(BruteForceFeatureMatcherTest, SymmetricMatches) {
  // Set up descriptors.
  KeypointsAndDescriptors features1, features2;
  features1.descriptors.resize(2);
  features2.descriptors.resize(2);

  features1.descriptors[0] =
      VectorXf::Constant(kNumDescriptorDimensions, 1).normalized();
  features1.descriptors[1] = VectorXf::Constant(kNumDescriptorDimensions, 0);
  features1.descriptors[1](0) = 1.0;

  // Set the two descriptors to be closer to features1.descriptors[0] so that
  // the symmetric matching produces only 1 match.
  features2.descriptors[0] = VectorXf::Constant(kNumDescriptorDimensions, 1);
  features2.descriptors[0](0) = 0;
  features2.descriptors[0].normalize();
  features2.descriptors[1] = VectorXf::Constant(kNumDescriptorDimensions, 1);
  features2.descriptors[1](1) = 0;
  features2.descriptors[1](2) = 0;
  features2.descriptors[1].normalize();

  // Set options.
  FeatureMatcherOptions options;
  options.min_num_feature_matches = 0;
  options.keep_only_symmetric_matches = true;
  options.use_lowes_ratio = false;
  options.perform_geometric_verification = false;

  // Add features.
  features1.keypoints.resize(features1.descriptors.size());
  features2.keypoints.resize(features2.descriptors.size());

  InMemoryFeaturesAndMatchesDatabase database;
  database.PutFeatures("1", features1);
  database.PutFeatures("2", features2);

  BruteForceFeatureMatcher matcher(options, &database);
  matcher.AddImage("1");
  matcher.AddImage("2");

  // Match features.
  matcher.MatchImages();

  // Check that the results are valid.
  EXPECT_EQ(database.NumMatches(), 1);
}

}  // namespace theia
