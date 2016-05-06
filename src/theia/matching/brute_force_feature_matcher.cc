// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include "theia/matching/brute_force_feature_matcher.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <algorithm>
#include <vector>

#include "theia/matching/distance.h"
#include "theia/matching/feature_matcher_utils.h"
#include "theia/matching/indexed_feature_match.h"

namespace theia {


bool BruteForceFeatureMatcher::MatchImagePair(
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    std::vector<FeatureCorrespondence>* matched_features) {

  const std::vector<Eigen::VectorXf>& descriptors1 = features1.descriptors;
  const std::vector<Keypoint>& keypoints1 = features1.keypoints;
  const std::vector<Eigen::VectorXf>& descriptors2 = features2.descriptors;
  const std::vector<Keypoint>& keypoints2 = features2.keypoints;

  const double sq_lowes_ratio =
      this->options_.lowes_ratio * this->options_.lowes_ratio;

  L2 distance;
  std::vector<IndexedFeatureMatch> matches;

  // Compute forward matches.
  std::vector<IndexedFeatureMatch> temp_matches(descriptors2.size());
  for (int i = 0; i < descriptors1.size(); i++) {
    for (int j = 0; j < descriptors2.size(); j++) {
      temp_matches[j] =
          IndexedFeatureMatch(i, j, distance(descriptors1[i], descriptors2[j]));
    }

    // Get the lowest distance matches.
    std::partial_sort(temp_matches.begin(),
                      temp_matches.begin() + 2,
                      temp_matches.end(),
                      CompareFeaturesByDistance);

    // Add to the matches vector if lowes ratio test is turned off or it is
    // turned on and passes the test.
    if (!this->options_.use_lowes_ratio ||
        temp_matches[0].distance < sq_lowes_ratio * temp_matches[1].distance) {
      matches.emplace_back(temp_matches[0]);
    }
  }

  if (matches.size() < this->options_.min_num_feature_matches) {
    return false;
  }

  // Compute the symmetric matches, if applicable.
  if (this->options_.keep_only_symmetric_matches) {
    std::vector<IndexedFeatureMatch> reverse_matches;
    temp_matches.resize(descriptors1.size());
    // Only compute the distances for the valid matches.
    for (int i = 0; i < descriptors2.size(); i++) {
      for (int j = 0; j < descriptors1.size(); j++) {
        temp_matches[j] = IndexedFeatureMatch(
            i, j, distance(descriptors2[i], descriptors1[j]));
      }

      // Get the lowest distance matches.
      std::partial_sort(temp_matches.begin(),
                        temp_matches.begin() + 2,
                        temp_matches.end(),
                        CompareFeaturesByDistance);

      // Add to the matches vector if lowes ratio test is turned off or it is
      // turned on and passes the test.
      if (!this->options_.use_lowes_ratio ||
          temp_matches[0].distance <
              sq_lowes_ratio * temp_matches[1].distance) {
        reverse_matches.emplace_back(temp_matches[0]);
      }
    }
    IntersectMatches(reverse_matches, &matches);
  }

  if (matches.size() < this->options_.min_num_feature_matches) {
    return false;
  }

  // Convert to FeatureCorrespondences and return true
  matched_features->resize(matches.size());
  for (int i = 0; i < matches.size(); i++) {
    const Keypoint& keypoint1 = keypoints1[matches[i].feature1_ind];
    const Keypoint& keypoint2 = keypoints2[matches[i].feature2_ind];
    matched_features->at(i).feature1 = Feature(keypoint1.x(), keypoint1.y());
    matched_features->at(i).feature2 = Feature(keypoint2.x(), keypoint2.y());
  }
  return true;
}

}  // namespace theia
