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

#include "theia/matching/cascade_hashing_feature_matcher.h"

#include <Eigen/Core>
#include <glog/logging.h>

#include <vector>

#include "theia/matching/cascade_hasher.h"
#include "theia/matching/feature_matcher.h"
#include "theia/matching/feature_matcher_utils.h"
#include "theia/matching/indexed_feature_match.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {

void CascadeHashingFeatureMatcher::AddImage(
    const std::string& image,
    const std::vector<Keypoint>& keypoints,
    const std::vector<Eigen::VectorXf>& descriptors) {
  // This will save the descriptors and keypoints to disk and set up our LRU
  // cache.
  this->template FeatureMatcher<L2>::AddImage(image, keypoints, descriptors);

  if (cascade_hasher_.get() == nullptr) {
    cascade_hasher_.reset(new CascadeHasher());
    CHECK(cascade_hasher_->Initialize(descriptors[0].size()))
        << "Could not initialize the cascade hasher.";
  }

  // Create the hashing information. NOTE: The HashedImage keeps a pointer to
  // the descriptors. This will become invalidated immediately after this
  // function exits so care must be taken going forward.
  hashed_images_[image] =
      cascade_hasher_->CreateHashedSiftDescriptors(descriptors);
}

void CascadeHashingFeatureMatcher::AddImage(
    const std::string& image,
    const std::vector<Keypoint>& keypoints,
    const std::vector<Eigen::VectorXf>& descriptors,
    const CameraIntrinsicsPrior& intrinsics) {
  // This will save the descriptors and keypoints to disk and set up our LRU
  // cache.
  this->template FeatureMatcher<L2>::AddImage(image, keypoints, descriptors, intrinsics);

  if (cascade_hasher_.get() == nullptr) {
    cascade_hasher_.reset(new CascadeHasher());
    CHECK(cascade_hasher_->Initialize(descriptors[0].size()))
        << "Could not initialize the cascade hasher.";
  }

  // Create the hashing information. NOTE: The HashedImage keeps a pointer to
  // the descriptors. This will become invalidated immediately after this
  // function exits so care must be taken going forward.
  hashed_images_[image] =
      cascade_hasher_->CreateHashedSiftDescriptors(descriptors);
}

bool CascadeHashingFeatureMatcher::MatchImagePair(
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    std::vector<FeatureCorrespondence>* matched_features) {
  // Get references to the hashed images for each set of features. NOTE: The
  // HashedImages keep pointers to the original SIFT descriptors but because we
  // use an LRUCache to load features, this pointer is very likely
  // invalidated. We explicitly set it here to ensure the pointer is valid.
  HashedImage& hashed_features1 =
      FindOrDie(hashed_images_, features1.image_name);
  hashed_features1.descriptors = &features1.descriptors;
  HashedImage& hashed_features2 =
      FindOrDie(hashed_images_, features2.image_name);
  hashed_features2.descriptors = &features2.descriptors;

  std::vector<IndexedFeatureMatch> matches;
  cascade_hasher_->MatchImages(hashed_features1,
                               hashed_features2,
                               this->matcher_options_.lowes_ratio,
                               &matches);
  // Only do symmetric matching if enough matches exist to begin with.
  if (matches.size() >= this->matcher_options_.min_num_feature_matches &&
      this->matcher_options_.keep_only_symmetric_matches) {
    std::vector<IndexedFeatureMatch> backwards_matches;
    cascade_hasher_->MatchImages(hashed_features2,
                                 hashed_features1,
                                 this->matcher_options_.lowes_ratio,
                                 &backwards_matches);
    IntersectMatches(backwards_matches, &matches);
  }

  if (matches.size() < this->matcher_options_.min_num_feature_matches) {
    return false;
  }

  // Convert to FeatureCorrespondences and return true;
  const std::vector<Keypoint>& keypoints1 = features1.keypoints;
  const std::vector<Keypoint>& keypoints2 = features2.keypoints;
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
