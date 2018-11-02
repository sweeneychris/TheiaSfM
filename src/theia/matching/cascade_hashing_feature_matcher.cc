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

#include <memory>
#include <mutex>  // NOLINT
#include <string>
#include <thread>  // NOLINT
#include <vector>

#include "theia/matching/cascade_hasher.h"
#include "theia/matching/feature_matcher.h"
#include "theia/matching/feature_matcher_utils.h"
#include "theia/matching/features_and_matches_database.h"
#include "theia/matching/indexed_feature_match.h"
#include "theia/util/lru_cache.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {
CascadeHashingFeatureMatcher::CascadeHashingFeatureMatcher(
    const FeatureMatcherOptions& options,
    FeaturesAndMatchesDatabase* features_and_matches_database)
    : FeatureMatcher(options, features_and_matches_database) {
  // Initialize the cache.
  const std::function<std::shared_ptr<HashedImage>(const std::string&)>
      fetch_hashed_images =
          std::bind(&CascadeHashingFeatureMatcher::FetchHashedImage,
                    this,
                    std::placeholders::_1);
  // TODO: Stop hardcoding this!
  static constexpr int kNumImagesInCache = 256;
  hashed_images_.reset(
      new HashedImageCache(fetch_hashed_images, kNumImagesInCache));
}

CascadeHashingFeatureMatcher::~CascadeHashingFeatureMatcher() {}

std::shared_ptr<HashedImage> CascadeHashingFeatureMatcher::FetchHashedImage(
    const std::string& image_name) {
  const auto features = this->feature_and_matches_db_->GetFeatures(image_name);
  return std::make_shared<HashedImage>(
      cascade_hasher_->CreateHashedSiftDescriptors(features.descriptors));
}

// Initializes the cascade hasher (only if needed).
void CascadeHashingFeatureMatcher::InitializeCascadeHasher(
    int descriptor_dimension) {
  CHECK_GT(descriptor_dimension, 0);
  // Initialize the cascade hasher
  cascade_hasher_.reset(new CascadeHasher());
  CHECK(cascade_hasher_->Initialize(descriptor_dimension))
      << "Could not initialize the cascade hasher.";
}

void CascadeHashingFeatureMatcher::AddImage(const std::string& image_name) {
  FeatureMatcher::AddImage(image_name);

  // Only initialize the cascade hasher if it is not initialized.
  if (cascade_hasher_) {
    return;
  }

  // Get the features from the db and create hashed descriptors.
  const KeypointsAndDescriptors& features =
      this->feature_and_matches_db_->GetFeatures(image_name);

  if (features.descriptors.size() == 0) {
    return;
  }

  // Initialize the cascade hasher if needed.
  InitializeCascadeHasher(features.descriptors[0].size());
}

void CascadeHashingFeatureMatcher::AddImages(
    const std::vector<std::string>& image_names) {
  CHECK_GT(image_names.size(), 0);
  image_names_.reserve(image_names.size() + image_names_.size());
  image_names_.insert(
      image_names_.end(), image_names.begin(), image_names.end());

  // Initialize cascade hasher (if needed).
  for (int i = 0; i < image_names.size(); i++) {
    const KeypointsAndDescriptors& init_features =
        this->feature_and_matches_db_->GetFeatures(image_names[i]);
    if (init_features.descriptors.size() > 0) {
      InitializeCascadeHasher(init_features.descriptors[0].size());
      return;
    }
  }
}

bool CascadeHashingFeatureMatcher::MatchImagePair(
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    std::vector<IndexedFeatureMatch>* matches) {
  // Get pointers to the hashed images for each set of features.
  auto hashed_features1 = hashed_images_->Fetch(features1.image_name);
  auto hashed_features2 = hashed_images_->Fetch(features2.image_name);

  // If no hashed features exist for either image, skip.
  if (!hashed_features1 || !hashed_features2) {
    return false;
  }

  // Match features between the images.
  const double lowes_ratio =
      (this->options_.use_lowes_ratio) ? this->options_.lowes_ratio : 1.0;
  cascade_hasher_->MatchImages(*hashed_features1,
                               features1.descriptors,
                               *hashed_features2,
                               features2.descriptors,
                               lowes_ratio,
                               matches);
  // Only do symmetric matching if enough matches exist to begin with.
  if (matches->size() >= this->options_.min_num_feature_matches &&
      this->options_.keep_only_symmetric_matches) {
    std::vector<IndexedFeatureMatch> backwards_matches;
    cascade_hasher_->MatchImages(*hashed_features2,
                                 features2.descriptors,
                                 *hashed_features1,
                                 features1.descriptors,
                                 lowes_ratio,
                                 &backwards_matches);
    IntersectMatches(backwards_matches, matches);
  }

  return matches->size() >= this->options_.min_num_feature_matches;
}

}  // namespace theia
