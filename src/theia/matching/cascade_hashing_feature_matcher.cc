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
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {

// Initializes the cascade hasher (only if needed).
void CascadeHashingFeatureMatcher::InitializeCascadeHasher(
    int descriptor_dimension) {
  // Initialize the cascade hasher if needed.
  if (cascade_hasher_.get() == nullptr && descriptor_dimension > 0) {
    cascade_hasher_.reset(new CascadeHasher());
    CHECK(cascade_hasher_->Initialize(descriptor_dimension))
        << "Could not initialize the cascade hasher.";
  }
}

void CascadeHashingFeatureMatcher::AddImage(const std::string& image_name) {
  FeatureMatcher::AddImage(image_name);

  // Get the features from the db and create hashed descriptors.
  const KeypointsAndDescriptors& features =
      this->feature_and_matches_db_->GetFeatures(image_name);

  if (features.descriptors.size() == 0) {
    return;
  }

  // Initialize the cascade hasher if needed.
  InitializeCascadeHasher(features.descriptors[0].size());

  // Create the hashing information.
  hashed_images_[image_name] =
      cascade_hasher_->CreateHashedSiftDescriptors(features.descriptors);
  VLOG(1) << "Created the hashed descriptors for image: " << image_name;
}

void CascadeHashingFeatureMatcher::AddImage(
    const std::string& image_name, const CameraIntrinsicsPrior& intrinsics) {
  FeatureMatcher::AddImage(image_name, intrinsics);

  // Get the features from the cache and create hashed descriptors.
  const KeypointsAndDescriptors& features =
      this->feature_and_matches_db_->GetFeatures(image_name);

  if (features.descriptors.size() == 0) {
    return;
  }

  // Initialize the cascade hasher if needed.
  InitializeCascadeHasher(features.descriptors[0].size());

  // Create the hashing information.
  hashed_images_[image_name] =
      cascade_hasher_->CreateHashedSiftDescriptors(features.descriptors);
  VLOG(1) << "Created the hashed descriptors for image: " << image_name;
}

void CascadeHashingFeatureMatcher::AddImages(
    const std::vector<std::string>& image_names,
    const std::vector<CameraIntrinsicsPrior>& intrinsics) {
  CHECK_EQ(image_names.size(), intrinsics.size())
      << "Number of images and intrinsic parameters mismatches.";
  image_names_.reserve(image_names.size() + image_names_.size());
  image_names_.insert(
      image_names_.end(), image_names.begin(), image_names.end());
  for (int i = 0; i < image_names.size(); ++i) {
    intrinsics_[image_names[i]] = intrinsics[i];
  }

  // Initialize cascade hasher (if needed).
  const KeypointsAndDescriptors& init_features =
      this->feature_and_matches_db_->GetFeatures(image_names[0]);
  InitializeCascadeHasher(init_features.descriptors[0].size());

  // Create the hashed images.
  ThreadPool thread_pool(options_.num_threads);
  for (int i = 0; i < image_names.size(); ++i) {
    thread_pool.Add(
        [&](const int i) {
          // Get the features from the cache and create hashed descriptors.
          const KeypointsAndDescriptors& features =
              this->feature_and_matches_db_->GetFeatures(image_names[i]);
          // Create the hashing information.
          const HashedImage hashed_image =
              cascade_hasher_->CreateHashedSiftDescriptors(
                  features.descriptors);

          std::lock_guard<std::mutex> lock(hashed_images_lock_);
          hashed_images_[image_names[i]] = hashed_image;
          VLOG(1) << "Created the hashed descriptors for image: "
                  << image_names[i];
        },
        i);
  }
}

bool CascadeHashingFeatureMatcher::MatchImagePair(
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    std::vector<IndexedFeatureMatch>* matches) {
  // Get pointers to the hashed images for each set of features.
  HashedImage* hashed_features1 =
      FindOrNull(hashed_images_, features1.image_name);

  HashedImage* hashed_features2 =
      FindOrNull(hashed_images_, features2.image_name);

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
