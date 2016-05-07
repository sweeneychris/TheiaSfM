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
#include "theia/matching/indexed_feature_match.h"
#include "theia/util/lru_cache.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {
namespace {
// An LRU cache that will manage the keypoints and descriptors of interest.
typedef LRUCache<std::string, std::shared_ptr<KeypointsAndDescriptors> >
KeypointAndDescriptorCache;

// A unique pointer to a keypoint and descriptor cache.
typedef std::unique_ptr<KeypointAndDescriptorCache>
KeypointAndDescriptorCachePtr;

void CreateHashedImage(
    const std::string& image_name,
    const std::string& feature_filename,
    const std::unique_ptr<CascadeHasher>& cascade_hasher,
    const KeypointAndDescriptorCachePtr& kpts_and_descriptors_cache,
    std::mutex* hashed_images_mutex,
    std::unordered_map<std::string, HashedImage>* hashed_images) {
  // Get the features from the cache and create hashed descriptors.
  std::shared_ptr<KeypointsAndDescriptors> features =
      kpts_and_descriptors_cache->Fetch(feature_filename);
  // Create the hashing information.
  const HashedImage& hashed_image =
      cascade_hasher->CreateHashedSiftDescriptors(features->descriptors);
  {
    std::lock_guard<std::mutex> lock(*hashed_images_mutex);
    (*hashed_images)[image_name] = hashed_image;
    VLOG(1) << "Created the hashed descriptors for image: " << image_name;
  }
}

}  // namespace

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

void CascadeHashingFeatureMatcher::AddImage(
    const std::string& image,
    const std::vector<Keypoint>& keypoints,
    const std::vector<Eigen::VectorXf>& descriptors) {
  // This will save the descriptors and keypoints to disk and set up our LRU
  // cache.
  FeatureMatcher::AddImage(image, keypoints, descriptors);

  if (cascade_hasher_.get() == nullptr && descriptors.size() > 0) {
    cascade_hasher_.reset(new CascadeHasher());
    CHECK(cascade_hasher_->Initialize(descriptors[0].size()))
        << "Could not initialize the cascade hasher.";
  }

  // Create the hashing information.
  if (!ContainsKey(hashed_images_, image)) {
    hashed_images_[image] =
      cascade_hasher_->CreateHashedSiftDescriptors(descriptors);
    VLOG(1) << "Created the hashed descriptors for image: " << image;
  }
}

void CascadeHashingFeatureMatcher::AddImage(
    const std::string& image,
    const std::vector<Keypoint>& keypoints,
    const std::vector<Eigen::VectorXf>& descriptors,
    const CameraIntrinsicsPrior& intrinsics) {
  // This will save the descriptors and keypoints to disk and set up our LRU
  // cache.
  FeatureMatcher::AddImage(image, keypoints, descriptors, intrinsics);

  if (cascade_hasher_.get() == nullptr && descriptors.size() > 0) {
    cascade_hasher_.reset(new CascadeHasher());
    CHECK(cascade_hasher_->Initialize(descriptors[0].size()))
        << "Could not initialize the cascade hasher.";
  }

  // Create the hashing information.
  if (!ContainsKey(hashed_images_, image)) {
    hashed_images_[image] =
      cascade_hasher_->CreateHashedSiftDescriptors(descriptors);
    VLOG(1) << "Created the hashed descriptors for image: " << image;
  }
}

void CascadeHashingFeatureMatcher::AddImage(const std::string& image_name) {
  image_names_.push_back(image_name);

  // Get the features from the cache and create hashed descriptors.
  std::shared_ptr<KeypointsAndDescriptors> features =
      this->keypoints_and_descriptors_cache_->Fetch(
          FeatureFilenameFromImage(image_name));

  // Initialize the cascade hasher if needed.
  InitializeCascadeHasher(features->descriptors[0].size());

  // Create the hashing information.
  hashed_images_[image_name] =
      cascade_hasher_->CreateHashedSiftDescriptors(features->descriptors);
  VLOG(1) << "Created the hashed descriptors for image: " << image_name;
}

void CascadeHashingFeatureMatcher::AddImage(
    const std::string& image_name, const CameraIntrinsicsPrior& intrinsics) {
  image_names_.push_back(image_name);
  intrinsics_[image_name] = intrinsics;
  // Get the features from the cache and create hashed descriptors.
  std::shared_ptr<KeypointsAndDescriptors> features =
      this->keypoints_and_descriptors_cache_->Fetch(
          FeatureFilenameFromImage(image_name));

  // Initialize the cascade hasher if needed.
  InitializeCascadeHasher(features->descriptors[0].size());

  // Create the hashing information.
  hashed_images_[image_name] =
      cascade_hasher_->CreateHashedSiftDescriptors(features->descriptors);
  VLOG(1) << "Created the hashed descriptors for image: " << image_name;
}

void CascadeHashingFeatureMatcher::CreateHashedImagesInParallel(
    const std::vector<std::string>& image_names,
    const std::vector<std::string>& feature_filenames) {
  ThreadPool thread_pool(options_.num_threads);
  for (int i = 0; i < image_names.size(); ++i) {
    thread_pool.Add(CreateHashedImage,
                    std::cref(image_names[i]),
                    std::cref(feature_filenames[i]),
                    std::cref(cascade_hasher_),
                    std::cref(this->keypoints_and_descriptors_cache_),
                    &hashed_images_lock_,
                    &hashed_images_);
  }
}

void CascadeHashingFeatureMatcher::AddImages(
    const std::vector<std::string>& image_names,
    const std::vector<CameraIntrinsicsPrior>& intrinsics) {
  CHECK_EQ(image_names.size(), intrinsics.size())
      << "Number of images and intrinsic parameters mismatches.";
  image_names_.reserve(image_names.size() + image_names_.size());
  image_names_.insert(image_names_.end(),
                      image_names.begin(),
                      image_names.end());
  std::vector<std::string> feature_filenames(image_names.size());
  for (int i = 0; i < image_names.size(); ++i) {
    intrinsics_[image_names[i]] = intrinsics[i];
    feature_filenames[i] = FeatureFilenameFromImage(image_names[i]);
  }
  // Initialize cascade hasher (if needed).
  std::shared_ptr<KeypointsAndDescriptors> features =
      this->keypoints_and_descriptors_cache_->Fetch(feature_filenames[0]);
  InitializeCascadeHasher(features->descriptors[0].size());
  // Create the hashed images.
  CreateHashedImagesInParallel(image_names, feature_filenames);
}

bool CascadeHashingFeatureMatcher::MatchImagePair(
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    std::vector<FeatureCorrespondence>* matched_features) {
  const double lowes_ratio =
      (this->options_.use_lowes_ratio) ? this->options_.lowes_ratio : 1.0;

  // Get references to the hashed images for each set of features.
  HashedImage& hashed_features1 =
      FindOrDie(hashed_images_, features1.image_name);

  HashedImage& hashed_features2 =
      FindOrDie(hashed_images_, features2.image_name);

  std::vector<IndexedFeatureMatch> matches;
  cascade_hasher_->MatchImages(hashed_features1, features1.descriptors,
                               hashed_features2, features2.descriptors,
                               lowes_ratio, &matches);
  // Only do symmetric matching if enough matches exist to begin with.
  if (matches.size() >= this->options_.min_num_feature_matches &&
      this->options_.keep_only_symmetric_matches) {
    std::vector<IndexedFeatureMatch> backwards_matches;
    cascade_hasher_->MatchImages(hashed_features2,
                                 features2.descriptors,
                                 hashed_features1,
                                 features1.descriptors,
                                 lowes_ratio,
                                 &backwards_matches);
    IntersectMatches(backwards_matches, &matches);
  }

  if (matches.size() < this->options_.min_num_feature_matches) {
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
