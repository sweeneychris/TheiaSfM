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

#ifndef THEIA_MATCHING_FEATURE_MATCHER_H_
#define THEIA_MATCHING_FEATURE_MATCHER_H_

#include <algorithm>
#include <glog/logging.h>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/matching/image_pair_match.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/verify_two_view_matches.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {

// Class for matching features between images. The intended use for these
// classes is for matching photos in image collections, so all pairwise matches
// are computed. Matching with geometric verification is also possible. Typical
// use case is:
//   FeatureMatcher matcher;
//   for (int i = 0; i < num_images_to_match; i++) {
//     matcher.AddImage(keypoints[i], descriptors[i]);
//     // Or, you could add the image with known intrinsics for use during
//     // geometric verification.
//     matcher.AddImage(keypoints[i], descriptors[i], intrinsics[i]);
//   }
//   std::vector<ImagePairMatch> matches;
//   FeatureMatcherOptions matcher_options;
//   matcher.MatchImages(matcher_options, &matches);
//   // Or, with geometric verification:
//   VerifyTwoViewMatchesOptions geometric_verification_options;
//   matcher.MatchImages(match_options,
//                       geometric_verification_options,
//                       &matches);
//
// The matches and match quality depend on the options passed to the feature
// matching.
template <class DistanceMetric> class FeatureMatcher {
 public:
  typedef typename DistanceMetric::DistanceType DistanceType;
  typedef typename DistanceMetric::DescriptorType DescriptorType;

  FeatureMatcher() : verify_image_pairs_(true) {}
  virtual ~FeatureMatcher() {}

  // Adds an image to the matcher with no known intrinsics for this image. The
  // caller still owns the keypoints and descriptors so they must remain valid
  // objects throughout the matching.
  virtual void AddImage(const std::vector<Keypoint>* keypoints,
                        const std::vector<DescriptorType>* descriptors);

  // Adds an image to the matcher with the known camera intrinsics. The
  // intrinsics (if known) are useful for geometric verification. The caller
  // still owns the keypoints and descriptors so they must remain valid objects
  // throughout the matching.
  virtual void AddImage(const std::vector<Keypoint>* keypoints,
                        const std::vector<DescriptorType>* descriptors,
                        const CameraIntrinsicsPrior& intrinsics);

  // Matches features between all images. No geometric verification is
  // performed. Only the matches which pass the have greater than
  // min_num_feature_matches are returned.
  virtual void MatchImages(const FeatureMatcherOptions& matcher_options,
                           std::vector<ImagePairMatch>* matches);

  // Matches features between all images. Only the matches that pass the
  // geometric verification are returned. Camera intrinsics are used for
  // geometric verification if the image was added with known intrinsics.
  virtual void MatchImagesWithGeometricVerification(
      const FeatureMatcherOptions& matcher_options,
      const VerifyTwoViewMatchesOptions& verification_options,
      std::vector<ImagePairMatch>* matches);

 protected:
  // NOTE: This method should be overridden in the subclass implementations!
  // Returns true if the image pair is a valid match.
  virtual bool MatchImagePair(
      const int image1_index,
      const int image2_index,
      std::vector<FeatureCorrespondence>* matched_features) = 0;

  // Performs matching and geometric verification (if desired) on the
  // pairs_to_match_ between the specified indices. This is useful for thread
  // pooling.
  virtual void MatchAndVerifyImagePairs(const int start_index,
                                        const int end_index,
                                        std::vector<ImagePairMatch>* matches);

  // Each Threadpool worker will perform matching on this many image pairs.  It
  // is more efficient to let each thread compute multiple matches at a time
  // than add each matching task to the pool. This is sort of like OpenMP's
  // dynamic schedule in that it is able to balance threads fairly efficiently.
  const int kMaxThreadingStepSize_ = 20;

  FeatureMatcherOptions matcher_options_;
  VerifyTwoViewMatchesOptions verification_options_;
  // Will be set to true if geometric verification is enabled.
  bool verify_image_pairs_;

  std::vector<const std::vector<Keypoint>*> keypoints_;
  std::vector<const std::vector<DescriptorType>*> descriptors_;
  std::unordered_map<int, CameraIntrinsicsPrior> intrinsics_;
  std::vector<std::pair<int, int> > pairs_to_match_;
  std::mutex mutex_;

 private:
  DISALLOW_COPY_AND_ASSIGN(FeatureMatcher);
};

// ---------------------- Implementation ------------------------ //

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::AddImage(
    const std::vector<Keypoint>* keypoints,
    const std::vector<DescriptorType>* descriptors) {
  CHECK_NOTNULL(keypoints);
  CHECK_NOTNULL(descriptors);

  keypoints_.push_back(keypoints);
  descriptors_.push_back(descriptors);
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::AddImage(
    const std::vector<Keypoint>* keypoints,
    const std::vector<DescriptorType>* descriptors,
    const CameraIntrinsicsPrior& intrinsics) {
  CHECK_NOTNULL(keypoints);
  CHECK_NOTNULL(descriptors);

  keypoints_.push_back(keypoints);
  descriptors_.push_back(descriptors);
  const int image_index = keypoints_.size() - 1;
  intrinsics_[image_index] = intrinsics;
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::MatchImages(
    const FeatureMatcherOptions& matcher_options,
    std::vector<ImagePairMatch>* matches) {
  // Set image verification to false so that it will be skipped.
  verify_image_pairs_ = false;
  VerifyTwoViewMatchesOptions verification_options;
  MatchImagesWithGeometricVerification(matcher_options,
                                       verification_options,
                                       matches);
  // Reset the value to true.
  verify_image_pairs_ = true;
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::MatchImagesWithGeometricVerification(
      const FeatureMatcherOptions& matcher_options,
      const VerifyTwoViewMatchesOptions& verification_options,
      std::vector<ImagePairMatch>* matches) {
  pairs_to_match_.clear();
  matcher_options_ = matcher_options;
  verification_options_ = verification_options;

  // Compute the total number of potential matches.
  const int num_pairs_to_match =
      keypoints_.size() * (keypoints_.size() - 1) / 2;
  matches->reserve(num_pairs_to_match);

  pairs_to_match_.reserve(num_pairs_to_match);

  // Create a list of all possible image pairs.
  for (int i = 0; i < keypoints_.size(); i++) {
    if (keypoints_[i]->size() == 0) {
      continue;
    }
    for (int j = i + 1; j < keypoints_.size(); j++) {
      if (keypoints_[j]->size() == 0) {
        continue;
      }
      pairs_to_match_.emplace_back(i, j);
    }
  }

  // Add workers for matchine. It is more efficient to let each thread compute
  // multiple matches at a time than add each matching task to the pool. This is
  // sort of like OpenMP's dynamic schedule in that it is able to balance
  // threads fairly efficiently.
  const int num_threads =
      std::min(matcher_options_.num_threads, num_pairs_to_match);
  std::unique_ptr<ThreadPool> pool(new ThreadPool(num_threads));
  const int interval_step =
      std::min(this->kMaxThreadingStepSize_, num_pairs_to_match / num_threads);
  for (int i = 0; i < pairs_to_match_.size(); i += interval_step) {
    const int end_interval = std::min(num_pairs_to_match, i + interval_step);
    pool->Add(&FeatureMatcher::MatchAndVerifyImagePairs,
              this,
              i,
              end_interval,
              matches);
  }
  // Wait for all threads to finish.
  pool.reset(nullptr);

  VLOG(1) << "Matched " << matches->size() << " image pairs out of "
          << num_pairs_to_match << " possible image pairs.";
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::MatchAndVerifyImagePairs(
    const int start_index,
    const int end_index,
    std::vector<ImagePairMatch>* matches) {
  for (int i = start_index; i < end_index; i++) {
    // Match the image pair. If the pair fails to match then continue to the
    // next match.
    ImagePairMatch image_pair_match;
    const int image1_index = pairs_to_match_[i].first;
    const int image2_index = pairs_to_match_[i].second;

    image_pair_match.image1_index = image1_index;
    image_pair_match.image2_index = image2_index;
    if (!MatchImagePair(image1_index,
                        image2_index,
                        &image_pair_match.correspondences)) {
      VLOG(2)
          << "Could not match a sufficient number of features between images "
          << image1_index << " and " << image2_index;
      continue;
    }

    // Add images to the valid matches if no geometric verification is required.
    if (!verify_image_pairs_) {
      VLOG(1) << image_pair_match.correspondences.size()
              << " putative matches between images " << image1_index << " and "
              << image2_index;
      std::lock_guard<std::mutex> lock(mutex_);
      matches->push_back(image_pair_match);
      continue;
    }


    const CameraIntrinsicsPrior intrinsics1 = FindWithDefault(
        intrinsics_, image1_index, CameraIntrinsicsPrior());
    const CameraIntrinsicsPrior intrinsics2 = FindWithDefault(
        intrinsics_, image2_index, CameraIntrinsicsPrior());
    // If the image pair passes two view verification then
    std::vector<int> inliers;
    // Do not add this image pair as a verified match if the verification does
    // not pass.
    if (!VerifyTwoViewMatches(verification_options_,
                              intrinsics1,
                              intrinsics2,
                              image_pair_match.correspondences,
                              &image_pair_match.twoview_info,
                              &inliers)) {
      VLOG(2) << "Geometric verification between images " << image1_index
              << " and " << image2_index << " failed.";
      continue;
    }

    // Output only the inliers.
    const std::vector<FeatureCorrespondence> old_correspondences =
        image_pair_match.correspondences;
    image_pair_match.correspondences.clear();
    image_pair_match.correspondences.reserve(inliers.size());
    for (int j = 0; j < inliers.size(); j++) {
      image_pair_match.correspondences.emplace_back(
          old_correspondences[inliers[j]]);
    }
    VLOG(1) << "Images " << image1_index << " and " << image2_index
            << " were matched with " << inliers.size()
            << " verified matches out of " << old_correspondences.size()
            << " putative matches.";
    {
      std::lock_guard<std::mutex> lock(mutex_);
      matches->push_back(image_pair_match);
    }
  }
}

}  // namespace theia

#endif  // THEIA_MATCHING_FEATURE_MATCHER_H_
