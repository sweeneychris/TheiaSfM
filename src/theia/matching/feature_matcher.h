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

#include <glog/logging.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <mutex>  // NOLINT
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/io/read_keypoints_and_descriptors.h"
#include "theia/io/write_keypoints_and_descriptors.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/matching/image_pair_match.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/verify_two_view_matches.h"
#include "theia/util/filesystem.h"
#include "theia/util/lru_cache.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {

// This struct is used by the internal cache to hold keypoints and descriptors
// when the are retrieved from the cache.
struct KeypointsAndDescriptors {
  std::string image_name;
  std::vector<Keypoint> keypoints;
  std::vector<Eigen::VectorXf> descriptors;
};

// Class for matching features between images. The intended use for these
// classes is for matching photos in image collections, so all pairwise matches
// are computed. Matching with geometric verification is also possible. Typical
// use case is:
//   FeatureMatcherOptions matcher_options;
//   FeatureMatcher matcher(matcher_options);
//   for (int i = 0; i < num_images_to_match; i++) {
//     matcher.AddImage(image_names[i], keypoints[i], descriptors[i]);
//     // Or, you could add the image with known intrinsics for use during
//     // geometric verification.
//     matcher.AddImage(image_names[i], keypoints[i],
//                      descriptors[i], intrinsics[i]);
//   }
//   std::vector<ImagePairMatch> matches;
//   matcher.MatchImages(&matches);
//   // Or, with geometric verification:
//   VerifyTwoViewMatchesOptions geometric_verification_options;
//   matcher.MatchImages(geometric_verification_options, &matches);
//
// The matches and match quality depend on the options passed to the feature
// matching.
template <class DistanceMetric> class FeatureMatcher {
 public:
  typedef typename DistanceMetric::DistanceType DistanceType;

  explicit FeatureMatcher(const FeatureMatcherOptions& matcher_options);
  virtual ~FeatureMatcher() {}

  // Adds an image to the matcher with no known intrinsics for this image. The
  // caller still owns the keypoints and descriptors so they must remain valid
  // objects throughout the matching. The image name must be a unique identifier
  // for the image.
  virtual void AddImage(const std::string& image_name,
                        const std::vector<Keypoint>& keypoints,
                        const std::vector<Eigen::VectorXf>& descriptors);

  // Adds an image to the matcher with the known camera intrinsics. The
  // intrinsics (if known) are useful for geometric verification. The caller
  // still owns the keypoints and descriptors so they must remain valid objects
  // throughout the matching. The image name must be a unique identifier for the
  // image.
  virtual void AddImage(const std::string& image_name,
                        const std::vector<Keypoint>& keypoints,
                        const std::vector<Eigen::VectorXf>& descriptors,
                        const CameraIntrinsicsPrior& intrinsics);

  // If features have been written to disk, the matcher can directly work with
  // them from the feature files so that you do not have to "add" them to the
  // matcher. This assumes that feature files have been written in the format:
  //
  //     //keypoints_and_descriptors_output_dir/image_name.features
  //
  // Care must be taken to ensure that the filenames are formatted correctly.
  virtual void AddImage(const std::string& image_name);
  virtual void AddImage(const std::string& image_name,
                        const CameraIntrinsicsPrior& intrinsics);

  // Matches features between all images. No geometric verification is
  // performed. Only the matches which pass the have greater than
  // min_num_feature_matches are returned.
  virtual void MatchImages(std::vector<ImagePairMatch>* matches);

  // Matches features between all images. Only the matches that pass the
  // geometric verification are returned. Camera intrinsics are used for
  // geometric verification if the image was added with known intrinsics.
  virtual void MatchImagesWithGeometricVerification(
      const VerifyTwoViewMatchesOptions& verification_options,
      std::vector<ImagePairMatch>* matches);

  // Set the image pairs that will be matched when MatchImages or
  // MatchImagesWithGeometricVerification is called. This is an optional method;
  // if it is not called, then all possible image-to-image pairs will be
  // matched. The vector should contain unique pairs of image names that should
  // be matched.
  virtual void SetImagePairsToMatch(
      const std::vector<std::pair<std::string, std::string> >& pairs_to_match);

 protected:
  // NOTE: This method should be overridden in the subclass implementations!
  // Returns true if the image pair is a valid match.
  virtual bool MatchImagePair(
      const KeypointsAndDescriptors& features1,
      const KeypointsAndDescriptors& features2,
      std::vector<FeatureCorrespondence>* matched_features) = 0;

  // Performs matching and geometric verification (if desired) on the
  // pairs_to_match_ between the specified indices. This is useful for thread
  // pooling.
  virtual void MatchAndVerifyImagePairs(const int start_index,
                                        const int end_index,
                                        std::vector<ImagePairMatch>* matches);

  // Fetches keypoints and descriptors from disk. This function is utilized by
  // the internal cache to preserve memory.
  static std::shared_ptr<KeypointsAndDescriptors>
  FetchKeypointsAndDescriptorsFromDisk(const std::string& features_file);

  // Returns the filepath of the feature file given the image name.
  std::string FeatureFilenameFromImage(const std::string& image);

  // Each Threadpool worker will perform matching on this many image pairs.  It
  // is more efficient to let each thread compute multiple matches at a time
  // than add each matching task to the pool. This is sort of like OpenMP's
  // dynamic schedule in that it is able to balance threads fairly efficiently.
  const int kMaxThreadingStepSize_ = 20;

  FeatureMatcherOptions matcher_options_;
  VerifyTwoViewMatchesOptions verification_options_;
  // Will be set to true if geometric verification is enabled.
  bool verify_image_pairs_;

  // A container for the image names.
  std::vector<std::string> image_names_;

  // An LRU cache that will manage the keypoints and descriptors of interest.
  typedef LRUCache<std::string, std::shared_ptr<KeypointsAndDescriptors> >
      KeypointAndDescriptorCache;
  std::unique_ptr<KeypointAndDescriptorCache> keypoints_and_descriptors_cache_;

  std::unordered_map<std::string, CameraIntrinsicsPrior> intrinsics_;
  std::vector<std::pair<std::string, std::string> > pairs_to_match_;
  std::mutex mutex_;

 private:
  DISALLOW_COPY_AND_ASSIGN(FeatureMatcher);
};

// ---------------------- Implementation ------------------------ //

template <class DistanceMetric>
FeatureMatcher<DistanceMetric>::FeatureMatcher(
    const FeatureMatcherOptions& options)
    : matcher_options_(options), verify_image_pairs_(true) {
  // Initialize the LRU cache. NOTE: even though the Fetch method will be set up
  // to retreive files from disk, it will only do so if
  // matcher_options_.match_out_of_core is set to true.
  keypoints_and_descriptors_cache_.reset(new KeypointAndDescriptorCache(
      &FeatureMatcher<DistanceMetric>::FetchKeypointsAndDescriptorsFromDisk,
      matcher_options_.cache_capacity));

  if (matcher_options_.match_out_of_core) {
    CHECK_GT(matcher_options_.cache_capacity, 2)
        << "The cache capacity must be greater than 2 in order to perform out "
           "of core matching.";
    // Determine if the directory for writing out feature exists. If not, try to
    // create it.
    if (!DirectoryExists(
            matcher_options_.keypoints_and_descriptors_output_dir)) {
      CHECK(CreateNewDirectory(
          matcher_options_.keypoints_and_descriptors_output_dir))
          << "Could not create the directory for storing features during "
             "matching: "
          << matcher_options_.keypoints_and_descriptors_output_dir;
    }
  } else {
    // If we want to perform all-in-memory matching then set the cache size to
    // the maximum.
    matcher_options_.cache_capacity = std::numeric_limits<int>::max();
  }
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::AddImage(
    const std::string& image_name,
    const std::vector<Keypoint>& keypoints,
    const std::vector<Eigen::VectorXf>& descriptors) {
  image_names_.push_back(image_name);

  // Write the features file to disk.
  const std::string features_file = FeatureFilenameFromImage(image_name);
  if (matcher_options_.match_out_of_core) {
    CHECK(WriteKeypointsAndDescriptors(features_file,
                                       keypoints,
                                       descriptors))
      << "Could not read features for image " << image_name << " from file "
      << features_file;
  }

  // Insert the features into the cache.
  std::shared_ptr<KeypointsAndDescriptors> keypoints_and_descriptors(
      new KeypointsAndDescriptors);
  keypoints_and_descriptors->image_name = image_name;
  keypoints_and_descriptors->keypoints = keypoints;
  keypoints_and_descriptors->descriptors = descriptors;
  keypoints_and_descriptors_cache_->Insert(features_file,
                                           keypoints_and_descriptors);
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::AddImage(
    const std::string& image_name,
    const std::vector<Keypoint>& keypoints,
    const std::vector<Eigen::VectorXf>& descriptors,
    const CameraIntrinsicsPrior& intrinsics) {
  AddImage(image_name, keypoints, descriptors);
  intrinsics_[image_name] = intrinsics;
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::AddImage(const std::string& image_name) {
  image_names_.push_back(image_name);
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::AddImage(
    const std::string& image_name, const CameraIntrinsicsPrior& intrinsics) {
  image_names_.push_back(image_name);
  intrinsics_[image_name] = intrinsics;
}

template <class DistanceMetric>
std::string FeatureMatcher<DistanceMetric>::FeatureFilenameFromImage(
    const std::string& image) {
  std::string output_dir =
      matcher_options_.keypoints_and_descriptors_output_dir;
  // Add a trailing slash if one does not exist.
  if (output_dir.back() != '/') {
    output_dir = output_dir + "/";
  }
  return output_dir + image + ".features";
}

template <class DistanceMetric>
std::shared_ptr<KeypointsAndDescriptors>
FeatureMatcher<DistanceMetric>::FetchKeypointsAndDescriptorsFromDisk(
    const std::string& features_file) {
  std::shared_ptr<KeypointsAndDescriptors> keypoints_and_descriptors(
      new KeypointsAndDescriptors);

  // Read in the features file from disk.
  CHECK(ReadKeypointsAndDescriptors(features_file,
                                    &keypoints_and_descriptors->keypoints,
                                    &keypoints_and_descriptors->descriptors))
      << "Could not read features from file " << features_file;
  return keypoints_and_descriptors;
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::SetImagePairsToMatch(
    const std::vector<std::pair<std::string, std::string> >& pairs_to_match) {
  pairs_to_match_ = pairs_to_match;
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::MatchImages(
    std::vector<ImagePairMatch>* matches) {
  // Set image verification to false so that it will be skipped.
  verify_image_pairs_ = false;
  VerifyTwoViewMatchesOptions verification_options;
  MatchImagesWithGeometricVerification(verification_options, matches);
  // Reset the value to true.
  verify_image_pairs_ = true;
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::MatchImagesWithGeometricVerification(
    const VerifyTwoViewMatchesOptions& verification_options,
    std::vector<ImagePairMatch>* matches) {
  verification_options_ = verification_options;

  // If SetImagePairsToMatch has not been called, match all image-to-image
  // pairs.
  if (pairs_to_match_.size() == 0) {
    // Compute the total number of potential matches.
    const int num_pairs_to_match =
      image_names_.size() * (image_names_.size() - 1) / 2;
    matches->reserve(num_pairs_to_match);

    pairs_to_match_.reserve(num_pairs_to_match);
    // Create a list of all possible image pairs.
    for (int i = 0; i < image_names_.size(); i++) {
      for (int j = i + 1; j < image_names_.size(); j++) {
        pairs_to_match_.emplace_back(image_names_[i], image_names_[j]);
      }
    }
  }

  // Add workers for matching. It is more efficient to let each thread compute
  // multiple matches at a time than add each matching task to the pool. This is
  // sort of like OpenMP's dynamic schedule in that it is able to balance
  // threads fairly efficiently.
  const int num_matches = pairs_to_match_.size();
  const int num_threads =
      std::min(matcher_options_.num_threads, static_cast<int>(num_matches));
  std::unique_ptr<ThreadPool> pool(new ThreadPool(num_threads));
  const int interval_step =
      std::min(this->kMaxThreadingStepSize_, num_matches / num_threads);
  for (int i = 0; i < num_matches; i += interval_step) {
    const int end_interval = std::min(num_matches, i + interval_step);
    pool->Add(&FeatureMatcher::MatchAndVerifyImagePairs,
              this,
              i,
              end_interval,
              matches);
  }
  // Wait for all threads to finish.
  pool.reset(nullptr);

  VLOG(1) << "Matched " << matches->size() << " image pairs out of "
          << num_matches << " possible image pairs.";
}

template <class DistanceMetric>
void FeatureMatcher<DistanceMetric>::MatchAndVerifyImagePairs(
    const int start_index,
    const int end_index,
    std::vector<ImagePairMatch>* matches) {
  for (int i = start_index; i < end_index; i++) {
    const std::string image1_name = pairs_to_match_[i].first;
    const std::string image2_name = pairs_to_match_[i].second;

    // Match the image pair. If the pair fails to match then continue to the
    // next match.
    ImagePairMatch image_pair_match;
    image_pair_match.image1 = image1_name;
    image_pair_match.image2 = image2_name;

    // Get the keypoints and descriptors from the cache. By using a shared_ptr
    // here we ensure that keypoints and descriptors will live in the cache as
    // long as they are currently being used in a matching thread, so the cache
    // will never evict these entries while they are still being used.
    std::shared_ptr<KeypointsAndDescriptors> features1 =
        keypoints_and_descriptors_cache_->Fetch(
            FeatureFilenameFromImage(image1_name));
    features1->image_name = image1_name;
    std::shared_ptr<KeypointsAndDescriptors> features2 =
        keypoints_and_descriptors_cache_->Fetch(
            FeatureFilenameFromImage(image2_name));
    features2->image_name = image2_name;

    if (!MatchImagePair(*features1,
                        *features2,
                        &image_pair_match.correspondences)) {
      VLOG(2)
          << "Could not match a sufficient number of features between images "
          << image1_name << " and " << image2_name;
      continue;
    }

    // Add images to the valid matches if no geometric verification is required.
    if (!verify_image_pairs_) {
      VLOG(1) << image_pair_match.correspondences.size()
              << " putative matches between images " << image1_name << " and "
              << image2_name;
      std::lock_guard<std::mutex> lock(mutex_);
      matches->push_back(image_pair_match);
      continue;
    }

    const CameraIntrinsicsPrior intrinsics1 =
        FindWithDefault(intrinsics_, image1_name, CameraIntrinsicsPrior());
    const CameraIntrinsicsPrior intrinsics2 =
        FindWithDefault(intrinsics_, image2_name, CameraIntrinsicsPrior());
    std::vector<int> inliers;
    // Do not add this image pair as a verified match if the verification does
    // not pass.
    if (!VerifyTwoViewMatches(verification_options_,
                              intrinsics1,
                              intrinsics2,
                              image_pair_match.correspondences,
                              &image_pair_match.twoview_info,
                              &inliers)) {
      VLOG(2) << "Geometric verification between images " << image1_name
              << " and " << image2_name << " failed.";
      continue;
    }

    // Output only the inliers.
    const std::vector<FeatureCorrespondence> old_correspondences =
        std::move(image_pair_match.correspondences);
    image_pair_match.correspondences.reserve(inliers.size());
    for (int j = 0; j < inliers.size(); ++j) {
      image_pair_match.correspondences.emplace_back(
          old_correspondences[inliers[j]]);
    }
    VLOG(1) << "Images " << image1_name << " and " << image2_name
            << " were matched with " << inliers.size()
            << " verified matches and "
            << image_pair_match.twoview_info.num_homography_inliers
            << " homography matches out of " << old_correspondences.size()
            << " putative matches.";
    {
      std::lock_guard<std::mutex> lock(mutex_);
      matches->push_back(image_pair_match);
    }
  }
}

}  // namespace theia

#endif  // THEIA_MATCHING_FEATURE_MATCHER_H_
