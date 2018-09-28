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

#include "theia/matching/feature_matcher.h"

#include <glog/logging.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <mutex>  // NOLINT
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/io/read_keypoints_and_descriptors.h"
#include "theia/io/write_keypoints_and_descriptors.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/matching/fisher_vector_extractor.h"
#include "theia/matching/global_descriptor_extractor.h"
#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/two_view_match_geometric_verification.h"
#include "theia/util/filesystem.h"
#include "theia/util/lru_cache.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/util/util.h"

namespace theia {
namespace {
void SelectAllPairs(
    const std::vector<std::string>& image_names,
    std::vector<std::pair<std::string, std::string>>* pairs_to_match) {
  // Compute the total number of potential matches.
  const int num_pairs_to_match =
      image_names.size() * (image_names.size() - 1) / 2;

  pairs_to_match->reserve(num_pairs_to_match);
  // Create a list of all possible image pairs.
  for (int i = 0; i < image_names.size(); i++) {
    for (int j = i + 1; j < image_names.size(); j++) {
      pairs_to_match->emplace_back(image_names[i], image_names[j]);
    }
  }
}
}  // namespace

FeatureMatcher::~FeatureMatcher() {}

FeatureMatcher::FeatureMatcher(const FeatureMatcherOptions& options)
    : options_(options) {
  if (options_.match_out_of_core) {
    CHECK_GT(options_.cache_capacity, 2)
        << "The cache capacity must be greater than 2 in order to perform out "
           "of core matching.";
    // Determine if the directory for writing out feature exists. If not, try to
    // create it.
    if (!DirectoryExists(options_.keypoints_and_descriptors_output_dir)) {
      CHECK(CreateNewDirectory(options_.keypoints_and_descriptors_output_dir))
          << "Could not create the directory for storing features during "
             "matching: "
          << options_.keypoints_and_descriptors_output_dir;
    }
  } else {
    // If we want to perform all-in-memory matching then set the cache size to
    // the maximum.
    options_.cache_capacity = std::numeric_limits<int>::max();
  }

  // Because the function that defines how the cache fetches features from disk
  // is a member function, we need to bind it to this instance of FeatureMatcher
  // and specify that it will take in 1 argument.
  std::function<std::shared_ptr<KeypointsAndDescriptors>(const std::string&)>
      fetch_features_from_cache =
          std::bind(&FeatureMatcher::FetchKeypointsAndDescriptorsFromDisk,
                    this,
                    std::placeholders::_1);

  // Initialize the LRU cache. NOTE: even though the Fetch method will be set up
  // to retreive files from disk, it will only do so if
  // options_.match_out_of_core is set to true.
  keypoints_and_descriptors_cache_.reset(new KeypointAndDescriptorCache(
      fetch_features_from_cache, options_.cache_capacity));

  // Initialize the global image descriptor extractor if desired.
  if (options_.select_image_pairs_with_global_image_descriptor_matching) {
    FisherVectorExtractor::Options fv_options;
    fv_options.num_gmm_clusters = options_.num_gmm_clusters_for_fisher_vector;
    fv_options.max_num_features_for_training =
        options_.max_num_features_for_fisher_vector_training;
    global_image_descriptor_extractor_.reset(
        new FisherVectorExtractor(fv_options));
  }
}

void FeatureMatcher::AddImage(const std::string& image_name,
                              const std::vector<Keypoint>& keypoints,
                              const std::vector<Eigen::VectorXf>& descriptors) {
  image_names_.push_back(image_name);

  // Write the features file to disk.
  const std::string features_file = FeatureFilenameFromImage(image_name);
  if (options_.match_out_of_core) {
    CHECK(WriteKeypointsAndDescriptors(features_file, keypoints, descriptors))
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

  // Add the descriptors to the global image descriptor extractor for training
  // if using a global image descriptor extractor.
  if (options_.select_image_pairs_with_global_image_descriptor_matching) {
    CHECK_GT(descriptors.size(), 0);
    global_image_descriptor_extractor_->AddFeaturesForTraining(descriptors);
  }
}

void FeatureMatcher::AddImage(const std::string& image_name,
                              const std::vector<Keypoint>& keypoints,
                              const std::vector<Eigen::VectorXf>& descriptors,
                              const CameraIntrinsicsPrior& intrinsics) {
  AddImage(image_name, keypoints, descriptors);
  intrinsics_[image_name] = intrinsics;
}

void FeatureMatcher::AddImage(const std::string& image_name) {
  image_names_.push_back(image_name);

  // Add the descriptors to the global image descriptor extractor for training
  // if using a global image descriptor extractor.
  if (options_.select_image_pairs_with_global_image_descriptor_matching) {
    std::shared_ptr<KeypointsAndDescriptors> features =
        keypoints_and_descriptors_cache_->Fetch(
            FeatureFilenameFromImage(image_name));
    CHECK_GT(features->descriptors.size(), 0);
    global_image_descriptor_extractor_->AddFeaturesForTraining(
        features->descriptors);
  }
}

void FeatureMatcher::AddImage(const std::string& image_name,
                              const CameraIntrinsicsPrior& intrinsics) {
  AddImage(image_name);
  intrinsics_[image_name] = intrinsics;
}

void FeatureMatcher::AddImages(
    const std::vector<std::string>& image_names,
    const std::vector<CameraIntrinsicsPrior>& intrinsics) {
  CHECK_EQ(image_names.size(), intrinsics.size());
  image_names_.reserve(image_names.size() + image_names_.size());

  for (int i = 0; i < image_names.size(); ++i) {
    AddImage(image_names[i], intrinsics[i]);
  }
}

std::string FeatureMatcher::FeatureFilenameFromImage(const std::string& image) {
  std::string output_dir = options_.keypoints_and_descriptors_output_dir;
  // Add a trailing slash if one does not exist.
  if (output_dir.back() != '/') {
    output_dir = output_dir + "/";
  }
  return output_dir + image + ".features";
}

std::shared_ptr<KeypointsAndDescriptors>
FeatureMatcher::FetchKeypointsAndDescriptorsFromDisk(
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

void FeatureMatcher::SetImagePairsToMatch(
    const std::vector<std::pair<std::string, std::string>>& pairs_to_match) {
  pairs_to_match_ = pairs_to_match;
}

void FeatureMatcher::SelectImagePairsWithGlobalDescriptorMatching() {
  // Train the global descriptor extractor based on the input features.
  VLOG(2) << "Training global image descriptor...";
  CHECK(global_image_descriptor_extractor_->Train());

  // Extract the global descriptors.
  std::unique_ptr<ThreadPool> pool(new ThreadPool(options_.num_threads));
  std::vector<Eigen::VectorXf> global_descriptors(image_names_.size());
  for (int i = 0; i < image_names_.size(); i++) {
    pool->Add(
        [&](const int i) {
          std::shared_ptr<KeypointsAndDescriptors> features =
              keypoints_and_descriptors_cache_->Fetch(
                  FeatureFilenameFromImage(image_names_[i]));
          // Extract the global descriptors
          global_descriptors[i] =
              global_image_descriptor_extractor_->ExtractGlobalDescriptor(
                  features->descriptors);
        },
        i);
  }
  // Wait for threads to finish.
  pool.reset();

  VLOG(2) << "Computing image-to-image similarity scores with global "
             "descriptors...";
  // Match all pairs of global descriptors. The map contains the matching scores
  // for each image_names index.
  std::unordered_map<int, std::vector<std::pair<float, int>>>
      global_matching_scores;
  global_matching_scores.reserve(image_names_.size());
  for (int i = 0; i < image_names_.size(); i++) {
    for (int j = i + 1; j < image_names_.size(); j++) {
      const float global_feature_match_score =
          (global_descriptors[i] - global_descriptors[j]).squaredNorm();
      // Add the global feature matching score to both images.
      global_matching_scores[i].emplace_back(global_feature_match_score, j);
      global_matching_scores[j].emplace_back(global_feature_match_score, i);
    }
  }

  VLOG(2) << "Selecting the k-Nearest Neighbors images to perform feature "
             "matching on...";
  // For each image, find the kNN and add those to our selection for matching.
  // We use an unordered_set with indices instead of the string image names to
  // help with performance and avoid duplicate pairs.
  const int num_nearest_neighbors =
      std::min(static_cast<int>(image_names_.size() - 1),
               options_.num_nearest_neighbors_for_global_descriptor_matching);
  std::unordered_set<std::pair<int, int>> pairs_to_match_indices(
      num_nearest_neighbors * image_names_.size());
  for (auto& global_matching_results : global_matching_scores) {
    const int first_id = global_matching_results.first;
    auto& matching_results = global_matching_results.second;

    // Sort the results to get the top K neighbors for this image.
    std::partial_sort(matching_results.begin(),
                      matching_results.begin() + num_nearest_neighbors,
                      matching_results.end());

    // Add each of the kNN to the output indices.
    for (int i = 0; i < num_nearest_neighbors; i++) {
      const int second_id = matching_results[i].second;
      // Sort the pairwise index such that the pair of indices is sorted.
      const auto id_pair = (first_id < second_id)
                               ? std::make_pair(first_id, second_id)
                               : std::make_pair(second_id, first_id);
      pairs_to_match_indices.emplace(id_pair);
    }
  }

  // Convert to the output type.
  pairs_to_match_.reserve(pairs_to_match_indices.size());
  for (const auto& pair_to_match_indices : pairs_to_match_indices) {
    pairs_to_match_.emplace_back(image_names_[pair_to_match_indices.first],
                                 image_names_[pair_to_match_indices.second]);
  }
}

void FeatureMatcher::MatchImages(std::vector<ImagePairMatch>* matches) {
  // If SetImagePairsToMatch has not been called, match all image-to-image
  // pairs.
  if (pairs_to_match_.size() == 0) {
    if (options_.select_image_pairs_with_global_image_descriptor_matching) {
      SelectImagePairsWithGlobalDescriptorMatching();
    } else {
      SelectAllPairs(image_names_, &pairs_to_match_);
    }
  }

  // Add workers for matching. It is more efficient to let each thread compute
  // multiple matches at a time than add each matching task to the pool. This is
  // sort of like OpenMP's dynamic schedule in that it is able to balance
  // threads fairly efficiently.
  const int num_matches = pairs_to_match_.size();
  matches->reserve(num_matches);
  const int num_threads =
      std::min(options_.num_threads, static_cast<int>(num_matches));
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

void FeatureMatcher::MatchAndVerifyImagePairs(
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

    // Compute the visual matches from feature descriptors.
    std::vector<IndexedFeatureMatch> putative_matches;
    if (!MatchImagePair(*features1, *features2, &putative_matches)) {
      VLOG(2)
          << "Could not match a sufficient number of features between images "
          << image1_name << " and " << image2_name;
      continue;
    }

    // Perform geometric verification if applicable.
    if (options_.perform_geometric_verification) {
      // If geometric verification fails, do not add the match to the output.
      if (!GeometricVerification(
              *features1, *features2, putative_matches, &image_pair_match)) {
        VLOG(2) << "Geometric verification between images " << image1_name
                << " and " << image2_name << " failed.";
        continue;
      }
    } else {
      // If no geometric verification is performed then the putative matches are
      // output.
      image_pair_match.correspondences.reserve(putative_matches.size());
      for (int i = 0; i < putative_matches.size(); i++) {
        const Keypoint& keypoint1 =
            features1->keypoints[putative_matches[i].feature1_ind];
        const Keypoint& keypoint2 =
            features2->keypoints[putative_matches[i].feature2_ind];
        image_pair_match.correspondences.emplace_back(
            Feature(keypoint1.x(), keypoint1.y()),
            Feature(keypoint2.x(), keypoint2.y()));
      }
    }

    // Log information about the matching results.
    VLOG(1) << "Images " << image1_name << " and " << image2_name
            << " were matched with " << image_pair_match.correspondences.size()
            << " verified matches and "
            << image_pair_match.twoview_info.num_homography_inliers
            << " homography matches out of " << putative_matches.size()
            << " putative matches.";
    {
      std::lock_guard<std::mutex> lock(mutex_);
      matches->push_back(image_pair_match);
    }
  }
}

bool FeatureMatcher::GeometricVerification(
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    const std::vector<IndexedFeatureMatch>& putative_matches,
    ImagePairMatch* image_pair_match) {
  const CameraIntrinsicsPrior intrinsics1 = FindWithDefault(
      intrinsics_, features1.image_name, CameraIntrinsicsPrior());
  const CameraIntrinsicsPrior intrinsics2 = FindWithDefault(
      intrinsics_, features2.image_name, CameraIntrinsicsPrior());

  TwoViewMatchGeometricVerification geometric_verification(
      options_.geometric_verification_options,
      intrinsics1,
      intrinsics2,
      features1,
      features2,
      putative_matches);

  // Return whether geometric verification succeeds.
  return geometric_verification.VerifyMatches(
      &image_pair_match->correspondences, &image_pair_match->twoview_info);
}

}  // namespace theia
