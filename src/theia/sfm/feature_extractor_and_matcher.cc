// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#include "theia/sfm/feature_extractor_and_matcher.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <algorithm>
#include <string>
#include <thread>  // NOLINT
#include <vector>

#include "theia/image/image.h"
#include "theia/image/descriptor/create_descriptor_extractor.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/create_feature_matcher.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/matching/image_pair_match.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/exif_reader.h"
#include "theia/sfm/verify_two_view_matches.h"
#include "theia/util/filesystem.h"
#include "theia/util/threadpool.h"

namespace theia {
namespace {

void ExtractFeatures(
    const FeatureExtractorAndMatcher::Options& options,
    const std::string& image_filepath,
    std::vector<Keypoint>* keypoints,
    std::vector<Eigen::VectorXf>* descriptors) {
  std::unique_ptr<FloatImage> image(new FloatImage(image_filepath));

  // We create these variable here instead of upon the construction of the
  // object so that they can be thread-safe. We *should* be able to use the
  // static thread_local keywords, but apparently Mac OS-X's version of clang
  // does not actually support it!
  //
  // TODO(cmsweeney): Change this so that each thread in the threadpool receives
  // exactly one object.
  CreateDescriptorExtractorOptions descriptor_options;
  descriptor_options.descriptor_extractor_type =
      options.descriptor_extractor_type;
  descriptor_options.sift_options = options.sift_parameters;

  std::unique_ptr<DescriptorExtractor> descriptor_extractor =
      CreateDescriptorExtractor(descriptor_options);

  // Exit if the descriptor extraction fails.
  if (!descriptor_extractor->DetectAndExtractDescriptors(*image,
                                                         keypoints,
                                                         descriptors)) {
    LOG(ERROR) << "Could not extract descriptors in image " << image_filepath;
    return;
  }

  if (keypoints->size() > options.max_num_features) {
    keypoints->resize(options.max_num_features);
    descriptors->resize(options.max_num_features);
  }

  VLOG(1) << "Successfully extracted " << descriptors->size()
          << " features from image " << image_filepath;
}

}  // namespace

FeatureExtractorAndMatcher::FeatureExtractorAndMatcher(
    const FeatureExtractorAndMatcher::Options& options)
    : options_(options) {
  // Create the feature matcher.
  FeatureMatcherOptions matcher_options = options_.feature_matcher_options;
  matcher_options.num_threads = options_.num_threads;
  matcher_options.min_num_feature_matches = options_.min_num_inlier_matches;
  matcher_options.perform_geometric_verification = true;
  matcher_options.geometric_verification_options =
      options.geometric_verification_options;
  matcher_options.geometric_verification_options.min_num_inlier_matches =
      options_.min_num_inlier_matches;

  matcher_ = CreateFeatureMatcher(options_.matching_strategy, matcher_options);
}

bool FeatureExtractorAndMatcher::AddImage(const std::string& image_filepath) {
  image_filepaths_.emplace_back(image_filepath);
  return true;
}

bool FeatureExtractorAndMatcher::AddImage(
    const std::string& image_filepath,
    const CameraIntrinsicsPrior& intrinsics) {
  if (!AddImage(image_filepath)) {
    return false;
  }
  intrinsics_[image_filepath] = intrinsics;
  return true;
}

// Performs feature matching between all images provided by the image
// filepaths. Features are extracted and matched between the images according to
// the options passed in. Only matches that have passed geometric verification
// are kept. EXIF data is parsed to determine the camera intrinsics if
// available.
void FeatureExtractorAndMatcher::ExtractAndMatchFeatures(
    std::vector<CameraIntrinsicsPrior>* intrinsics,
    std::vector<ImagePairMatch>* matches) {
  CHECK_NOTNULL(intrinsics)->resize(image_filepaths_.size());
  CHECK_NOTNULL(matches);
  CHECK_NOTNULL(matcher_.get());

  // For each image, process the features and add it to the matcher.
  const int num_threads =
      std::min(options_.num_threads, static_cast<int>(image_filepaths_.size()));
  std::unique_ptr<ThreadPool> thread_pool(new ThreadPool(num_threads));
  for (int i = 0; i < image_filepaths_.size(); i++) {
    if (!FileExists(image_filepaths_[i])) {
      LOG(ERROR) << "Could not extract features for " << image_filepaths_[i]
                 << " because the file cannot be found.";
      continue;
    }
    thread_pool->Add(&FeatureExtractorAndMatcher::ProcessImage, this, i);
  }
  // This forces all tasks to complete before proceeding.
  thread_pool.reset(nullptr);

  // After all threads complete feature extraction, perform matching.

  // Perform the matching.
  LOG(INFO) << "Matching images...";
  matcher_->MatchImages(matches);

  // Add the intrinsics to the output.
  for (int i = 0; i < image_filepaths_.size(); i++) {
    (*intrinsics)[i] = FindOrDie(intrinsics_, image_filepaths_[i]);
  }
}

void FeatureExtractorAndMatcher::ProcessImage(
    const int i) {
  const std::string& image_filepath = image_filepaths_[i];

  // Extract Exif if it wasn't provided.
  if (!ContainsKey(intrinsics_, image_filepath)) {
    CameraIntrinsicsPrior intrinsics;
    CHECK(exif_reader_.ExtractEXIFMetadata(image_filepath, &intrinsics));
    intrinsics_mutex_.lock();
    intrinsics_.emplace(image_filepath, intrinsics);
    intrinsics_mutex_.unlock();
  }

  // Early exit if no EXIF calibration exists and we are only processing
  // calibration views.
  const CameraIntrinsicsPrior& intrinsics =
      FindOrDie(intrinsics_, image_filepath);
  if (intrinsics.focal_length.is_set) {
    LOG(INFO) << "Image " << image_filepath
              << " contained an EXIF focal length: "
              << intrinsics.focal_length.value;
  } else if (!options_.only_calibrated_views) {
    LOG(INFO) << "Image " << image_filepath
              << " did not contain an EXIF focal length.";
  } else {
    LOG(INFO) << "Image " << image_filepath
              << " did not contain an EXIF focal length. Skipping this image.";
    return;
  }

  // Extract Features.
  std::vector<Keypoint> keypoints;
  std::vector<Eigen::VectorXf> descriptors;
  ExtractFeatures(options_, image_filepath, &keypoints, &descriptors);

  // Add the relevant image and feature data to the feature matcher. This allows
  // the feature matcher to control fine-grained things like multi-threading and
  // caching. For instance, the matcher may choose to write the descriptors to
  // disk and read them back as needed.
  std::string image_filename;
  CHECK(GetFilenameFromFilepath(image_filepath, true, &image_filename));
  matcher_mutex_.lock();
  matcher_->AddImage(image_filename, keypoints, descriptors, intrinsics);
  matcher_mutex_.unlock();
}

}  // namespace theia
