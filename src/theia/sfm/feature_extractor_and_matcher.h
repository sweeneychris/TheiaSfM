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

#ifndef THEIA_SFM_FEATURE_EXTRACTOR_AND_MATCHER_H_
#define THEIA_SFM_FEATURE_EXTRACTOR_AND_MATCHER_H_

#include <string>
#include <thread>  // NOLINT
#include <vector>

#include "theia/image/descriptor/create_descriptor_extractor.h"
#include "theia/matching/create_feature_matcher.h"
#include "theia/matching/feature_matcher_options.h"
#include "theia/matching/image_pair_match.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/exif_reader.h"
#include "theia/sfm/two_view_match_geometric_verification.h"

namespace theia {

class FeatureExtractorAndMatcher {
 public:
  struct Options {
    // Number of threads for multithreading.
    int num_threads = 1;

    // If true, only images that contain EXIF focal length values have features
    // extracted and matched, and images that do not contain EXIF focal length
    // are not considered for the feature extraction and matching.
    bool only_calibrated_views = false;

    // The type of feature to use for feature extraction.
    DescriptorExtractorType descriptor_extractor_type =
        DescriptorExtractorType::SIFT;

    // The density of features to extract. DENSE means more features are
    // extracted per image and SPARSE means fewer features per image are
    // extracted.
    FeatureDensity feature_density = FeatureDensity::NORMAL;

    // The features returned will be no larger than this size.
    int max_num_features = 16384;

    // Minimum number of inliers to consider the matches a good match.
    int min_num_inlier_matches = 30;

    // Matching strategy to use for establishing feature correspondences.
    MatchingStrategy matching_strategy;

    // Matching options for determining which feature matches are good matches.
    FeatureMatcherOptions feature_matcher_options;
  };

  explicit FeatureExtractorAndMatcher(const Options& options);

  // Add an image to the image matcher queue.
  bool AddImage(const std::string& image_filepath);

  // Add an image with known camera intrinsics to the image matcher queue.
  bool AddImage(const std::string& image_filepath,
                const CameraIntrinsicsPrior& intrinsics);

  // Sort and fill the vector of masks in order to make a correspondance between
  // the idexes of 'image_filepaths_' and 'mask_filepaths_'.
  bool SetMasksForFeaturesExtraction(
      const std::vector<std::string>& mask_filepaths);

  // Performs feature matching between all images provided by the image
  // filepaths. Features are extracted and matched between the images according
  // to the options passed in. Only matches that have passed geometric
  // verification are kept. EXIF data is parsed to determine the camera
  // intrinsics if available.
  void ExtractAndMatchFeatures(
      std::vector<CameraIntrinsicsPrior>* intrinsics,
      std::vector<ImagePairMatch>* matches);

 private:
  // Processes a single image by extracting EXIF information, extracting
  // features and descriptors, and adding the image to the matcher.
  void ProcessImage(const int i);

  const Options options_;

  // Local copies of the images to be matches, masks for use and any priors on
  // the camera intrinsics.
  std::vector<std::string> image_filepaths_;
  std::unordered_map<std::string, CameraIntrinsicsPrior> intrinsics_;

  // Local copies of the sorted masks (indexed like 'image_filepaths_').
  std::vector<std::string> mask_filepaths_;

  // Exif reader for loading exif information. This object is created once so
  // that the EXIF focal length database does not have to be loaded multiple
  // times.
  ExifReader exif_reader_;

  // Feature matcher and mutex for thread-safe access.
  std::unique_ptr<FeatureMatcher> matcher_;
  std::mutex intrinsics_mutex_, matcher_mutex_;
};

}  // namespace theia

#endif  // THEIA_SFM_FEATURE_EXTRACTOR_AND_MATCHER_H_
