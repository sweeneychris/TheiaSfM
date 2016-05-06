// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <time.h>
#include <theia/theia.h>
#include <chrono>
#include <string>
#include <vector>

#include "applications/command_line_helpers.h"

DEFINE_string(input_features, "",
              "Filepath of the features that you want to perform matching "
              "against. This should be a wildcard that encapsulates all "
              "features files");
DEFINE_string(matching_strategy, "CASCADE_HASHING",
              "Strategy used to match features. Must be BRUTE_FORCE "
              " or CASCADE_HASHING");
DEFINE_bool(geometrically_verifiy_matches, true,
            "Set to true to perform geometric verification on the matches.");
DEFINE_string(calibration_file, "",
              "Optional calibration file specifying camera calibration for the "
              "images. This file may be created from "
              "create_calibration_file_from_exif.cc");

DEFINE_bool(match_out_of_core, true,
            "Perform matching out of core by saving features to disk and "
            "reading them as needed. Set to false to perform matching all in "
            "memory.");
DEFINE_int32(matching_max_num_images_in_cache, 64,
             "Maximum number of images to store in the LRU cache during "
             "feature matching. The higher this number is the more memory is "
             "consumed during matching.");
DEFINE_double(lowes_ratio, 0.75, "Lowes ratio used for feature matching.");
DEFINE_double(
    max_sampson_error_for_verified_match, 4.0,
    "Maximum sampson error for a match to be considered geometrically valid.");
DEFINE_int32(min_num_inliers_for_valid_match, 30,
             "Minimum number of geometrically verified inliers that a pair on "
             "images must have in order to be considered a valid two-view "
             "match.");
DEFINE_bool(bundle_adjust_two_view_geometry, true,
            "Set to false to turn off 2-view BA.");
DEFINE_bool(keep_only_symmetric_matches, true,
            "Performs two-way matching and keeps symmetric matches.");
DEFINE_int32(num_threads, 1,
             "Number of threads to use for feature extraction and matching.");
DEFINE_string(output_matches_file, "",
              "Filepath that the matches file should be written to.");

void SetMatchingOptions(theia::FeatureMatcherOptions* matching_options) {
  matching_options->match_out_of_core = FLAGS_match_out_of_core;
  theia::GetDirectoryFromFilepath(
      FLAGS_input_features,
      &matching_options->keypoints_and_descriptors_output_dir);
  matching_options->cache_capacity = FLAGS_matching_max_num_images_in_cache;
  matching_options->lowes_ratio = FLAGS_lowes_ratio;
  matching_options->keep_only_symmetric_matches =
      FLAGS_keep_only_symmetric_matches;
  matching_options->num_threads = FLAGS_num_threads;
  matching_options->min_num_feature_matches =
      FLAGS_min_num_inliers_for_valid_match;
  matching_options->perform_geometric_verification =
      FLAGS_geometrically_verifiy_matches;
  matching_options->geometric_verification_options.estimate_twoview_info_options
      .max_sampson_error_pixels = FLAGS_max_sampson_error_for_verified_match;
  matching_options->geometric_verification_options.bundle_adjustment =
      FLAGS_bundle_adjust_two_view_geometry;
}

// Gets the image filenames and filepaths from the features filepaths.
void GetImageFilesAndFilenames(
    const std::vector<std::string>& features_filepaths,
    std::vector<std::string>* image_filenames) {
  image_filenames->resize(features_filepaths.size());
  for (int i = 0; i < features_filepaths.size(); i++) {
    const std::size_t features_pos = features_filepaths[i].find(".features");
    const std::string image_filepath =
        features_filepaths[i].substr(0, features_pos);
    CHECK(theia::GetFilenameFromFilepath(image_filepath,
                                         true,
                                         &image_filenames->at(i)));
  }
}

// Read the camera intrinsics from a file.
void ReadIntrinsicsFromCalibrationFile(
    const std::vector<std::string>& image_filenames,
    std::vector<theia::CameraIntrinsicsPrior>* intrinsics) {
  intrinsics->resize(image_filenames.size());

  std::unordered_map<std::string, theia::CameraIntrinsicsPrior> intrinsics_map;
  if (!theia::ReadCalibration(FLAGS_calibration_file, &intrinsics_map)) {
    return;
  }

  for (int i = 0; i < image_filenames.size(); i++) {
    if (theia::ContainsKey(intrinsics_map, image_filenames[i])) {
      intrinsics->at(i) = theia::FindOrDie(intrinsics_map, image_filenames[i]);
    }
  }
}

int main(int argc, char *argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  // Get the filepaths of the features files.
  std::vector<std::string> features_filepaths;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_input_features,
                                        &features_filepaths));
  CHECK_GT(features_filepaths.size(), 0)
      << "No images found in: " << FLAGS_input_features;

  // Get the image files and image filenames.
  std::vector<std::string> image_filenames;
  GetImageFilesAndFilenames(features_filepaths, &image_filenames);

  // Create the feature matcher.
  theia::FeatureMatcherOptions matching_options;
  SetMatchingOptions(&matching_options);
  const theia::MatchingStrategy matching_strategy =
      StringToMatchingStrategyType(FLAGS_matching_strategy);
  std::unique_ptr<theia::FeatureMatcher> matcher =
      CreateFeatureMatcher(matching_strategy, matching_options);

  // Optionally read the intrinsics from a calibration file.
  std::vector<theia::CameraIntrinsicsPrior> intrinsics;
  ReadIntrinsicsFromCalibrationFile(image_filenames, &intrinsics);

  // Add all the features to the matcher.
  for (int i = 0; i < features_filepaths.size(); i++) {
    matcher->AddImage(image_filenames[i], intrinsics[i]);
  }

  // Match the images with optional geometric verification.
  std::vector<theia::ImagePairMatch> matches;
  matcher->MatchImages(&matches);

  // Write the matches out.
  LOG(INFO) << "Writing matches to file: " << FLAGS_output_matches_file;
  CHECK(theia::WriteMatchesAndGeometry(FLAGS_output_matches_file,
                                       image_filenames, intrinsics, matches))
      << "Could not write the matches to " << FLAGS_output_matches_file;
}
