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

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <time.h>
#include <theia/theia.h>
#include <chrono>  // NOLINT
#include <string>
#include <vector>

// Input/output files.
DEFINE_string(images, "", "Wildcard of images to reconstruct.");
DEFINE_string(matches_file, "", "Filename of the matches file.");
DEFINE_string(calibration_file, "",
              "Calibration file containing image calibration data.");
DEFINE_string(
    output_matches_file, "",
    "File to write the two-view matches to. This file can be used in "
    "future iterations as input to the reconstruction builder. Leave empty if "
    "you do not want to output matches.");
DEFINE_string(
    output_reconstruction, "",
    "Filename to write reconstruction to. The filename will be appended with "
    "the reconstruction number if multiple reconstructions are created.");

// Multithreading.
DEFINE_int32(num_threads, 1,
             "Number of threads to use for feature extraction and matching.");

// Feature and matching options.
DEFINE_string(
    descriptor, "SIFT",
    "Type of feature descriptor to use. Must be one of the following: "
    "SIFT, BRIEF, BRISK, FREAK");
DEFINE_string(matching_strategy, "BRUTE_FORCE",
              "Strategy used to match features. Must be BRUTE_FORCE "
              " or CASCADE_HASHING");
DEFINE_double(lowes_ratio, 0.8, "Lowes ratio used for feature matching.");
DEFINE_double(
    max_sampson_error_for_verified_match, 4.0,
    "Maximum sampson error for a match to be considered geometrically valid.");
DEFINE_int32(min_num_inliers_for_valid_match, 30,
             "Minimum number of geometrically verified inliers that a pair on "
             "images must have in order to be considered a valid two-view "
             "match.");
DEFINE_bool(bundle_adjust_two_view_geometry, true,
            "Set to false to turn off 2-view BA.");

// Reconstruction building options.
DEFINE_string(reconstruction_estimator, "NONLINEAR",
              "Type of global SfM reconstruction estimation to use.");
DEFINE_bool(reconstruct_largest_connected_component, false,
            "If set to true, only the single largest connected component is "
            "reconstructed. Otherwise, as many models as possible are "
            "estimated.");
DEFINE_bool(only_calibrated_views, false,
            "Set to true to only reconstruct the views where calibration is "
            "provided or can be extracted from EXIF");
DEFINE_int32(max_track_length, 20, "Maximum length of a track.");
DEFINE_bool(initialize_focal_lengths_from_median_estimate, false,
            "If false, initialize focal lengths from a median viewing angle if "
            "EXIF is not available. If true, then use the median focal length "
            "from fundamental matrix decompositions.");
DEFINE_bool(constant_camera_intrinsics, false,
            "Set to true to keep camera intrinsic parameters constant during "
            "bundle adjustment.");
DEFINE_double(max_reprojection_error_pixels, 4.0,
              "Maximum reprojection error for a correspondence to be "
              "considered an inlier. This is used for absolute pose estimation "
              "as well as outlier filtering after bundle adjustment.");
DEFINE_int32(num_retriangulation_iterations, 1,
             "Number of times to retriangulate any unestimated tracks. Bundle "
             "adjustment is performed after retriangulation.");

// View graph filtering options.
DEFINE_bool(refine_relative_translations_after_rotation_estimation, true,
            "Refine the relative translation estimation after computing the "
            "absolute rotations. This can help improve the accuracy of the "
            "position estimation.");
DEFINE_double(post_rotation_filtering_degrees, 5.0,
              "Max degrees difference in relative rotation and rotation "
              "estimates for rotation filtering.");

// Position estimation options.
DEFINE_int32(
    position_estimation_min_num_tracks_per_view, 0,
    "Minimum number of point to camera constraints for position estimation.");
DEFINE_int32(position_estimation_max_num_reweighted_iterations, 100,
             "Maximum number of reweighted iterations to perform for position "
             "estimation. Set to zero if only a single robust optimization is "
             "desired.");

// Triangulation options.
DEFINE_double(min_triangulation_angle_degrees, 4.0,
              "Minimum angle between views for triangulation.");
DEFINE_double(
    triangulation_reprojection_error_pixels, 15.0,
    "Max allowable reprojection error on initial triangulation of points.");
DEFINE_bool(bundle_adjust_tracks, true,
            "Set to true to optimize tracks immediately upon estimation.");

using theia::DescriptorExtractorType;
using theia::MatchingStrategy;
using theia::Reconstruction;
using theia::ReconstructionBuilder;
using theia::ReconstructionBuilderOptions;
using theia::ReconstructionEstimatorType;

DescriptorExtractorType GetDescriptorExtractorType(
    const std::string& descriptor) {
  if (descriptor == "SIFT") {
    return DescriptorExtractorType::SIFT;
  } else if (descriptor == "BRIEF") {
    return DescriptorExtractorType::BRIEF;
  } else if (descriptor == "BRISK") {
    return DescriptorExtractorType::BRISK;
  } else if (descriptor == "FREAK") {
    return DescriptorExtractorType::FREAK;
  } else {
    LOG(FATAL) << "Invalid DescriptorExtractor specified. Using SIFT instead.";
    return DescriptorExtractorType::SIFT;
  }
}

MatchingStrategy GetMatchingStrategyType(
    const std::string& matching_strategy) {
  if (matching_strategy == "BRUTE_FORCE") {
    return MatchingStrategy::BRUTE_FORCE;
  } else if (matching_strategy == "CASCADE_HASHING") {
    return MatchingStrategy::CASCADE_HASHING;
  } else {
    LOG(FATAL)
        << "Invalid matching strategy specified. Using BRUTE_FORCE instead.";
    return MatchingStrategy::BRUTE_FORCE;
  }
}

ReconstructionEstimatorType GetReconstructionEstimatorType(
    const std::string& reconstruction_estimator) {
  if (reconstruction_estimator == "NONLINEAR") {
    return ReconstructionEstimatorType::NONLINEAR;
  } else {
    LOG(FATAL)
        << "Invalid reconstruction estimator type. Using NONLINEAR instead.";
    return ReconstructionEstimatorType::NONLINEAR;
  }
}

// Sets the feature extraction, matching, and reconstruction options based on
// the command line flags. There are many more options beside just these located
// in //theia/vision/sfm/reconstruction_builder.h
ReconstructionBuilderOptions SetReconstructionBuilderOptions() {
  ReconstructionBuilderOptions options;
  options.num_threads = FLAGS_num_threads;
  options.reconstruction_estimator_options.num_threads = FLAGS_num_threads;

  options.descriptor_type = GetDescriptorExtractorType(FLAGS_descriptor);
  options.matching_strategy = GetMatchingStrategyType(FLAGS_matching_strategy);
  options.matching_options.lowes_ratio = FLAGS_lowes_ratio;
  options.min_num_inlier_matches = FLAGS_min_num_inliers_for_valid_match;
  options.reconstruction_estimator_options.min_num_two_view_inliers =
      FLAGS_min_num_inliers_for_valid_match;
  options.geometric_verification_options.estimate_twoview_info_options
      .max_sampson_error_pixels = FLAGS_max_sampson_error_for_verified_match;
  options.geometric_verification_options.bundle_adjustment =
      FLAGS_bundle_adjust_two_view_geometry;

  options.max_track_length = FLAGS_max_track_length;
  options.reconstruction_estimator_options.constant_camera_intrinsics =
      FLAGS_constant_camera_intrinsics;

  options.reconstruction_estimator_options.reconstruction_estimator_type =
      GetReconstructionEstimatorType(FLAGS_reconstruction_estimator);
  options.reconstruct_largest_connected_component =
      FLAGS_reconstruct_largest_connected_component;
  options.only_calibrated_views = FLAGS_only_calibrated_views;
  options.reconstruction_estimator_options
      .initialize_focal_lengths_from_median_estimate =
      FLAGS_initialize_focal_lengths_from_median_estimate;
  options.reconstruction_estimator_options.max_reprojection_error_in_pixels =
      FLAGS_max_reprojection_error_pixels;
  options.reconstruction_estimator_options.num_retriangulation_iterations =
      FLAGS_num_retriangulation_iterations;

  options.reconstruction_estimator_options
      .refine_relative_translations_after_rotation_estimation =
      FLAGS_refine_relative_translations_after_rotation_estimation;
  options.reconstruction_estimator_options
      .rotation_filtering_max_difference_degrees =
      FLAGS_post_rotation_filtering_degrees;
  options.reconstruction_estimator_options
      .position_estimation_min_num_tracks_per_view =
      FLAGS_position_estimation_min_num_tracks_per_view;
  options.reconstruction_estimator_options
      .position_estimation_max_reweighted_iterations =
      FLAGS_position_estimation_max_num_reweighted_iterations;
  options.output_matches_file = FLAGS_output_matches_file;

  options.reconstruction_estimator_options.min_triangulation_angle_degrees =
      FLAGS_min_triangulation_angle_degrees;
  options.reconstruction_estimator_options
      .triangulation_max_reprojection_error_in_pixels =
      FLAGS_triangulation_reprojection_error_pixels;
  options.reconstruction_estimator_options.bundle_adjust_tracks =
      FLAGS_bundle_adjust_tracks;
  return options;
}

void AddMatchesToReconstructionBuilder(
    ReconstructionBuilder* reconstruction_builder) {
  // Load matches from file.
  std::vector<std::string> image_files;
  std::vector<theia::CameraIntrinsicsPrior> camera_intrinsics_prior;
  std::vector<theia::ImagePairMatch> image_matches;

  // Read in match file.
  theia::ReadMatchesAndGeometry(FLAGS_matches_file,
                                &image_files,
                                &camera_intrinsics_prior,
                                &image_matches);

  // Add all the views.
  for (int i = 0; i < image_files.size(); i++) {
    reconstruction_builder->AddImageWithCameraIntrinsicsPrior(
        image_files[i], camera_intrinsics_prior[i]);
  }

  // Add the matches.
  for (const auto& match : image_matches) {
    const std::string& image1 = image_files[match.image1_index];
    const std::string& image2 = image_files[match.image2_index];
    CHECK(reconstruction_builder->AddTwoViewMatch(image1, image2, match));
  }
}

void AddImagesToReconstructionBuilder(
    ReconstructionBuilder* reconstruction_builder) {
  std::vector<std::string> image_files;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_images, &image_files))
      << "Could not find images that matched the filepath: " << FLAGS_images
      << ". NOTE that the ~ filepath is not supported.";

  CHECK_GT(image_files.size(), 0) << "No images found in: " << FLAGS_images;

  // Load calibration file if it is provided.
  std::unordered_map<std::string, theia::CameraIntrinsicsPrior>
      camera_intrinsics_prior;
  if (FLAGS_calibration_file.size() != 0) {
    CHECK(theia::ReadCalibration(FLAGS_calibration_file,
                                 &camera_intrinsics_prior))
        << "Could not read calibration file.";
  }

  // Add images with possible calibration.
  for (const std::string& image_file : image_files) {
    std::string image_filename;
    CHECK(theia::GetFilenameFromFilepath(image_file, true, &image_filename));

    const theia::CameraIntrinsicsPrior* image_camera_intrinsics_prior =
        FindOrNull(camera_intrinsics_prior, image_filename);
    if (image_camera_intrinsics_prior != nullptr) {
      CHECK(reconstruction_builder->AddImageWithCameraIntrinsicsPrior(
          image_file, *image_camera_intrinsics_prior));
    } else {
      CHECK(reconstruction_builder->AddImage(image_file));
    }
  }

  // Extract and match features.
  CHECK(reconstruction_builder->ExtractAndMatchFeatures());
}

int main(int argc, char *argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  CHECK_GT(FLAGS_output_reconstruction.size(), 0)
      << "Must specify a filepath to output the reconstruction.";

  const ReconstructionBuilderOptions options =
      SetReconstructionBuilderOptions();

  ReconstructionBuilder reconstruction_builder(options);
  // If matches are provided, load matches otherwise load images.
  if (FLAGS_matches_file.size() != 0) {
    AddMatchesToReconstructionBuilder(&reconstruction_builder);
  } else if (FLAGS_images.size() != 0) {
    AddImagesToReconstructionBuilder(&reconstruction_builder);
  } else {
    LOG(FATAL)
        << "You must specifiy either images to reconstruct or a match file.";
  }

  std::vector<Reconstruction*> reconstructions;
  CHECK(reconstruction_builder.BuildReconstruction(&reconstructions))
      << "Could not create a reconstruction.";

  for (int i = 0; i < reconstructions.size(); i++) {
    const std::string output_file =
        theia::StringPrintf("%s-%d", FLAGS_output_reconstruction.c_str(), i);
    LOG(INFO) << "Writing reconstruction " << i << " to " << output_file;
    CHECK(theia::WriteReconstruction(*reconstructions[i], output_file))
        << "Could not write reconstruction to file.";
  }
}
