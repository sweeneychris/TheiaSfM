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

#include <chrono>  // NOLINT
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <string>
#include <theia/theia.h>
#include <time.h>
#include <vector>

#include "applications/command_line_helpers.h"
#include "applications/print_reconstruction_statistics.h"

// Number of iterations to run the calibration.
DEFINE_int32(num_calibration_iterations,
             2,
             "Number of iterations to run the calibration. Typically 2-3 "
             "iterations is enough to get a stable result.");
DEFINE_string(camera_model,
              "PINHOLE",
              "The type of camera model to use for calibration. The most "
              "common camera model types are PINHOLE and FISHEYE but a full "
              "list may be found at "
              "//theia/sfm/camera/camera_intrinsic_model_type.h");

// Input/output files.
DEFINE_string(images, "", "Wildcard of images to reconstruct.");

// Multithreading.
DEFINE_int32(num_threads,
             1,
             "Number of threads to use for feature extraction and matching.");

// Feature and matching options.
DEFINE_string(
    descriptor,
    "SIFT",
    "Type of feature descriptor to use. Must be one of the following: "
    "SIFT");
DEFINE_string(feature_density,
              "NORMAL",
              "Set to SPARSE, NORMAL, or DENSE to extract fewer or more "
              "features from each image.");
DEFINE_string(matching_strategy,
              "CASCADE_HASHING",
              "Strategy used to match features. Must be BRUTE_FORCE "
              " or CASCADE_HASHING");
DEFINE_bool(match_out_of_core,
            true,
            "Perform matching out of core by saving features to disk and "
            "reading them as needed. Set to false to perform matching all in "
            "memory.");
DEFINE_string(matching_working_directory,
              "",
              "Directory used during matching to store features for "
              "out-of-core matching.");
DEFINE_int32(matching_max_num_images_in_cache,
             128,
             "Maximum number of images to store in the LRU cache during "
             "feature matching. The higher this number is the more memory is "
             "consumed during matching.");
DEFINE_double(lowes_ratio, 0.8, "Lowes ratio used for feature matching.");
DEFINE_double(max_sampson_error_for_verified_match,
              6.0,
              "Maximum sampson error for a match to be considered "
              "geometrically valid. This threshold is relative to an image "
              "with a width of 1024 pixels and will be appropriately scaled "
              "for images with different resolutions.");
DEFINE_int32(min_num_inliers_for_valid_match,
             30,
             "Minimum number of geometrically verified inliers that a pair on "
             "images must have in order to be considered a valid two-view "
             "match.");
DEFINE_bool(bundle_adjust_two_view_geometry,
            true,
            "Set to false to turn off 2-view BA.");
DEFINE_bool(keep_only_symmetric_matches,
            true,
            "Performs two-way matching and keeps symmetric matches.");

// Reconstruction building options.
DEFINE_string(intrinsics_to_optimize,
              "NONE",
              "Set to control which intrinsics parameters are optimized during "
              "bundle adjustment.");
DEFINE_double(max_reprojection_error_pixels,
              4.0,
              "Maximum reprojection error for a correspondence to be "
              "considered an inlier after bundle adjustment.");

// Incremental SfM options.
DEFINE_double(absolute_pose_reprojection_error_threshold,
              4.0,
              "The inlier threshold for absolute pose estimation. This "
              "threshold is relative to an image with a width of 1024 pixels "
              "and will be appropriately scaled based on the input image "
              "resolutions.");
DEFINE_int32(min_num_absolute_pose_inliers,
             30,
             "Minimum number of inliers in order for absolute pose estimation "
             "to be considered successful.");
DEFINE_double(full_bundle_adjustment_growth_percent,
              5.0,
              "Full BA is only triggered for incremental SfM when the "
              "reconstruction has growth by this percent since the last time "
              "full BA was used.");
DEFINE_int32(partial_bundle_adjustment_num_views,
             20,
             "When full BA is not being run, partial BA is executed on a "
             "constant number of views specified by this parameter.");

// Triangulation options.
DEFINE_double(min_triangulation_angle_degrees,
              4.0,
              "Minimum angle between views for triangulation.");
DEFINE_double(
    triangulation_reprojection_error_pixels,
    15.0,
    "Max allowable reprojection error on initial triangulation of points.");
DEFINE_bool(bundle_adjust_tracks,
            true,
            "Set to true to optimize tracks immediately upon estimation.");

// Bundle adjustment parameters.
DEFINE_string(bundle_adjustment_robust_loss_function,
              "CAUCHY",
              "By setting this to an option other than NONE, a robust loss "
              "function will be used during bundle adjustment which can "
              "improve robustness to outliers. Options are NONE, HUBER, "
              "SOFTLONE, CAUCHY, ARCTAN, and TUKEY.");
DEFINE_double(bundle_adjustment_robust_loss_width,
              5.0,
              "If the BA loss function is not NONE, then this value controls "
              "where the robust loss begins with respect to reprojection error "
              "in pixels.");

using theia::Reconstruction;
using theia::ReconstructionBuilder;
using theia::ReconstructionBuilderOptions;

// Sets the feature extraction, matching, and reconstruction options based on
// the command line flags. There are many more options beside just these located
// in //theia/vision/sfm/reconstruction_builder.h
ReconstructionBuilderOptions SetReconstructionBuilderOptions() {
  ReconstructionBuilderOptions options;
  options.num_threads = FLAGS_num_threads;

  options.descriptor_type = StringToDescriptorExtractorType(FLAGS_descriptor);
  options.feature_density = StringToFeatureDensity(FLAGS_feature_density);
  options.match_out_of_core = FLAGS_match_out_of_core;
  options.features_and_matches_database_directory =
      FLAGS_matching_working_directory;
  options.cache_capacity = FLAGS_matching_max_num_images_in_cache;
  options.matching_strategy =
      StringToMatchingStrategyType(FLAGS_matching_strategy);
  options.matching_options.lowes_ratio = FLAGS_lowes_ratio;
  options.matching_options.keep_only_symmetric_matches =
      FLAGS_keep_only_symmetric_matches;
  options.min_num_inlier_matches = FLAGS_min_num_inliers_for_valid_match;
  options.matching_options.perform_geometric_verification = true;
  options.matching_options.geometric_verification_options
      .estimate_twoview_info_options.max_sampson_error_pixels =
      FLAGS_max_sampson_error_for_verified_match;
  options.matching_options.geometric_verification_options.bundle_adjustment =
      FLAGS_bundle_adjust_two_view_geometry;
  options.matching_options.geometric_verification_options
      .triangulation_max_reprojection_error =
      FLAGS_triangulation_reprojection_error_pixels;
  options.matching_options.geometric_verification_options
      .min_triangulation_angle_degrees = FLAGS_min_triangulation_angle_degrees;
  options.matching_options.geometric_verification_options
      .final_max_reprojection_error = FLAGS_max_reprojection_error_pixels;

  // Reconstruction Estimator Options.
  theia::ReconstructionEstimatorOptions& reconstruction_estimator_options =
      options.reconstruction_estimator_options;
  reconstruction_estimator_options.min_num_two_view_inliers =
      FLAGS_min_num_inliers_for_valid_match;
  reconstruction_estimator_options.num_threads = FLAGS_num_threads;
  reconstruction_estimator_options.intrinsics_to_optimize =
      StringToOptimizeIntrinsicsType(FLAGS_intrinsics_to_optimize);
  options.reconstruct_largest_connected_component = true;
  options.only_calibrated_views = false;
  reconstruction_estimator_options.max_reprojection_error_in_pixels =
      FLAGS_max_reprojection_error_pixels;

  // Which type of SfM pipeline to use (e.g., incremental, global, etc.);
  reconstruction_estimator_options.reconstruction_estimator_type =
      theia::ReconstructionEstimatorType::INCREMENTAL;

  // Incremental SfM Options.
  reconstruction_estimator_options.absolute_pose_reprojection_error_threshold =
      FLAGS_absolute_pose_reprojection_error_threshold;
  reconstruction_estimator_options.min_num_absolute_pose_inliers =
      FLAGS_min_num_absolute_pose_inliers;
  reconstruction_estimator_options.full_bundle_adjustment_growth_percent =
      FLAGS_full_bundle_adjustment_growth_percent;
  reconstruction_estimator_options.partial_bundle_adjustment_num_views =
      FLAGS_partial_bundle_adjustment_num_views;

  // Triangulation options.
  reconstruction_estimator_options.min_triangulation_angle_degrees =
      FLAGS_min_triangulation_angle_degrees;
  reconstruction_estimator_options
      .triangulation_max_reprojection_error_in_pixels =
      FLAGS_triangulation_reprojection_error_pixels;
  reconstruction_estimator_options.bundle_adjust_tracks =
      FLAGS_bundle_adjust_tracks;

  // Bundle adjustment options.
  reconstruction_estimator_options.bundle_adjustment_loss_function_type =
      StringToLossFunction(FLAGS_bundle_adjustment_robust_loss_function);
  reconstruction_estimator_options.bundle_adjustment_robust_loss_width =
      FLAGS_bundle_adjustment_robust_loss_width;
  return options;
}

void AddImagesToReconstructionBuilder(
    ReconstructionBuilder* reconstruction_builder,
    const theia::CameraIntrinsicsPrior& prior) {
  std::vector<std::string> image_files;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_images, &image_files))
      << "Could not find images that matched the filepath: " << FLAGS_images
      << ". NOTE that the ~ filepath is not supported.";

  CHECK_GT(image_files.size(), 0) << "No images found in: " << FLAGS_images;

  // Add images with possible calibration. When the intrinsics group id is
  // invalid, the reconstruction builder will assume that the view does not
  // share its intrinsics with any other views.
  theia::CameraIntrinsicsGroupId intrinsics_group_id = 0;

  for (const std::string& image_file : image_files) {
    std::string image_filename;
    CHECK(theia::GetFilenameFromFilepath(image_file, true, &image_filename));
    CHECK(reconstruction_builder->AddImageWithCameraIntrinsicsPrior(
        image_file, prior, intrinsics_group_id));
  }

  // Extract and match features.
  CHECK(reconstruction_builder->ExtractAndMatchFeatures());
}

int main(int argc, char* argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  LOG(INFO) << "This calibration technique simply runs incremental SfM several "
               "times while refining the intrinsic parameters. After each "
               "iteration, the intrinsic parameters output are then used as "
               "the initialization point for the next round of SfM. Each "
               "iteration recomputes matching and the reconstruction from "
               "scratch to obtain the best possible results. It is recommended "
               "to use at least 5 images and up to 100 images with a wide "
               "range of motion for calibration.";

  const ReconstructionBuilderOptions options =
      SetReconstructionBuilderOptions();

  // You can set the camera intrinsics prior explicitly here, or let Theia try
  // to initialize the values for you. Be sure to set the "is_set" field to
  // "true" for any intrinsics values that you would like to explicitly set. If
  // the "is_set" field is not set to true, the value(s) for that intrinsic
  // parameter will not be used for calibration.
  theia::CameraIntrinsicsPrior prior;
  prior.camera_intrinsics_model_type = FLAGS_camera_model;

  // prior.image_width = 1920;
  // prior.image_height = 1080;
  // prior.camera_intrinsics_model_type = "PINHOLE";

  // prior.focal_length.is_set = true;
  // prior.focal_length.value[0] = 1784.759;

  // prior.principal_point.is_set = true;
  // prior.principal_point.value[0] = 960.0;
  // prior.principal_point.value[1] = 540.0;

  // prior.aspect_ratio.is_set = true;
  // prior.aspect_ratio.value[0] = 1.0;

  // prior.skew.is_set = true;
  // prior.skew.value[0] = 0.0;

  // prior.radial_distortion.is_set = true;
  // prior.radial_distortion.value[0] = 0.840712;
  // prior.radial_distortion.value[1] = 0.0;
  // prior.radial_distortion.value[2] = 0.0;
  // prior.radial_distortion.value[3] = 0.0;

  // prior.tangential_distortion.is_set = true;
  // prior.tangential_distortion.value[0] = 0.840712;
  // prior.tangential_distortion.value[1] = 0.0;

  std::vector<Reconstruction*> reconstructions;
  for (int i = 0; i < FLAGS_num_calibration_iterations; i++) {
    reconstructions.clear();

    ReconstructionBuilder reconstruction_builder(options);
    AddImagesToReconstructionBuilder(&reconstruction_builder, prior);

    CHECK(reconstruction_builder.BuildReconstruction(&reconstructions))
        << "Could not create a reconstruction.";

    // Print the final reconstruction statistics so the user may assess the
    // quality of the calibration.
    PrintReprojectionErrors(*reconstructions[0]);
    PrintTrackLengthHistogram(*reconstructions[0]);

    // Print the output camera parameters.
    const auto& view_ids = reconstructions[0]->ViewIds();
    const theia::View* view = reconstructions[0]->View(view_ids[0]);
    CHECK(view != nullptr)
        << "Calibration failed! Could not reconstruct the scene at iteration "
        << i + 1;
    view->Camera().PrintCameraIntrinsics();

    // Use the camera intrinsics from the optimized reconstruction as the
    // initialization for the next iteration.
    prior = view->Camera().CameraIntrinsicsPriorFromIntrinsics();
  }
}
