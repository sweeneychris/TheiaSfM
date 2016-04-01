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

#include "applications/command_line_helpers.h"

// Input/output files.
DEFINE_string(1dsfm_dataset_directory, "",
              "Dataset where the 1dSFM dataset is located. Do not include a "
              "trailing slash.");
DEFINE_string(
    output_reconstruction, "",
    "Filename to write reconstruction to. The filename will be appended with "
    "the reconstruction number if multiple reconstructions are created.");

// Multithreading.
DEFINE_int32(num_threads, 1,
             "Number of threads to use for feature extraction and matching.");

// Reconstruction building options.
DEFINE_string(reconstruction_estimator, "GLOBAL",
              "Type of SfM reconstruction estimation to use.");
DEFINE_int32(min_num_inliers_for_valid_match, 30,
             "Minimum number of geometrically verified inliers that a pair on "
             "images must have in order to be considered a valid two-view "
             "match.");
DEFINE_bool(reconstruct_largest_connected_component, false,
            "If set to true, only the single largest connected component is "
            "reconstructed. Otherwise, as many models as possible are "
            "estimated.");
DEFINE_bool(only_calibrated_views, false,
            "Set to true to only reconstruct the views where calibration is "
            "provided or can be extracted from EXIF");
DEFINE_int32(max_track_length, 20, "Maximum length of a track.");
DEFINE_string(intrinsics_to_optimize,
              "NONE",
              "Set to control which intrinsics parameters are optimized during "
              "bundle adjustment.");
DEFINE_double(max_reprojection_error_pixels, 4.0,
              "Maximum reprojection error for a correspondence to be "
              "considered an inlier after bundle adjustment.");

// Global SfM options.
DEFINE_string(global_rotation_estimator, "ROBUST_L1L2",
              "Type of global rotation estimation to use for global SfM.");
DEFINE_string(global_position_estimator, "NONLINEAR",
              "Type of global position estimation to use for global SfM.");
DEFINE_bool(refine_relative_translations_after_rotation_estimation, true,
            "Refine the relative translation estimation after computing the "
            "absolute rotations. This can help improve the accuracy of the "
            "position estimation.");
DEFINE_double(post_rotation_filtering_degrees, 5.0,
              "Max degrees difference in relative rotation and rotation "
              "estimates for rotation filtering.");
DEFINE_bool(extract_maximal_rigid_subgraph, false,
            "If true, only cameras that are well-conditioned for position "
            "estimation will be used for global position estimation.");
DEFINE_bool(filter_relative_translations_with_1dsfm, true,
            "Filter relative translation estimations with the 1DSfM algorithm "
            "to potentially remove outlier relativep oses for position "
            "estimation.");
DEFINE_int32(num_retriangulation_iterations, 1,
             "Number of times to retriangulate any unestimated tracks. Bundle "
             "adjustment is performed after retriangulation.");

// Nonlinear position estimation options.
DEFINE_int32(
    position_estimation_min_num_tracks_per_view, 0,
    "Minimum number of point to camera constraints for position estimation.");
DEFINE_double(position_estimation_robust_loss_width, 0.1,
              "Robust loss width to use for position estimation.");

// Incremental SfM options.
DEFINE_double(absolute_pose_reprojection_error_threshold, 8.0,
              "The inlier threshold for absolute pose estimation.");
DEFINE_int32(min_num_absolute_pose_inliers, 30,
             "Minimum number of inliers in order for absolute pose estimation "
             "to be considered successful.");
DEFINE_double(full_bundle_adjustment_growth_percent, 5.0,
              "Full BA is only triggered for incremental SfM when the "
              "reconstruction has growth by this percent since the last time "
              "full BA was used.");
DEFINE_int32(partial_bundle_adjustment_num_views, 20,
             "When full BA is not being run, partial BA is executed on a "
             "constant number of views specified by this parameter.");


// Triangulation options.
DEFINE_double(min_triangulation_angle_degrees, 4.0,
              "Minimum angle between views for triangulation.");
DEFINE_double(
    triangulation_reprojection_error_pixels, 15.0,
    "Max allowable reprojection error on initial triangulation of points.");
DEFINE_bool(bundle_adjust_tracks, true,
            "Set to true to optimize tracks immediately upon estimation.");

// Bundle adjustment parameters.
DEFINE_string(bundle_adjustment_robust_loss_function, "NONE",
              "By setting this to an option other than NONE, a robust loss "
              "function will be used during bundle adjustment which can "
              "improve robustness to outliers. Options are NONE, HUBER, "
              "SOFTLONE, CAUCHY, ARCTAN, and TUKEY.");
DEFINE_double(bundle_adjustment_robust_loss_width, 10.0,
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
  options.max_track_length = FLAGS_max_track_length;

  // Reconstruction Estimator Options.
  theia::ReconstructionEstimatorOptions& reconstruction_estimator_options =
      options.reconstruction_estimator_options;
  reconstruction_estimator_options.min_num_two_view_inliers =
      FLAGS_min_num_inliers_for_valid_match;
  reconstruction_estimator_options.num_threads = FLAGS_num_threads;
  reconstruction_estimator_options.intrinsics_to_optimize =
    StringToOptimizeIntrinsicsType(FLAGS_intrinsics_to_optimize);
  options.reconstruct_largest_connected_component =
      FLAGS_reconstruct_largest_connected_component;
  options.only_calibrated_views = FLAGS_only_calibrated_views;
  reconstruction_estimator_options.max_reprojection_error_in_pixels =
      FLAGS_max_reprojection_error_pixels;

  // Which type of SfM pipeline to use (e.g., incremental, global, etc.);
  reconstruction_estimator_options.reconstruction_estimator_type =
      StringToReconstructionEstimatorType(FLAGS_reconstruction_estimator);

  // Global SfM Options.
  reconstruction_estimator_options.global_rotation_estimator_type =
      StringToRotationEstimatorType(FLAGS_global_rotation_estimator);
  reconstruction_estimator_options.global_position_estimator_type =
      StringToPositionEstimatorType(FLAGS_global_position_estimator);
  reconstruction_estimator_options.num_retriangulation_iterations =
      FLAGS_num_retriangulation_iterations;
  reconstruction_estimator_options
      .refine_relative_translations_after_rotation_estimation =
      FLAGS_refine_relative_translations_after_rotation_estimation;
  reconstruction_estimator_options.extract_maximal_rigid_subgraph =
      FLAGS_extract_maximal_rigid_subgraph;
  reconstruction_estimator_options.filter_relative_translations_with_1dsfm =
      FLAGS_filter_relative_translations_with_1dsfm;
  reconstruction_estimator_options
      .rotation_filtering_max_difference_degrees =
      FLAGS_post_rotation_filtering_degrees;
  reconstruction_estimator_options.nonlinear_position_estimator_options
      .min_num_points_per_view =
      FLAGS_position_estimation_min_num_tracks_per_view;

  // Incremental SfM Options.
  reconstruction_estimator_options
      .absolute_pose_reprojection_error_threshold =
      FLAGS_absolute_pose_reprojection_error_threshold;
  reconstruction_estimator_options.min_num_absolute_pose_inliers =
      FLAGS_min_num_absolute_pose_inliers;
  reconstruction_estimator_options
      .full_bundle_adjustment_growth_percent =
      FLAGS_full_bundle_adjustment_growth_percent;
  reconstruction_estimator_options.partial_bundle_adjustment_num_views =
      FLAGS_partial_bundle_adjustment_num_views;

  // Triangulation options (used by all SfM pipelines).
  reconstruction_estimator_options.min_triangulation_angle_degrees =
      FLAGS_min_triangulation_angle_degrees;
  reconstruction_estimator_options
      .triangulation_max_reprojection_error_in_pixels =
      FLAGS_triangulation_reprojection_error_pixels;
  reconstruction_estimator_options.bundle_adjust_tracks =
      FLAGS_bundle_adjust_tracks;

  // Bundle adjustment options (used by all SfM pipelines).
  reconstruction_estimator_options.bundle_adjustment_loss_function_type =
      StringToLossFunction(FLAGS_bundle_adjustment_robust_loss_function);
  reconstruction_estimator_options.bundle_adjustment_robust_loss_width =
      FLAGS_bundle_adjustment_robust_loss_width;
  return options;
}

void InitializeFrom1DSFM(ReconstructionBuilder* reconstruction_builder) {
  std::unique_ptr<Reconstruction> reconstruction(new Reconstruction);
  std::unique_ptr<theia::ViewGraph> view_graph(new theia::ViewGraph);
  CHECK(Read1DSFM(FLAGS_1dsfm_dataset_directory,
                  reconstruction.get(),
                  view_graph.get()))
      << "Could not read 1dsfm dataset from " << FLAGS_1dsfm_dataset_directory;
  LOG(INFO) << "Initializing reconstruction builder from 1dsfm.";
  reconstruction_builder->InitializeReconstructionAndViewGraph(
      reconstruction.release(), view_graph.release());
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
  if (FLAGS_1dsfm_dataset_directory.size() == 0) {
    LOG(FATAL)
        << "You must specifiy the directory of the 1dsfm dataset.";
  }

  InitializeFrom1DSFM(&reconstruction_builder);
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
