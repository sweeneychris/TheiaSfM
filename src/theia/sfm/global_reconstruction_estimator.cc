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

#include "theia/sfm/global_reconstruction_estimator.h"

#include <Eigen/Core>
#include <memory>
#include <sstream>  // NOLINT

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/estimate_track.h"
#include "theia/sfm/extract_maximally_parallel_rigid_subgraph.h"
#include "theia/sfm/filter_view_graph_cycles_by_rotation.h"
#include "theia/sfm/filter_view_pairs_from_orientation.h"
#include "theia/sfm/filter_view_pairs_from_relative_translation.h"
#include "theia/sfm/global_pose_estimation/least_unsquared_deviation_position_estimator.h"
#include "theia/sfm/global_pose_estimation/linear_position_estimator.h"
#include "theia/sfm/global_pose_estimation/linear_rotation_estimator.h"
#include "theia/sfm/global_pose_estimation/nonlinear_position_estimator.h"
#include "theia/sfm/global_pose_estimation/nonlinear_rotation_estimator.h"
#include "theia/sfm/global_pose_estimation/position_estimator.h"
#include "theia/sfm/global_pose_estimation/robust_rotation_estimator.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/set_camera_intrinsics_from_priors.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/view_graph/orientations_from_maximum_spanning_tree.h"
#include "theia/sfm/view_graph/remove_disconnected_view_pairs.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/random.h"
#include "theia/util/timer.h"

namespace theia {

using Eigen::Vector3d;

namespace {

// All times are given in seconds.
struct GlobalReconstructionEstimatorTimings {
  double initial_view_graph_filtering_time = 0.0;
  double camera_intrinsics_calibration_time = 0.0;
  double rotation_estimation_time = 0.0;
  double rotation_filtering_time = 0.0;
  double relative_translation_optimization_time = 0.0;
  double relative_translation_filtering_time = 0.0;
  double position_estimation_time = 0.0;
};

FilterViewPairsFromRelativeTranslationOptions
SetRelativeTranslationFilteringOptions(
    const ReconstructionEstimatorOptions& options) {
  FilterViewPairsFromRelativeTranslationOptions fvpfrt_options;
  fvpfrt_options.num_threads = options.num_threads;
  fvpfrt_options.num_iterations = options.translation_filtering_num_iterations;
  fvpfrt_options.translation_projection_tolerance =
      options.translation_filtering_projection_tolerance;
  return fvpfrt_options;
}

void SetUnderconstrainedAsUnestimated(Reconstruction* reconstruction) {
  int num_underconstrained_views = -1;
  int num_underconstrained_tracks = -1;
  while (num_underconstrained_views != 0 && num_underconstrained_tracks != 0) {
    num_underconstrained_views =
        SetUnderconstrainedViewsToUnestimated(reconstruction);
    num_underconstrained_tracks =
        SetUnderconstrainedTracksToUnestimated(reconstruction);
  }
}

}  // namespace

GlobalReconstructionEstimator::GlobalReconstructionEstimator(
    const ReconstructionEstimatorOptions& options) {
  options_ = options;
  translation_filter_options_ = SetRelativeTranslationFilteringOptions(options);
  options_.nonlinear_position_estimator_options.num_threads =
      options_.num_threads;
  options_.linear_triplet_position_estimator_options.num_threads =
      options_.num_threads;
  ransac_params_ = SetRansacParameters(options);
}

// The pipeline for estimating camera poses and structure is as follows:
//   1) Filter potentially bad pairwise geometries by enforcing a loop
//      constaint on rotations that form a triplet.
//   2) Initialize focal lengths.
//   3) Estimate the global rotation for each camera.
//   4) Remove any pairwise geometries where the relative rotation is not
//      consistent with the global rotation.
//   5) Optimize the relative translation given the known rotations.
//   6) Filter potentially bad relative translations.
//   7) Estimate positions.
//   8) Estimate structure.
//   9) Bundle adjustment.
//   10) Retriangulate, and bundle adjust.
//
// After each filtering step we remove any views which are no longer connected
// to the largest connected component in the view graph.
ReconstructionEstimatorSummary GlobalReconstructionEstimator::Estimate(
    ViewGraph* view_graph, Reconstruction* reconstruction) {
  CHECK_NOTNULL(reconstruction);
  reconstruction_ = reconstruction;
  view_graph_ = view_graph;
  orientations_.clear();
  positions_.clear();

  ReconstructionEstimatorSummary summary;
  GlobalReconstructionEstimatorTimings global_estimator_timings;
  Timer total_timer;
  Timer timer;

  // Step 1. Filter the initial view graph and remove any bad two view
  // geometries.
  LOG(INFO) << "Filtering the intial view graph.";
  timer.Reset();
  if (!FilterInitialViewGraph()) {
    LOG(INFO) << "Insufficient view pairs to perform estimation.";
    return summary;
  }
  global_estimator_timings.initial_view_graph_filtering_time =
      timer.ElapsedTimeInSeconds();

  // Step 2. Calibrate any uncalibrated cameras.
  LOG(INFO) << "Calibrating any uncalibrated cameras.";
  timer.Reset();
  CalibrateCameras();
  summary.camera_intrinsics_calibration_time = timer.ElapsedTimeInSeconds();

  // Step 3. Estimate global rotations.
  LOG(INFO) << "Estimating the global rotations of all cameras.";
  timer.Reset();
  if (!EstimateGlobalRotations()) {
    LOG(WARNING) << "Rotation estimation failed!";
    summary.success = false;
    return summary;
  }
  global_estimator_timings.rotation_estimation_time =
      timer.ElapsedTimeInSeconds();

  // Step 4. Filter bad rotations.
  LOG(INFO) << "Filtering any bad rotation estimations.";
  timer.Reset();
  FilterRotations();
  global_estimator_timings.rotation_filtering_time =
      timer.ElapsedTimeInSeconds();

  // Step 5. Optimize relative translations.
  LOG(INFO) << "Optimizing the pairwise translation estimations.";
  timer.Reset();
  OptimizePairwiseTranslations();
  global_estimator_timings.relative_translation_optimization_time =
      timer.ElapsedTimeInSeconds();

  // Step 6. Filter bad relative translations.
  LOG(INFO) << "Filtering any bad relative translations.";
  timer.Reset();
  FilterRelativeTranslation();
  global_estimator_timings.relative_translation_filtering_time =
      timer.ElapsedTimeInSeconds();

  // Step 7. Estimate global positions.
  LOG(INFO) << "Estimating the positions of all cameras.";
  timer.Reset();
  if (!EstimatePosition()) {
    LOG(WARNING) << "Position estimation failed!";
    summary.success = false;
    return summary;
  }
  LOG(INFO) << positions_.size()
            << " camera positions were estimated successfully.";
  global_estimator_timings.position_estimation_time =
      timer.ElapsedTimeInSeconds();

  summary.pose_estimation_time =
      global_estimator_timings.rotation_estimation_time +
      global_estimator_timings.rotation_filtering_time +
      global_estimator_timings.relative_translation_optimization_time +
      global_estimator_timings.relative_translation_filtering_time +
      global_estimator_timings.position_estimation_time;

  // Set the poses in the reconstruction object.
  SetReconstructionFromEstimatedPoses(orientations_,
                                      positions_,
                                      reconstruction_);

  // Always triangulate once, then retriangulate and remove outliers depending
  // on the reconstruciton estimator options.
  for (int i = 0; i < options_.num_retriangulation_iterations + 1; i++) {
    // Step 8. Triangulate features.
    LOG(INFO) << "Triangulating all features.";
    timer.Reset();
    EstimateStructure();
    summary.triangulation_time += timer.ElapsedTimeInSeconds();

    SetUnderconstrainedAsUnestimated(reconstruction_);

    // Step 9. Bundle Adjustment.
    LOG(INFO) << "Performing bundle adjustment.";
    timer.Reset();
    if (!BundleAdjustment()) {
      summary.success = false;
      LOG(WARNING) << "Bundle adjustment failed!";
      return summary;
    }
    summary.bundle_adjustment_time += timer.ElapsedTimeInSeconds();

    int num_points_removed = RemoveOutlierFeatures(
        options_.max_reprojection_error_in_pixels,
        options_.min_triangulation_angle_degrees,
        reconstruction_);
    LOG(INFO) << num_points_removed << " outlier points were removed.";
  }

  // Set the output parameters.
  GetEstimatedViewsFromReconstruction(*reconstruction_,
                                      &summary.estimated_views);
  GetEstimatedTracksFromReconstruction(*reconstruction_,
                                       &summary.estimated_tracks);
  summary.success = true;
  summary.total_time = total_timer.ElapsedTimeInSeconds();

  // Output some timing statistics.
  std::ostringstream string_stream;
  string_stream
      << "Global Reconstruction Estimator timings:"
      << "\n\tInitial view graph filtering time = "
      << global_estimator_timings.initial_view_graph_filtering_time
      << "\n\tCamera intrinsic calibration time = "
      << summary.camera_intrinsics_calibration_time
      << "\n\tRotation estimation time = "
      << global_estimator_timings.rotation_estimation_time
      << "\n\tRotation filtering time = "
      << global_estimator_timings.rotation_filtering_time
      << "\n\tRelative translation optimization time = "
      << global_estimator_timings.relative_translation_optimization_time
      << "\n\tRelative translation filtering time = "
      << global_estimator_timings.relative_translation_filtering_time
      << "\n\tPosition estimation time = "
      << global_estimator_timings.position_estimation_time;
  summary.message = string_stream.str();

  return summary;
}

bool GlobalReconstructionEstimator::FilterInitialViewGraph() {
  // Remove any view pairs that do not have a sufficient number of inliers.
  std::unordered_set<ViewIdPair> view_pairs_to_remove;
  const auto& view_pairs = view_graph_->GetAllEdges();
  for (const auto& view_pair : view_pairs) {
    if (view_pair.second.num_verified_matches <
        options_.min_num_two_view_inliers) {
      view_pairs_to_remove.insert(view_pair.first);
    }
  }
  for (const ViewIdPair view_id_pair : view_pairs_to_remove) {
    view_graph_->RemoveEdge(view_id_pair.first, view_id_pair.second);
  }

  // Only reconstruct the largest connected component.
  RemoveDisconnectedViewPairs(view_graph_);
  return view_graph_->NumEdges() >= 1;
}

void GlobalReconstructionEstimator::CalibrateCameras() {
  SetCameraIntrinsicsFromPriors(reconstruction_);
}

bool GlobalReconstructionEstimator::EstimateGlobalRotations() {
  const auto& view_pairs = view_graph_->GetAllEdges();

  // Choose the global rotation estimation type.
  std::unique_ptr<RotationEstimator> rotation_estimator;
  switch (options_.global_rotation_estimator_type) {
    case GlobalRotationEstimatorType::ROBUST_L1L2: {
      // Initialize the orientation estimations by walking along the maximum
      // spanning tree.
      OrientationsFromMaximumSpanningTree(*view_graph_, &orientations_);
      RobustRotationEstimator::Options robust_rotation_estimator_options;
      rotation_estimator.reset(
          new RobustRotationEstimator(robust_rotation_estimator_options));
      break;
    }
    case GlobalRotationEstimatorType::NONLINEAR: {
      // Initialize the orientation estimations by walking along the maximum
      // spanning tree.
      OrientationsFromMaximumSpanningTree(*view_graph_, &orientations_);
      rotation_estimator.reset(new NonlinearRotationEstimator());
      break;
    }
    case GlobalRotationEstimatorType::LINEAR: {
      // Set the constructor variable to true to weigh each term by the inlier
      // count.
      rotation_estimator.reset(new LinearRotationEstimator(false));
      break;
    }
    default: {
      LOG(FATAL) << "Invalid type of global rotation estimation chosen.";
      break;
    }
  }

  return rotation_estimator->EstimateRotations(view_pairs, &orientations_);
}

void GlobalReconstructionEstimator::FilterRotations() {
  // Filter view pairs based on the relative rotation and the estimated global
  // orientations.
  FilterViewPairsFromOrientation(
      orientations_,
      options_.rotation_filtering_max_difference_degrees,
      view_graph_);
  // Remove any disconnected views from the estimation.
  const std::unordered_set<ViewId> removed_views =
      RemoveDisconnectedViewPairs(view_graph_);
  for (const ViewId removed_view : removed_views) {
    orientations_.erase(removed_view);
  }
}

void GlobalReconstructionEstimator::OptimizePairwiseTranslations() {
  if (options_.refine_relative_translations_after_rotation_estimation) {
    RefineRelativeTranslationsWithKnownRotations(*reconstruction_,
                                                 orientations_,
                                                 options_.num_threads,
                                                 view_graph_);
  }
}

void GlobalReconstructionEstimator::FilterRelativeTranslation() {
  if (options_.extract_maximal_rigid_subgraph) {
    LOG(INFO) << "Extracting maximal rigid component of viewing graph to "
                 "determine which cameras are well-constrained for position "
                 "estimation.";
    ExtractMaximallyParallelRigidSubgraph(orientations_, view_graph_);
  }

  // Filter potentially bad relative translations.
  if (options_.filter_relative_translations_with_1dsfm) {
    LOG(INFO) << "Filtering relative translations with 1DSfM filter.";
    FilterViewPairsFromRelativeTranslation(translation_filter_options_,
                                           orientations_,
                                           view_graph_);
  }
  // Remove any disconnected views from the estimation.
  const std::unordered_set<ViewId> removed_views =
      RemoveDisconnectedViewPairs(view_graph_);
  for (const ViewId removed_view : removed_views) {
    orientations_.erase(removed_view);
  }
}

bool GlobalReconstructionEstimator::EstimatePosition() {
  // Estimate position.
  const auto& view_pairs = view_graph_->GetAllEdges();
  std::unique_ptr<PositionEstimator> position_estimator;

  // Choose the global position estimation type.
  switch (options_.global_position_estimator_type) {
    case GlobalPositionEstimatorType::NONLINEAR: {
      position_estimator.reset(new NonlinearPositionEstimator(
          options_.nonlinear_position_estimator_options, *reconstruction_));
      break;
    }
    case GlobalPositionEstimatorType::LINEAR_TRIPLET: {
      position_estimator.reset(new LinearPositionEstimator(
          options_.linear_triplet_position_estimator_options,
          *reconstruction_));
      break;
    }
    case GlobalPositionEstimatorType::LEAST_UNSQUARED_DEVIATION: {
      position_estimator.reset(new LeastUnsquaredDeviationPositionEstimator(
          options_.least_unsquared_deviation_position_estimator_options));
      break;
    }
    default: {
      LOG(FATAL) << "Invalid type of global position estimation chosen.";
      break;
    }
  }

  return position_estimator->EstimatePositions(view_pairs,
                                               orientations_,
                                               &positions_);
}

void GlobalReconstructionEstimator::EstimateStructure() {
  // Estimate all tracks.
  TrackEstimator::Options triangulation_options;
  triangulation_options.max_acceptable_reprojection_error_pixels =
      options_.triangulation_max_reprojection_error_in_pixels;
  triangulation_options.min_triangulation_angle_degrees =
      options_.min_triangulation_angle_degrees;
  triangulation_options.bundle_adjustment = options_.bundle_adjust_tracks;
  triangulation_options.ba_options = SetBundleAdjustmentOptions(options_, 0);
  triangulation_options.ba_options.verbose = false;
  triangulation_options.num_threads = options_.num_threads;
  TrackEstimator track_estimator(triangulation_options, reconstruction_);
  const TrackEstimator::Summary summary = track_estimator.EstimateAllTracks();
}

bool GlobalReconstructionEstimator::BundleAdjustment() {
  // Bundle adjustment.
  bundle_adjustment_options_ =
      SetBundleAdjustmentOptions(options_, positions_.size());
  const auto& bundle_adjustment_summary =
      BundleAdjustReconstruction(bundle_adjustment_options_, reconstruction_);
  return bundle_adjustment_summary.success;
}

}  // namespace theia
