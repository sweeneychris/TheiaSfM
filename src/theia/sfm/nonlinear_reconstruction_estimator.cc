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

#include "theia/sfm/nonlinear_reconstruction_estimator.h"

#include <Eigen/Core>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/estimate_track.h"
#include "theia/sfm/filter_view_graph_cycles_by_rotation.h"
#include "theia/sfm/filter_view_pairs_from_orientation.h"
#include "theia/sfm/filter_view_pairs_from_relative_translation.h"
#include "theia/sfm/pose/estimate_positions_nonlinear.h"
#include "theia/sfm/pose/estimate_rotations_robust.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/set_camera_intrinsics_from_priors.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/view_graph/orientations_from_view_graph.h"
#include "theia/sfm/view_graph/remove_disconnected_view_pairs.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/random.h"
#include "theia/util/timer.h"

namespace theia {

using Eigen::Vector3d;

namespace {

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

// Sets the nonlinear position estimation options from the reconstruction
// estimator options.
NonlinearPositionEstimatorOptions SetNonlinearPositionEstimatorOptions(
    const ReconstructionEstimatorOptions& options) {
  NonlinearPositionEstimatorOptions npe_options;
  npe_options.num_threads = options.num_threads;
  npe_options.robust_loss_width =
      options.position_estimation_robust_loss_scale;
  npe_options.min_num_points_per_view =
      options.position_estimation_min_num_tracks_per_view;
  npe_options.point_to_camera_weight =
      options.position_estimation_point_to_camera_weight;
  return npe_options;
}

ViewId RandomViewId(const ViewGraph& view_graph) {
  const auto& view_pairs = view_graph.GetAllEdges();

  // Collect all view ids.
  std::unordered_set<ViewId> views;
  for (const auto& view_pair : view_pairs) {
    views.insert(view_pair.first.first);
    views.insert(view_pair.first.second);
  }

  // Find a random view id. TODO(cmsweeney): Choose the "best" random view by
  // some criterion such as highest connectivity.
  InitRandomGenerator();
  const int num_advances = RandInt(0, views.size() - 1);
  auto it = views.begin();
  std::advance(it, num_advances);
  return *it;
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

NonlinearReconstructionEstimator::NonlinearReconstructionEstimator(
    const ReconstructionEstimatorOptions& options) {
  options_ = options;
  translation_filter_options_ = SetRelativeTranslationFilteringOptions(options);
  position_estimator_options_ = SetNonlinearPositionEstimatorOptions(options);
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
//   6) Filter potentially bad relative translations with 1D SfM.
//   7) Estimate positions.
//   8) Estimate structure.
//   9) Bundle adjustment.
//   10) TODO: localize any unestimated cameras, retriangulate, and bundle
//      adjust.
//
// After each filtering step we remove any views which are no longer connected
// to the largest connected component in the view graph.
ReconstructionEstimatorSummary NonlinearReconstructionEstimator::Estimate(
    ViewGraph* view_graph, Reconstruction* reconstruction) {
  CHECK_NOTNULL(reconstruction);
  reconstruction_ = reconstruction;
  view_graph_ = view_graph;
  orientations_.clear();
  positions_.clear();

  ReconstructionEstimatorSummary summary;
  Timer timer;

  // Step 1. Filter the initial view graph and remove any bad two view
  // geometries.
  LOG(INFO) << "Filtering the intial view graph.";
  timer.Reset();
  if (!FilterInitialViewGraph()) {
    LOG(INFO) << "Insufficient view pairs to perform estimation.";
    return summary;
  }
  summary.initial_view_graph_filtering_time = timer.ElapsedTimeInSeconds();

  // Step 2. Calibrate any uncalibrated cameras.
  LOG(INFO) << "Calibrating any uncalibrated cameras.";
  timer.Reset();
  CalibrateCameras();
  summary.camera_intrinsics_calibration_time = timer.ElapsedTimeInSeconds();

  // Step 3. Estimate global rotations.
  LOG(INFO) << "Estimating the global rotations of all cameras.";
  timer.Reset();
  EstimateGlobalRotations();
  summary.rotation_estimation_time = timer.ElapsedTimeInSeconds();

  // Step 4. Filter bad rotations.
  LOG(INFO) << "Filtering any bad rotation estimations.";
  timer.Reset();
  FilterRotations();
  summary.rotation_filtering_time = timer.ElapsedTimeInSeconds();

  // Step 5. Optimize relative translations.
  LOG(INFO) << "Optimizing the pairwise translation estimations.";
  timer.Reset();
  OptimizePairwiseTranslations();
  summary.relative_translation_optimization_time = timer.ElapsedTimeInSeconds();

  // Step 6. Filter bad relative translations.
  LOG(INFO) << "Filtering any bad relative translations.";
  timer.Reset();
  FilterRelativeTranslation();
  summary.relative_translation_filtering_time = timer.ElapsedTimeInSeconds();

  // Step 7. Estimate global positions.
  LOG(INFO) << "Estimating the positions of all cameras.";
  timer.Reset();
  EstimatePosition();
  summary.position_estimation_time = timer.ElapsedTimeInSeconds();

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
    LOG(INFO) << "Performing bundle adjustment again.";
    timer.Reset();
    BundleAdjustment();
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
  return summary;
}

bool NonlinearReconstructionEstimator::FilterInitialViewGraph() {
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
  return view_graph_->NumEdges() >= 2;
}

void NonlinearReconstructionEstimator::CalibrateCameras() {
  SetCameraIntrinsicsFromPriors(reconstruction_);
}

void NonlinearReconstructionEstimator::EstimateGlobalRotations() {
  const ViewId random_starting_view = RandomViewId(*view_graph_);
  OrientationsFromViewGraph(*view_graph_, random_starting_view, &orientations_);
  const auto& relative_rotations = RelativeRotationsFromViewGraph(*view_graph_);
  RobustRotationEstimator::Options rotation_estimator_options;
  RobustRotationEstimator rotation_estimator(rotation_estimator_options,
                                             relative_rotations);
  CHECK(rotation_estimator.EstimateRotations(&orientations_))
      << "Could not estimate rotations.";
}

void NonlinearReconstructionEstimator::FilterRotations() {
  // Filter view pairs based on the relative rotation and the estimated global
  // orientations.
  FilterViewPairsFromOrientation(
      orientations_,
      options_.rotation_filtering_max_difference_degrees,
      view_graph_);
  RemoveDisconnectedViewPairs(view_graph_);
}

void NonlinearReconstructionEstimator::OptimizePairwiseTranslations() {
  if (options_.refine_relative_translations_after_rotation_estimation) {
    RefineRelativeTranslationsWithKnownRotations(*reconstruction_,
                                                 orientations_,
                                                 options_.num_threads,
                                                 view_graph_);
  }
}

void NonlinearReconstructionEstimator::FilterRelativeTranslation() {
  // Filter potentially bad relative translations.
  FilterViewPairsFromRelativeTranslation(translation_filter_options_,
                                         orientations_,
                                         view_graph_);
  RemoveDisconnectedViewPairs(view_graph_);
}

void NonlinearReconstructionEstimator::EstimatePosition() {
  // Estimate position.
  const auto& view_pairs = view_graph_->GetAllEdges();
  NonlinearPositionEstimator position_estimator(position_estimator_options_,
                                                *reconstruction_,
                                                view_pairs);
  CHECK(position_estimator.EstimatePositions(orientations_, &positions_))
      << "Position estimation failed!";
  LOG(INFO) << positions_.size()
            << " camera positions were estimated successfully.";
}

void NonlinearReconstructionEstimator::EstimateStructure() {
  // Estimate all tracks.
  EstimateTrackOptions triangulation_options;
  triangulation_options.max_acceptable_reprojection_error_pixels =
      options_.triangulation_max_reprojection_error_in_pixels;
  triangulation_options.min_triangulation_angle_degrees =
      options_.min_triangulation_angle_degrees;
  triangulation_options.bundle_adjustment = options_.bundle_adjust_tracks;
  EstimateAllTracks(triangulation_options,
                    options_.num_threads,
                    reconstruction_);
}

void NonlinearReconstructionEstimator::BundleAdjustment() {
  // Bundle adjustment.
  bundle_adjustment_options_ =
      SetBundleAdjustmentOptions(options_, positions_.size());
  const auto& bundle_adjustment_summary =
      BundleAdjustReconstruction(bundle_adjustment_options_, reconstruction_);
  CHECK(bundle_adjustment_summary.success) << "Could not perform BA.";
}

}  // namespace theia
