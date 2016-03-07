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

#include "theia/sfm/incremental_reconstruction_estimator.h"

#include <ceres/rotation.h>
#include <glog/logging.h>

#include <algorithm>
#include <functional>
#include <sstream>  // NOLINT
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/math/util.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/find_common_tracks_in_views.h"
#include "theia/sfm/localize_view_to_reconstruction.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/set_camera_intrinsics_from_priors.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/map_util.h"
#include "theia/util/stringprintf.h"
#include "theia/util/timer.h"
#include "theia/util/util.h"

namespace theia {
namespace {

void SetReconstructionAsUnestimated(Reconstruction* reconstruction) {
  // Set tracks as unestimated.
  const auto& track_ids = reconstruction->TrackIds();
  for (const TrackId track_id : track_ids) {
    Track* track = reconstruction->MutableTrack(track_id);
    if (track->IsEstimated()) {
      track->SetEstimated(false);
    }
  }

  // Set views as unestimated.
  const auto& view_ids = reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    View* view = reconstruction->MutableView(view_id);
    if (view->IsEstimated()) {
      view->SetEstimated(false);
    }
  }
}

}  // namespace

IncrementalReconstructionEstimator::IncrementalReconstructionEstimator(
    const ReconstructionEstimatorOptions& options) {
  CHECK_LE(options.multiple_view_localization_ratio, 1.0)
      << "The multiple view localization ratio must be between 0 and 1.0.";
  CHECK_GE(options.multiple_view_localization_ratio, 0.0)
      << "The multiple view localization ratio must be between 0 and 1.0.";
  CHECK_GE(options.full_bundle_adjustment_growth_percent, 0.0)
      << "The bundle adjustment growth percent must be greater than 0 percent.";
  CHECK_GE(options.partial_bundle_adjustment_num_views, 0)
      << "The bundle adjustment growth percent must be greater than 0 percent.";

  options_ = options;
  ransac_params_ = SetRansacParameters(options);

  // Triangulation options.
  triangulation_options_.max_acceptable_reprojection_error_pixels =
      options_.triangulation_max_reprojection_error_in_pixels;
  triangulation_options_.min_triangulation_angle_degrees =
      options_.min_triangulation_angle_degrees;
  triangulation_options_.bundle_adjustment = options_.bundle_adjust_tracks;
  triangulation_options_.ba_options = SetBundleAdjustmentOptions(options_, 0);
  triangulation_options_.ba_options.verbose = false;
  triangulation_options_.num_threads = options_.num_threads;

  // Localization options.
  localization_options_.reprojection_error_threshold_pixels =
      options_.absolute_pose_reprojection_error_threshold;
  localization_options_.ransac_params = ransac_params_;
  localization_options_.bundle_adjust_view = false;
  localization_options_.ba_options = SetBundleAdjustmentOptions(options_, 0);
  localization_options_.ba_options.verbose = false;
  localization_options_.min_num_inliers =
      options_.min_num_absolute_pose_inliers;

  num_optimized_views_ = 0;
}

// Estimates the camera position and 3D structure of the scene using an
// incremental Structure from Motion approach. The method begins by first
// estimating the 3D structure and camera poses of 2 cameras based on their
// relative pose. Then additional cameras are added on sequentially and new 3D
// structure is estimated as new parts of the scene are observed. Bundle
// adjustment is repeatedly performed as more cameras are added to ensure high
// quality reconstructions and to avoid drift.
//
// The incremental SfM pipeline is as follows:
//   1) Choose an initial camera pair to reconstruct.
//   2) Estimate 3D structure of the scene.
//   3) Bundle adjustment on the 2-view reconstruction.
//   4) Localize a new camera to the current 3D points. Choose the camera that
//      observes the most 3D points currently in the scene.
//   5) Estimate new 3D structure.
//   6) Bundle adjustment if the model has grown by more than 5% since the last
//      bundle adjustment.
//   7) Repeat steps 4-6 until all cameras have been added.
//
// Incremental SfM is generally considered to be more robust than global SfM
// methods; however, it requires many more instances of bundle adjustment (which
// is very costly) and so incremental SfM is not as efficient or scalable.
ReconstructionEstimatorSummary IncrementalReconstructionEstimator::Estimate(
    ViewGraph* view_graph, Reconstruction* reconstruction) {
  reconstruction_ = reconstruction;
  view_graph_ = view_graph;

  // Initialize the views_to_localize_ variable.
  const auto& view_ids = reconstruction_->ViewIds();
  views_to_localize_.reserve(view_ids.size());
  for (const ViewId view_id : view_ids) {
    views_to_localize_.insert(view_id);
  }

  Timer total_timer;
  Timer timer;
  double time_to_find_initial_seed = 0;

  // Set the known camera intrinsics.
  timer.Reset();
  SetCameraIntrinsicsFromPriors(reconstruction_);
  summary_.camera_intrinsics_calibration_time = timer.ElapsedTimeInSeconds();

  // Steps 1 - 3: Choose an initial camera pair to reconstruct.
  timer.Reset();
  if (!ChooseInitialViewPair()) {
    LOG(ERROR) << "Could not find a suitable initial pair for starting "
                  "incremental SfM!";
    summary_.success = false;
    return summary_;
  }
  time_to_find_initial_seed = timer.ElapsedTimeInSeconds();

  // Try to add as many views as possible to the reconstruction until no more
  // views can be localized.
  RansacSummary unused_ransac_summary;
  std::vector<ViewId> views_to_localize;
  int failed_localization_attempts = -1;
  while (!views_to_localize_.empty() &&
         failed_localization_attempts != views_to_localize.size()) {
    failed_localization_attempts = 0;
    views_to_localize.clear();

    // Step 4: Localize new views.
    // Compute the 2D-3D point count to determine which views should be
    // localized.
    timer.Reset();
    FindViewsToLocalize(&views_to_localize);
    summary_.pose_estimation_time += timer.ElapsedTimeInSeconds();

    // Attempt to localize all candidate views and estimate new 3D
    // points. Bundle Adjustment is run as either partial or full BA depending
    // on the current state of the reconstruction.
    for (int i = 0; i < views_to_localize.size(); i++) {
      timer.Reset();
      if (!LocalizeViewToReconstruction(views_to_localize[i],
                                        localization_options_,
                                        reconstruction_,
                                        &unused_ransac_summary)) {
        ++failed_localization_attempts;
        continue;
      }
      summary_.pose_estimation_time += timer.ElapsedTimeInSeconds();

      reconstructed_views_.push_back(views_to_localize[i]);
      views_to_localize_.erase(views_to_localize[i]);

      // Remove any tracks that have very bad 3D point reprojections after the
      // new view has been merged. This can happen when a new observation of a
      // 3D point has a very high reprojection error in the newly localized
      // view.
      const auto& tracks_in_new_view_vec =
          reconstruction_->View(reconstructed_views_.back())->TrackIds();
      const std::unordered_set<TrackId> tracks_in_new_view(
          tracks_in_new_view_vec.begin(), tracks_in_new_view_vec.end());
      RemoveOutlierTracks(
          tracks_in_new_view,
          triangulation_options_.max_acceptable_reprojection_error_pixels);

      // Step 5: Estimate new 3D points. and Step 6: Bundle adjustment.
      bool ba_success = false;
      if (UnoptimizedGrowthPercentage() <
          options_.full_bundle_adjustment_growth_percent) {
        // Step 5: Perform triangulation on the most recent view.
        timer.Reset();
        EstimateStructure(reconstructed_views_.back());
        summary_.triangulation_time += timer.ElapsedTimeInSeconds();

        // Step 6: Then perform partial Bundle Adjustment.
        timer.Reset();
        ba_success = PartialBundleAdjustment();
        summary_.bundle_adjustment_time += timer.ElapsedTimeInSeconds();
      } else {
        // Step 5: Perform triangulation on all views.
        timer.Reset();
        TrackEstimator track_estimator(triangulation_options_, reconstruction_);
        const TrackEstimator::Summary triangulation_summary =
            track_estimator.EstimateAllTracks();
        summary_.triangulation_time += timer.ElapsedTimeInSeconds();

        // Step 6: Full Bundle Adjustment.
        timer.Reset();
        ba_success = FullBundleAdjustment();
        summary_.bundle_adjustment_time += timer.ElapsedTimeInSeconds();
      }

      SetUnderconstrainedAsUnestimated();

      if (!ba_success) {
        LOG(WARNING) << "Bundle adjustment failed!";
        summary_.success = false;
        return summary_;
      }

    }
  }

  // Set the output parameters.
  GetEstimatedViewsFromReconstruction(*reconstruction_,
                                      &summary_.estimated_views);
  GetEstimatedTracksFromReconstruction(*reconstruction_,
                                       &summary_.estimated_tracks);
  summary_.success = true;
  summary_.total_time = total_timer.ElapsedTimeInSeconds();

  std::ostringstream string_stream;
  string_stream << "Incremental Reconstruction Estimator timings:"
                << "\n\tTime to find an initial seed for the reconstruction: "
                << time_to_find_initial_seed;
  summary_.message = string_stream.str();
  return summary_;
}

void IncrementalReconstructionEstimator::InitializeCamerasFromTwoViewInfo(
    const ViewIdPair& view_ids) {
  View* view1 = reconstruction_->MutableView(view_ids.first);
  View* view2 = reconstruction_->MutableView(view_ids.second);
  const TwoViewInfo* info =
      view_graph_->GetEdge(view_ids.first, view_ids.second);

  // The pose of camera 1 will be the identity pose, so we only need to set
  // camera 2's pose.
  Camera* camera1 = view1->MutableCamera();
  camera1->SetFocalLength(info->focal_length_1);
  Camera* camera2 = view2->MutableCamera();
  camera2->SetOrientationFromAngleAxis(info->rotation_2);
  camera2->SetPosition(info->position_2);
  camera2->SetFocalLength(info->focal_length_2);

  view1->SetEstimated(true);
  view2->SetEstimated(true);
}

bool IncrementalReconstructionEstimator::ChooseInitialViewPair() {
  static const int kMinNumInitialTracks = 100;

  // Sort the view pairs by the number of geometrically verified matches.
  std::vector<ViewIdPair> candidate_initial_view_pairs;
  OrderViewPairsByInitializationCriterion(kMinNumInitialTracks,
                                          &candidate_initial_view_pairs);

  if (candidate_initial_view_pairs.size() == 0) {
    return false;
  }

  // Find the k view pairs that contain the highest number of verified matches
  // and contain a sufficient baseline between them.
  OrderViewPairsByInitializationCriterion(kMinNumInitialTracks,
                                          &candidate_initial_view_pairs);

  // Try to initialize the reconstruction from the candidate view pairs. An
  // initial seed is only considered valid if the baseline relative to the 3D
  // point depths is sufficient. This robustness is measured by the angle of all
  // 3D points.
  for (const ViewIdPair view_id_pair : candidate_initial_view_pairs) {
    // Set all values as unestimated and try to use the next candidate pair.
    SetReconstructionAsUnestimated(reconstruction_);

    // Initialize the camera poses of the intiial views and set the two views to
    // estimated.
    InitializeCamerasFromTwoViewInfo(view_id_pair);

    // Estimate 3D structure of the scene.
    EstimateStructure(view_id_pair.first);

    // If we did not triangulate enough tracks then skip this view and try
    // another.
    std::unordered_set<TrackId> estimated_tracks;
    GetEstimatedTracksFromReconstruction(*reconstruction_, &estimated_tracks);

    if (estimated_tracks.size() < kMinNumInitialTracks) {
      continue;
    }

    // Bundle adjustment on the 2-view reconstruction.
    if (!FullBundleAdjustment()) {
      continue;
    }

    // If we triangulated enough 3D points then return.  Otherwise, try the next
    // view pair as the seed for the initial reconstruction.
    estimated_tracks.clear();
    GetEstimatedTracksFromReconstruction(*reconstruction_, &estimated_tracks);

    if (estimated_tracks.size() > kMinNumInitialTracks) {
      reconstructed_views_.push_back(view_id_pair.first);
      reconstructed_views_.push_back(view_id_pair.second);
      views_to_localize_.erase(view_id_pair.first);
      views_to_localize_.erase(view_id_pair.second);

      return true;
    }
  }

  return false;
}

void IncrementalReconstructionEstimator::
    OrderViewPairsByInitializationCriterion(
        const int min_num_verified_matches,
        std::vector<ViewIdPair>* view_id_pairs) {
  const auto& view_pairs = view_graph_->GetAllEdges();
  view_id_pairs->reserve(view_pairs.size());

  // Collect the number of inliers for each view pair. The tuples store:
  //     # homography inliers, negative of essential matrix inliers, ViewIdPair
  //
  // We store the negative of the essential matrix inliers because we want a
  // view pair with the fewest homography inliers but the greatest number of
  // essential matrix inliers. This situation only arises if there were a
  // tiebreaker in the number of homography inliers, or if (for some unknown
  // reason) the number of homography inliers is set to 0 for all view pairs.
  std::vector<std::tuple<int, int, ViewIdPair> >
      initialization_criterion_for_view_pairs;
  initialization_criterion_for_view_pairs.reserve(view_pairs.size());
  for (const auto& view_pair : view_pairs) {
    // TODO(cmsweeney): Prefer view pairs with known intrinsics.
    if (view_pair.second.num_verified_matches > min_num_verified_matches) {
      initialization_criterion_for_view_pairs.emplace_back(
          view_pair.second.num_homography_inliers,
          -view_pair.second.num_verified_matches,
          view_pair.first);
    }
  }

  // Sort the views to find the ones that are least well-modelled by a
  // homography.
  std::sort(initialization_criterion_for_view_pairs.begin(),
            initialization_criterion_for_view_pairs.end());
  for (int i = 0; i < initialization_criterion_for_view_pairs.size(); i++) {
    view_id_pairs->emplace_back(
        std::get<2>(initialization_criterion_for_view_pairs[i]));
  }
}

void IncrementalReconstructionEstimator::FindViewsToLocalize(
    std::vector<ViewId>* views_to_localize) {
  // We localize all views that observe 75% or more than the number of 3D points
  // observed by the view with the largest number of observed 3D points.
  static const double kObserved3dPointsRatio = 0.75;

  // Determine the number of estimated tracks that each view observes.
  std::vector<std::pair<int, ViewId> > track_count_for_view;
  track_count_for_view.reserve(views_to_localize_.size());
  for (const ViewId view_id : views_to_localize_) {
    // Do not consider estimated views since they have already been localized.
    const View* view = reconstruction_->View(view_id);

    // Count the number of estimated tracks for this view.
    const auto& track_ids = view->TrackIds();
    int num_estimated_tracks = 0;
    for (const TrackId track_id : track_ids) {
      if (reconstruction_->Track(track_id)->IsEstimated()) {
        ++num_estimated_tracks;
      }
    }

    track_count_for_view.emplace_back(num_estimated_tracks, view_id);
  }

  // Sort the track count so that the view with the most tracks is at the front.
  std::sort(track_count_for_view.begin(),
            track_count_for_view.end(),
            std::greater<std::pair<int, ViewId> >());
  const int min_3d_points_observed = static_cast<int>(
      track_count_for_view.begin()->first * kObserved3dPointsRatio);
  for (const auto& track_count : track_count_for_view) {
    if (track_count.first < min_3d_points_observed ||
        track_count.first < options_.min_num_absolute_pose_inliers) {
      break;
    }
    views_to_localize->emplace_back(track_count.second);
  }
}

void IncrementalReconstructionEstimator::EstimateStructure(
    const ViewId view_id) {
  // Estimate all tracks.
  TrackEstimator track_estimator(triangulation_options_, reconstruction_);
  const std::vector<TrackId>& tracks_in_view =
      reconstruction_->View(view_id)->TrackIds();
  const std::unordered_set<TrackId> tracks_to_triangulate(
      tracks_in_view.begin(), tracks_in_view.end());
  const TrackEstimator::Summary summary =
      track_estimator.EstimateTracks(tracks_to_triangulate);
}

double IncrementalReconstructionEstimator::UnoptimizedGrowthPercentage() {
  return 100.0 * (reconstructed_views_.size() - num_optimized_views_) /
         static_cast<double>(num_optimized_views_);
}

bool IncrementalReconstructionEstimator::FullBundleAdjustment() {
  // Full bundle adjustment.
  LOG(INFO) << "Running full bundle adjustment on the entire reconstruction.";

  // Set up the BA options.
  bundle_adjustment_options_ =
      SetBundleAdjustmentOptions(options_, reconstructed_views_.size());

  // Inner iterations are not really needed for incremental SfM because we are
  // *hopefully* already starting at a good local minima. Inner iterations are
  // disabled because they slow down BA a lot.
  bundle_adjustment_options_.use_inner_iterations = false;

  const BundleAdjustmentSummary ba_summary =
    BundleAdjustReconstruction(bundle_adjustment_options_, reconstruction_);
  num_optimized_views_ = reconstructed_views_.size();

  const auto& track_ids = reconstruction_->TrackIds();
  const std::unordered_set<TrackId> all_tracks(track_ids.begin(),
                                               track_ids.end());
  RemoveOutlierTracks(all_tracks, options_.max_reprojection_error_in_pixels);

  return ba_summary.success;
}

bool IncrementalReconstructionEstimator::PartialBundleAdjustment() {
// Partial bundle adjustment only only the k most recently added views that
  // have not been optimized by full BA.
  const int partial_ba_size =
    std::min(static_cast<int>(reconstructed_views_.size()),
             options_.partial_bundle_adjustment_num_views);
  LOG(INFO) << "Running partial bundle adjustment on " << partial_ba_size
            << " views.";

  // Set up the BA options.
  bundle_adjustment_options_ =
      SetBundleAdjustmentOptions(options_, partial_ba_size);

  // Inner iterations are not really needed for incremental SfM because we are
  // *hopefully* already starting at a good local minima. Inner iterations are
  // disabled because they slow down BA a lot.
  bundle_adjustment_options_.use_inner_iterations = false;
  bundle_adjustment_options_.verbose = VLOG_IS_ON(2);

  // If the model has grown sufficiently then run BA on the entire
  // model. Otherwise, run partial BA.
  BundleAdjustmentSummary ba_summary;

  // Get the views to optimize for partial BA.
  std::unordered_set<ViewId> views_to_optimize(
      reconstructed_views_.end() - partial_ba_size,
      reconstructed_views_.end());
  // Get the tracks observed in these views.
  std::unordered_set<TrackId> tracks_to_optimize;
  for (const ViewId view_to_optimize : views_to_optimize) {
    const View* view = reconstruction_->View(view_to_optimize);
    const auto& tracks_in_view = view->TrackIds();
    for (const TrackId track_in_view : tracks_in_view) {
      tracks_to_optimize.insert(track_in_view);
    }
  }

  ba_summary = BundleAdjustPartialReconstruction(bundle_adjustment_options_,
                                                 views_to_optimize,
                                                 tracks_to_optimize,
                                                 reconstruction_);

  RemoveOutlierTracks(tracks_to_optimize,
                      options_.max_reprojection_error_in_pixels);
  return ba_summary.success;
}

void IncrementalReconstructionEstimator::RemoveOutlierTracks(
    const std::unordered_set<TrackId>& tracks_to_check,
    const double max_reprojection_error_in_pixels) {
  // Remove the outlier points based on the reprojection error and how
  // well-constrained the 3D points are.
  int num_points_removed = RemoveOutlierFeatures(
      tracks_to_check,
      max_reprojection_error_in_pixels,
      options_.min_triangulation_angle_degrees,
      reconstruction_);
  LOG(INFO) << num_points_removed << " outlier points were removed.";
}

void IncrementalReconstructionEstimator::SetUnderconstrainedAsUnestimated() {
  int num_underconstrained_views = -1;
  int num_underconstrained_tracks = -1;
  while (num_underconstrained_views != 0 && num_underconstrained_tracks != 0) {
    num_underconstrained_views =
        SetUnderconstrainedViewsToUnestimated(reconstruction_);
    num_underconstrained_tracks =
        SetUnderconstrainedTracksToUnestimated(reconstruction_);
  }

  // If any views were removed then we need to update the localization container
  // so that we can try to re-estimate the view.
  if (num_underconstrained_views > 0) {
    const auto& view_ids = reconstruction_->ViewIds();
    for (const ViewId view_id : view_ids) {
      if (!reconstruction_->View(view_id)->IsEstimated() &&
          !ContainsKey(views_to_localize_, view_id)) {
        views_to_localize_.insert(view_id);

        // Remove the view from the list of localized views.
        auto view_to_remove = std::find(reconstructed_views_.begin(),
                                        reconstructed_views_.end(), view_id);
        reconstructed_views_.erase(view_to_remove);
        --num_optimized_views_;
      }
    }
  }
}

}  // namespace theia
