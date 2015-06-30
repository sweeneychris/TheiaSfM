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

#include <glog/logging.h>

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/bundle_adjustment/bundle_adjust_view.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/estimate_homography.h"
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
#include "theia/util/timer.h"
#include "theia/util/util.h"

namespace theia {
namespace {

// Given a viewing graph, extract the top k edges that have the most inliers
// between them.
void GetCandidateViewPairsByInlierCount(
    const ViewGraph& view_graph,
    const int num_candidate_view_pairs,
    std::vector<ViewIdPair>* candidate_view_pairs) {
  candidate_view_pairs->reserve(num_candidate_view_pairs);

  const auto& view_pairs = view_graph.GetAllEdges();

  // Collect the number of inliers for each view pair.
  std::vector<std::pair<int, ViewIdPair> > num_inliers_for_view_pair;
  num_inliers_for_view_pair.reserve(view_pairs.size());
  for (const auto& view_pair : view_pairs) {
    // TODO(cmsweeney): Prefer view pairs with known intrinsics.
    num_inliers_for_view_pair.emplace_back(
          view_pair.second.num_verified_matches, view_pair.first);
  }

  // Obtain the top k view pairs by number of inliers.
  const int num_elements_to_sort =
      std::min(num_candidate_view_pairs,
               static_cast<int>(num_inliers_for_view_pair.size()));
  std::partial_sort(num_inliers_for_view_pair.begin(),
                    num_inliers_for_view_pair.begin() + num_elements_to_sort,
                    num_inliers_for_view_pair.end(),
                    std::greater<std::pair<int, ViewIdPair> >());

  // Output the top k candidates
  for (int i = 0; i < num_candidate_view_pairs; i++) {
    candidate_view_pairs->emplace_back(num_inliers_for_view_pair[i].second);
  }
}

void NormalizeFeature(const Camera& camera, Feature* feature) {
  feature->y() =
      (feature->y() - camera.PrincipalPointY()) / camera.AspectRatio();
  feature->x() =
      feature->x() - camera.Skew() * feature->y() - camera.PrincipalPointX();
}

// Normalize a feature by its camera intrinsics.
void GetNormalizedFeatures(const Camera& camera1,
                           const Camera& camera2,
                           std::vector<FeatureCorrespondence>* matches) {
  for (FeatureCorrespondence& match : *matches) {
    NormalizeFeature(camera1, &match.feature1);
    NormalizeFeature(camera2, &match.feature2);
  }
}

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
  CHECK_GE(options.bundle_adjustment_growth_percent, 0.0)
      << "The bundle adjustment growth percent must be greater than 0 percent.";

  options_ = options;
  ransac_params_ = SetRansacParameters(options);

  // Triangulation options.
  triangulation_options_.max_acceptable_reprojection_error_pixels =
      options_.triangulation_max_reprojection_error_in_pixels;
  triangulation_options_.min_triangulation_angle_degrees =
      options_.min_triangulation_angle_degrees;
  triangulation_options_.bundle_adjustment = options_.bundle_adjust_tracks;

  // Localization options.
  localization_options_.reprojection_error_threshold_pixels =
      options_.absolute_pose_reprojection_error_threshold;
  localization_options_.ransac_params = ransac_params_;
  localization_options_.bundle_adjust_view = true;
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
// methods; hwoever, it requires many more instances of bundle adjustment (which
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

  // Set the known camera intrinsics.
  timer_.Reset();
  SetCameraIntrinsicsFromPriors(reconstruction_);
  summary_.camera_intrinsics_calibration_time = timer_.ElapsedTimeInSeconds();

  // Steps 1 - 3: Choose an initial camera pair to reconstruct.
  if (!ChooseInitialViewPair()) {
    LOG(ERROR) << "Could not find a suitable initial pair for starting "
                  "incremental SfM!";
    summary_.success = false;
    return summary_;
  }

  // Repeatedly try to localize new views until no more can be successfully
  // localized.
  const double ba_growth_ratio =
      1.0 + options_.bundle_adjustment_growth_percent / 100.0;

  // Try to add as many views as possible to the reconstruction until no more
  // views can be localized.
  while (!views_to_localize_.empty()) {
    // Step 4: Localize new views.
    if (!AddViewsToReconstruction()) {
      LOG(INFO) << "Failed to localize any remaining views.";
      break;
    }

    // Step 5: Estimate new 3D structure.
    EstimateStructure();

    // Step 6: Perform bundle adjustment on the reconstruction if the model has
    // grown by more than a given ratio.
    const int num_estimated_views =
        reconstruction_->NumViews() - views_to_localize_.size();
    const double reconstruction_growth =
        static_cast<double>(num_estimated_views) /
        static_cast<double>(num_optimized_views_);

    LOG(INFO) << "Num estimated views = " << num_estimated_views;
    LOG(INFO) << "Reconstruction growth = " << reconstruction_growth;
    if (reconstruction_growth > ba_growth_ratio) {
      // Set underconstrained tracks and views as unestimated and update the
      // points observed per view.
      SetUnderconstrainedAsUnestimated(reconstruction_);
      BundleAdjustment();
      num_optimized_views_ = num_estimated_views;
    }
  }

  // Run a final BA if any views are still unoptimized.
  if (reconstruction_->NumViews() - views_to_localize_.size() >
      num_optimized_views_) {
    BundleAdjustment();
  }

  // Set the output parameters.
  GetEstimatedViewsFromReconstruction(*reconstruction_,
                                      &summary_.estimated_views);
  GetEstimatedTracksFromReconstruction(*reconstruction_,
                                       &summary_.estimated_tracks);
  summary_.success = true;
  summary_.total_time = total_timer.ElapsedTimeInSeconds();
  return summary_;
}

void IncrementalReconstructionEstimator::InitializeCamerasFromTwoViewInfo(
    const ViewIdPair& view_ids) {
  View* view1 = reconstruction_->MutableView(view_ids.first);
  View* view2 = reconstruction_->MutableView(view_ids.second);
  const TwoViewInfo info =
      *view_graph_->GetEdge(view_ids.first, view_ids.second);

  // The pose of camera 1 will be the identity pose, so we only need to set
  // camera 2's pose.
  Camera* camera2 = view2->MutableCamera();
  camera2->SetOrientationFromAngleAxis(info.rotation_2);
  camera2->SetPosition(info.position_2);

  view1->SetEstimated(true);
  view2->SetEstimated(true);
}

bool IncrementalReconstructionEstimator::ChooseInitialViewPair() {
  static const int kNumCandidateViewPairs = 10;

  // Find the top k view pairs by the number of geometrically verified matches.
  std::vector<ViewIdPair> candidate_initial_view_pairs;
  GetCandidateViewPairsByInlierCount(*view_graph_,
                                     kNumCandidateViewPairs,
                                     &candidate_initial_view_pairs);

  OrderViewPairsByHomographyRatio(&candidate_initial_view_pairs);

  // Try to initialize the reconstruction from the candidate view pairs. An
  // initial seed is only considered valid if the baseline relative to the 3D
  // point depths is sufficient. This robustness is measured by the angle of all
  // 3D points.
  for (const ViewIdPair view_id_pair : candidate_initial_view_pairs) {
    // Initialize the camera poses of the intiial views and set the two views to
    // estimated.
    InitializeCamerasFromTwoViewInfo(view_id_pair);

    // Estimate 3D structure of the scene.
    EstimateStructure();

    // If we triangulated enough 3D points then run bundle adjustment and
    // return. Otherwise, try the next view pair as the seed for the initial
    // reconstruction.
    std::unordered_set<TrackId> estimated_tracks;
    GetEstimatedTracksFromReconstruction(*reconstruction_, &estimated_tracks);

    if (estimated_tracks.size() > options_.min_num_absolute_pose_inliers) {
      // Bundle adjustment on the 2-view reconstruction.
      BundleAdjustment();
      num_optimized_views_ = 2;
      views_to_localize_.erase(view_id_pair.first);
      views_to_localize_.erase(view_id_pair.second);
      return true;
    }

    // Set all values as unestimated and try to use the next candidate pair.
    SetReconstructionAsUnestimated(reconstruction_);
  }

  return false;
}

void IncrementalReconstructionEstimator::OrderViewPairsByHomographyRatio(
    std::vector<ViewIdPair>* view_id_pairs) {
  RansacParameters homography_ransac_params = ransac_params_;
  homography_ransac_params.error_thresh = options_.max_homography_error_pixels *
                                          options_.max_homography_error_pixels;

  // This container holds the ratio:
  //    num homography inliers / num essential matrix inliers
  // which is the ratio we want to minimize when finding an initial pair for
  // incremental SfM.
  std::vector<std::pair<double, ViewIdPair> > homography_ratios;
  homography_ratios.reserve(view_id_pairs->size());
  for (const ViewIdPair& view_id_pair : *view_id_pairs) {
    const View* view1 = reconstruction_->View(view_id_pair.first);
    const View* view2 = reconstruction_->View(view_id_pair.second);

    // Get tracks common to both views.
    const std::vector<ViewId> view_ids = {view_id_pair.first,
                                          view_id_pair.second};
    const std::vector<TrackId> track_ids =
        FindCommonTracksInViews(*reconstruction_, view_ids);

    // Normalize the features for both views.
    std::vector<FeatureCorrespondence> matches(track_ids.size());
    for (int i = 0; i < matches.size(); i++) {
      matches[i].feature1 = *view1->GetFeature(track_ids[i]);
      matches[i].feature2 = *view2->GetFeature(track_ids[i]);
    }
    GetNormalizedFeatures(view1->Camera(), view2->Camera(), &matches);

    // Estimate homography.
    Eigen::Matrix3d unused_homography;
    RansacSummary summary;
    EstimateHomography(homography_ransac_params,
                       RansacType::RANSAC,
                       matches,
                       &unused_homography,
                       &summary);

    const int num_ematrix_inliers =
        view_graph_->GetEdge(view_id_pair.first, view_id_pair.second)
            ->num_verified_matches;
    const double homography_ratio =
        static_cast<double>(summary.inliers.size()) /
        static_cast<double>(num_ematrix_inliers);
    homography_ratios.emplace_back(homography_ratio, view_id_pair);
  }

  // Sort the ratios in ascending order. This puts the most desirable view id
  // pair at the front.
  std::sort(homography_ratios.begin(), homography_ratios.end());

  // Output the view ids in the new ordering.
  for (int i = 0; i < view_id_pairs->size(); i++) {
    (*view_id_pairs)[i] = homography_ratios[i].second;
  }
}

bool IncrementalReconstructionEstimator::AddViewsToReconstruction() {
  RansacSummary unused_ransac_summary;

  // Compute the 2D-3D point count to determine which views should be localized.
  std::vector<ViewId> views_to_localize;
  FindViewsToLocalize(&views_to_localize);

  // Localize each view to determine the camera poses.
  if (reconstruction_->NumViews() - views_to_localize_.size() >
      num_optimized_views_) {
    BundleAdjustment();
  }


  int num_successful_localizations = 0;
  for (int i = 0; i < views_to_localize.size(); i++) {
    if (LocalizeViewToReconstruction(views_to_localize[i],
                                     localization_options_,
                                     reconstruction_,
                                     &unused_ransac_summary)) {
      views_to_localize_.erase(views_to_localize[i]);
      ++num_successful_localizations;
    }
  }

  return num_successful_localizations > 0;
}

void IncrementalReconstructionEstimator::FindViewsToLocalize(
    std::vector<ViewId>* views_to_localize) {
  // Determine the number of estimated tracks that each view observes.
  std::vector<std::pair<int, ViewId> > track_count_for_view;
  track_count_for_view.reserve(views_to_localize_.size());
  for (const ViewId view_id : views_to_localize_) {
    // Do not consider estimated views since they have already been localized.
    const View* view = reconstruction_->View(view_id);
    if (view->IsEstimated()) {
      continue;
    }

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

  // If M is the maximum number of 3D points observed by a view, we want to
  // localize all views that observe > M * localization_ratio 3D points.
  const int min_num_points_for_localization =
      static_cast<int>(options_.multiple_view_localization_ratio *
                       track_count_for_view[0].first);

  for (const auto& track_count : track_count_for_view) {
    if (track_count.first < min_num_points_for_localization) {
      break;
    }
    views_to_localize->emplace_back(track_count.second);
  }
}

void IncrementalReconstructionEstimator::EstimateStructure() {
  // Estimate all tracks.
  timer_.Reset();
  EstimateAllTracks(triangulation_options_,
                    options_.num_threads,
                    reconstruction_);
  summary_.triangulation_time += timer_.ElapsedTimeInSeconds();
}

void IncrementalReconstructionEstimator::BundleAdjustment() {
  timer_.Reset();
  const int num_estimated_views =
      reconstruction_->NumViews() - views_to_localize_.size();
  bundle_adjustment_options_ =
      SetBundleAdjustmentOptions(options_, num_estimated_views);

  // Inner iterations are not really needed for incremental SfM because we are
  // *hopefully* already starting at a good local minima. Inner iterations are
  // disabled because they slow down BA a lot.
  bundle_adjustment_options_.use_inner_iterations = false;

  const auto& bundle_adjustment_summary =
      BundleAdjustReconstruction(bundle_adjustment_options_, reconstruction_);
  CHECK(bundle_adjustment_summary.success) << "Could not perform BA.";
  summary_.bundle_adjustment_time += timer_.ElapsedTimeInSeconds();

  // Remove the outlier points based on the reprojection error and how
  // well-constrained the 3D points are.
  int num_points_removed = RemoveOutlierFeatures(
      options_.max_reprojection_error_in_pixels,
      options_.min_triangulation_angle_degrees,
      reconstruction_);
  LOG(INFO) << num_points_removed << " outlier points were removed.";

}

void IncrementalReconstructionEstimator::SetUnderconstrainedAsUnestimated(
    Reconstruction* reconstruction) {
  int num_underconstrained_views = -1;
  int num_underconstrained_tracks = -1;
  while (num_underconstrained_views != 0 && num_underconstrained_tracks != 0) {
    num_underconstrained_views =
        SetUnderconstrainedViewsToUnestimated(reconstruction);
    num_underconstrained_tracks =
        SetUnderconstrainedTracksToUnestimated(reconstruction);
  }

  // If any views were removed then we need to update the localization container
  // so that we can try to re-estimate the view.
  if (num_underconstrained_views > 0) {
    const auto& view_ids = reconstruction_->ViewIds();
    for (const ViewId view_id : view_ids) {
      if (!reconstruction_->View(view_id)->IsEstimated()) {
        views_to_localize_.insert(view_id);
      }
    }
  }
}

}  // namespace theia
