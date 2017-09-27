// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#include "theia/sfm/hybrid_reconstruction_estimator.h"

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
#include "theia/sfm/estimators/estimate_relative_pose_with_known_orientation.h"
#include "theia/sfm/find_common_tracks_in_views.h"
#include "theia/sfm/global_pose_estimation/linear_rotation_estimator.h"
#include "theia/sfm/global_pose_estimation/nonlinear_rotation_estimator.h"
#include "theia/sfm/global_pose_estimation/robust_rotation_estimator.h"
#include "theia/sfm/global_pose_estimation/rotation_estimator.h"
#include "theia/sfm/localize_view_to_reconstruction.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/select_good_tracks_for_bundle_adjustment.h"
#include "theia/sfm/set_camera_intrinsics_from_priors.h"
#include "theia/sfm/set_outlier_tracks_to_unestimated.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/orientations_from_maximum_spanning_tree.h"
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

HybridReconstructionEstimator::HybridReconstructionEstimator(
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
  triangulation_options_.ba_options.num_threads = 1;
  triangulation_options_.ba_options.verbose = false;
  triangulation_options_.num_threads = options_.num_threads;

  // Localization options.
  localization_options_.reprojection_error_threshold_pixels =
      options_.absolute_pose_reprojection_error_threshold;
  localization_options_.ransac_params = ransac_params_;
  localization_options_.bundle_adjust_view = false;
  localization_options_.ba_options = SetBundleAdjustmentOptions(options_, 0);
  localization_options_.ba_options.constant_camera_orientation = true;
  localization_options_.ba_options.verbose = false;
  localization_options_.min_num_inliers =
      options_.min_num_absolute_pose_inliers;

  num_optimized_views_ = 0;
}

ReconstructionEstimatorSummary HybridReconstructionEstimator::Estimate(
    ViewGraph* view_graph, Reconstruction* reconstruction) {
  reconstruction_ = reconstruction;
  view_graph_ = view_graph;

  // Initialize the unlocalized_views_ variable.
  const auto& view_ids = reconstruction_->ViewIds();
  unlocalized_views_.reserve(view_ids.size());
  for (const ViewId view_id : view_ids) {
    const View* view = reconstruction_->View(view_id);
    if (!view->IsEstimated()) {
      unlocalized_views_.insert(view_id);
    }
  }

  Timer total_timer;
  Timer timer;
  double time_to_find_initial_seed = 0;

  // Set the known camera intrinsics.
  timer.Reset();
  SetCameraIntrinsicsFromPriors(reconstruction_);
  summary_.camera_intrinsics_calibration_time = timer.ElapsedTimeInSeconds();

  // Step 1: Estimate camera orientations using a global rotation estimator.
  if (!EstimateCameraOrientations()) {
    LOG(ERROR) << "Could not estimate camera rotations for Hybrid SfM.";
    summary_.success = false;
    return summary_;
  }

  // Steps 2 - 3: Choose an initial camera pair to reconstruct if the
  // reconstruction is not already initialized.
  const int num_estimated_tracks = 0;
  const int num_estimated_views = 0;
  timer.Reset();
  LOG(INFO) << "Attempting to initialize the camera translation estimation.";
  if (!ChooseInitialViewPair()) {
    LOG(ERROR) << "Could not find a suitable initial pair for starting "
        "hybrid SfM!";
    summary_.success = false;
    return summary_;
  }
  time_to_find_initial_seed = timer.ElapsedTimeInSeconds();

  // Try to add as many views as possible to the reconstruction until no more
  // views can be localized.
  std::vector<ViewId> views_to_localize;
  int failed_localization_attempts = -1;
  while (!unlocalized_views_.empty() &&
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

      // Localize the view to the reconstruciton. If the orientation was
      // estimated from the global algorithm, this will first try to use a
      // simplified solver to estimate the camera position assuming the known
      // orientation.
      if (!LocalizeView(views_to_localize[i])) {
        ++failed_localization_attempts;
        continue;
      }
      summary_.pose_estimation_time += timer.ElapsedTimeInSeconds();

      reconstructed_views_.push_back(views_to_localize[i]);
      unlocalized_views_.erase(views_to_localize[i]);

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
  string_stream << "Hybrid Reconstruction Estimator timings:"
                << "\n\tTime to find an initial seed for the reconstruction: "
                << time_to_find_initial_seed;
  summary_.message = string_stream.str();
  return summary_;
}

bool HybridReconstructionEstimator::LocalizeView(const ViewId view_id) {
  if (ContainsKey(orientations_, view_id)) {
    localization_options_.assume_known_orientation = true;
    RansacSummary unused_ransac_summary;
    if (LocalizeViewToReconstruction(view_id,
                                     localization_options_,
                                     reconstruction_,
                                     &unused_ransac_summary)) {
      return true;
    }
  }

  // If we reached here, then either the orientation of this view was not
  // computed during global orientation estimation or the localization of
  // only the position failed.
  localization_options_.assume_known_orientation = false;
  RansacSummary unused_ransac_summary;
  return LocalizeViewToReconstruction(view_id,
                                      localization_options_,
                                      reconstruction_,
                                      &unused_ransac_summary);
}

bool HybridReconstructionEstimator::EstimateCameraOrientations() {
  // TODO(csweeney): Currently we use all view pairs to estimate the orientation
  // for all possible cameras. This ignores any information about views that are
  // already estimated, which should instead be exposed to improve the
  // estimation of the unknown camera orientations.
  const auto& view_pairs = view_graph_->GetAllEdges();

  // Choose the global rotation estimation type.
  std::unique_ptr<RotationEstimator> rotation_estimator;
  switch (options_.global_rotation_estimator_type) {
    case GlobalRotationEstimatorType::ROBUST_L1L2: {
      // Initialize the orientation estimations by walking along the maximum
      // spanning tree.
      CHECK(OrientationsFromMaximumSpanningTree(*view_graph_, &orientations_))
          << "Could not estimate orientations from a spanning tree.";
      RobustRotationEstimator::Options robust_rotation_estimator_options;
      rotation_estimator.reset(
          new RobustRotationEstimator(robust_rotation_estimator_options));
      break;
    }
    case GlobalRotationEstimatorType::NONLINEAR: {
      // Initialize the orientation estimations by walking along the maximum
      // spanning tree.
      CHECK(OrientationsFromMaximumSpanningTree(*view_graph_, &orientations_))
          << "Could not estimate orientations from a spanning tree.";
      rotation_estimator.reset(new NonlinearRotationEstimator());
      break;
    }
    case GlobalRotationEstimatorType::LINEAR: {
      // Set the constructor variable to true to weigh each term by the inlier
      // count.
      rotation_estimator.reset(new LinearRotationEstimator());
      break;
    }
    default: {
      LOG(FATAL) << "Invalid type of global rotation estimation chosen.";
      break;
    }
  }

  // Return false if the rotation estimation does not succeed.
  if (!rotation_estimator->EstimateRotations(view_pairs, &orientations_)) {
    return false;
  }

  // Set the camera orientations of all views that were successfully estimated.
  for (const auto& orientation : orientations_) {
    View* view = reconstruction_->MutableView(orientation.first);
    // Do not update the camera orientation of this view if the pose was already
    // estimated.
    if (view == nullptr || view->IsEstimated()) {
      continue;
    }
    Camera* camera = view->MutableCamera();
    camera->SetOrientationFromAngleAxis(orientation.second);
  }

  return orientations_.size() > 0;
}

double HybridReconstructionEstimator::ComputeMedianTriangulationAngle(
    const ViewIdPair& view_ids) {
  const std::vector<ViewId> views = {view_ids.first, view_ids.second};
  const std::vector<TrackId> common_tracks =
      FindCommonTracksInViews(*reconstruction_, views);

  // Fetch the view and camera objects for this pair.
  const View* view1 = reconstruction_->View(view_ids.first);
  const View* view2 = reconstruction_->View(view_ids.second);
  const Camera& camera1 = view1->Camera();
  const Camera& camera2 = view2->Camera();

  // Compute the angle between the viewing rays.
  const Eigen::Vector2d camera1_principal_point(camera1.PrincipalPointX(),
                                                camera1.PrincipalPointY());
  const Eigen::Vector3d view1_ray =
      camera1.PixelToUnitDepthRay(camera1_principal_point).normalized();
  const Eigen::Vector2d camera2_principal_point(camera2.PrincipalPointX(),
                                                camera2.PrincipalPointY());
  const Eigen::Vector3d view2_ray =
      camera2.PixelToUnitDepthRay(camera2_principal_point).normalized();
  return std::abs(RadToDeg(std::acos(view1_ray.dot(view2_ray))));

  // NOTE: Originally, csweeney used the median triangulation angle between
  // matched features but this is fairly slow so instead I chose to use the
  // angle between principal viewing rays as the angular constraint for the view
  // pair. The code for using the median triangulation angle is below.
  //
  // std::vector<double> triangulation_angles;
  // triangulation_angles.reserve(common_tracks.size());
  // // Compute the triangulation angle for each view. We store the cosine of the
  // // dot product to save computation.
  // for (const TrackId track_id : common_tracks) {
  //   const Feature* feature1 = view1->GetFeature(track_id);
  //   const Eigen::Vector3d ray1 =
  //       camera1.PixelToUnitDepthRay(*feature1).normalized();
  //   const Feature* feature2 = view2->GetFeature(track_id);
  //   const Eigen::Vector3d ray2 =
  //       camera2.PixelToUnitDepthRay(*feature2).normalized();

  //   // Store the dot product of the rays which is the cosine of triangulation
  //   // angle.
  //   triangulation_angles.emplace_back(ray1.dot(ray2));
  // }

  // // Return the median triangulation angle.
  // const int median_index = triangulation_angles.size() / 2;
  // std::nth_element(triangulation_angles.begin(),
  //                  triangulation_angles.begin() + median_index,
  //                  triangulation_angles.end());
  // return RadToDeg(std::acos(triangulation_angles[median_index]));
}

bool HybridReconstructionEstimator::InitializeCamerasFromTwoViewInfo(
    const ViewIdPair& view_ids) {
  if (!ContainsKey(orientations_, view_ids.first)) {
    return false;
  }

  View* view1 = reconstruction_->MutableView(view_ids.first);
  View* view2 = reconstruction_->MutableView(view_ids.second);
  const TwoViewInfo* info =
      view_graph_->GetEdge(view_ids.first, view_ids.second);

  // The pose of camera 1 will be the identity pose, so we only need to set
  // camera 2's pose.
  Camera* camera1 = view1->MutableCamera();
  camera1->SetFocalLength(info->focal_length_1);
  Camera* camera2 = view2->MutableCamera();
  // The two view info position is in camera 1's coordinate system. We need to
  // transform this to the world coordinate system defined by the orientation.
  const Eigen::Matrix3d camera1_rotation =
      camera1->GetOrientationAsRotationMatrix();
  camera2->SetPosition(camera1_rotation.transpose() * info->position_2);
  camera2->SetFocalLength(info->focal_length_2);

  view1->SetEstimated(true);
  view2->SetEstimated(true);
}

bool HybridReconstructionEstimator::InitializeCamerasWithKnownOrientation(
    const ViewIdPair& view_ids) {
  if (!ContainsKey(orientations_, view_ids.first) &&
      !ContainsKey(orientations_, view_ids.second)) {
    return false;
  }

  View* view1 = reconstruction_->MutableView(view_ids.first);
  const Camera& camera1 = view1->Camera();
  View* view2 = reconstruction_->MutableView(view_ids.second);
  Camera* camera2 = view2->MutableCamera();
  const TwoViewInfo* info =
      view_graph_->GetEdge(view_ids.first, view_ids.second);

  // Get the rotated feature correspondences. These will be used to estimate an
  // improved relative position between the cameras given the known
  // orientations.
  const std::vector<ViewId> view_pair = {view_ids.first, view_ids.second};
  const std::vector<TrackId> common_tracks =
      FindCommonTracksInViews(*reconstruction_, view_pair);
  std::vector<FeatureCorrespondence> rotated_correspondences;
  rotated_correspondences.reserve(common_tracks.size());
  for (const TrackId track_id : common_tracks) {
    // Retrieve the rotated and normalized correspondences.
    FeatureCorrespondence match;
    match.feature1 =
        camera1.PixelToUnitDepthRay(*view1->GetFeature(track_id)).hnormalized();
    match.feature2 = camera2->PixelToUnitDepthRay(*view2->GetFeature(track_id))
                         .hnormalized();
    rotated_correspondences.emplace_back(match);
  }

  // Set up the ransac parameters for estimating the relative position. Compute
  // the sampson error threshold to account for the resolution of the images.
  RansacParameters relative_position_estimator_ransac_params = ransac_params_;
  const double max_sampson_error_pixels1 = ComputeResolutionScaledThreshold(
      options_.relative_position_estimation_max_sampson_error_pixels,
      camera1.ImageWidth(),
      camera1.ImageHeight());
  const double max_sampson_error_pixels2 = ComputeResolutionScaledThreshold(
      options_.relative_position_estimation_max_sampson_error_pixels,
      camera2->ImageWidth(),
      camera2->ImageHeight());
  relative_position_estimator_ransac_params.error_thresh =
      max_sampson_error_pixels1 * max_sampson_error_pixels2 /
      (info->focal_length_1 * info->focal_length_2);

  // Estimate the relative position given the known orientation. This should
  // improve the relative position estimation and give a better initial
  // configuration for the reconstruction.
  Eigen::Vector3d relative_camera2_position;
  RansacType ransac_type = RansacType::RANSAC;
  RansacSummary ransac_summary;
  if (!EstimateRelativePoseWithKnownOrientation(
          relative_position_estimator_ransac_params,
          ransac_type,
          rotated_correspondences,
          &relative_camera2_position,
          &ransac_summary)){
    return false;
  }

  // Set the position to the new position estimation.
  camera2->SetPosition(relative_camera2_position);
  view1->SetEstimated(true);
  view2->SetEstimated(true);

  return ransac_summary.inliers.size() > options_.min_num_two_view_inliers;
}

bool HybridReconstructionEstimator::ChooseInitialViewPair() {
  static const int kMinNumInitialTracks = 100;

  // Sort the view pairs by the number of geometrically verified matches.
  std::vector<ViewIdPair> candidate_initial_view_pairs;
  OrderViewPairsByInitializationCriterion(kMinNumInitialTracks,
                                          &candidate_initial_view_pairs);

  if (candidate_initial_view_pairs.size() == 0) {
    return false;
  }

  // Try to initialize the reconstruction from the candidate view pairs. An
  // initial seed is only considered valid if the baseline relative to the 3D
  // point depths is sufficient. This robustness is measured by the angle of all
  // 3D points.
  for (const ViewIdPair view_id_pair : candidate_initial_view_pairs) {
    // Set all values as unestimated and try to use the next candidate pair.
    SetReconstructionAsUnestimated(reconstruction_);

    // Initialize the camera poses of the initial views and set the two views to
    // estimated. First try to use the relative pose solver that utilizes the
    // known orientation to improve the relative translation estimate. If that
    // fails, use the two view info estimated during matching.
    if (!InitializeCamerasWithKnownOrientation(view_id_pair) &&
        !InitializeCamerasFromTwoViewInfo(view_id_pair)) {
      continue;
    }

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
      unlocalized_views_.erase(view_id_pair.first);
      unlocalized_views_.erase(view_id_pair.second);

      return true;
    }
  }

  return false;
}

void HybridReconstructionEstimator::
    OrderViewPairsByInitializationCriterion(
        const int min_num_verified_matches,
        std::vector<ViewIdPair>* view_id_pairs) {
  static const double kMaxTriangulationAngleDegrees = 45;

  const auto& view_pairs = view_graph_->GetAllEdges();
  view_id_pairs->reserve(view_pairs.size());

  // Choose the initialization criterion based on the estimated triangulation
  // angle between the views and the number of features matched between the view
  // pairs.
  std::vector<std::tuple<int, double, ViewIdPair> >
      initialization_criterion_for_view_pairs;
  initialization_criterion_for_view_pairs.reserve(view_pairs.size());
  for (const auto& view_pair : view_pairs) {
    // Estimate the triangulation angle if the orientations are known. If either
    // one of the orientations are not known then set the angle to zero. This
    // will push this view pair to the back of the list of camera
    // initialization, allowing it to still be used if all pairs with known
    // orientation fail.
    double median_triangulation_angle = 0;
    if (ContainsKey(orientations_, view_pair.first.first) &&
        ContainsKey(orientations_, view_pair.first.second)) {
      median_triangulation_angle =
          ComputeMedianTriangulationAngle(view_pair.first);
    }

    // Take the scaled sqrt of the triangulation angle and round it to the
    // nearest integer. This essentially buckets the triangulation angles in a
    // geometric sequence. We additionally cap the triangulation angle at 45
    // degrees to prevent favoring view pairs with too wide of a baseline that
    // do not have sufficiently many matched features. This allows us to choose
    // the view pair with the most number of features that have a sufficiently
    // large triangulation angle.
    const int normalized_triangulation_angle =
        std::round(2.0 * std::sqrt(std::min(median_triangulation_angle,
                                            kMaxTriangulationAngleDegrees)));

    // TODO(csweeney): Prefer view pairs with known intrinsics.
    if (view_pair.second.num_verified_matches > min_num_verified_matches) {
      // Insert negative values so that the largest triangulation angles and
      // highest number of matches appear at the front.
      initialization_criterion_for_view_pairs.emplace_back(
          -normalized_triangulation_angle,
          -view_pair.second.num_verified_matches,
          view_pair.first);
    }
  }

  // Sort the views such that the view pairs with the largest triangulation
  // angles and the most matched features appear at the front.
  std::sort(initialization_criterion_for_view_pairs.begin(),
            initialization_criterion_for_view_pairs.end());
  for (int i = 0; i < initialization_criterion_for_view_pairs.size(); i++) {
    view_id_pairs->emplace_back(
        std::get<2>(initialization_criterion_for_view_pairs[i]));
  }
}

void HybridReconstructionEstimator::FindViewsToLocalize(
    std::vector<ViewId>* views_to_localize) {
  // We localize all views that observe 75% or more than the number of 3D points
  // observed by the view with the largest number of observed 3D points.
  static const double kObserved3dPointsRatio = 0.75;

  // Determine the number of estimated tracks that each view observes.
  std::vector<std::pair<int, ViewId> > track_count_for_view;
  track_count_for_view.reserve(unlocalized_views_.size());
  for (const ViewId view_id : unlocalized_views_) {
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

void HybridReconstructionEstimator::EstimateStructure(
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

double HybridReconstructionEstimator::UnoptimizedGrowthPercentage() {
  return 100.0 * (reconstructed_views_.size() - num_optimized_views_) /
         static_cast<double>(num_optimized_views_);
}

bool HybridReconstructionEstimator::FullBundleAdjustment() {
  // Full bundle adjustment.
  LOG(INFO) << "Running full bundle adjustment on the entire reconstruction.";

  // Set up the BA options.
  bundle_adjustment_options_ =
      SetBundleAdjustmentOptions(options_, reconstructed_views_.size());

  // Do not optimize the camera orientations since they are considered to be
  // optimized.
  bundle_adjustment_options_.constant_camera_orientation = true;

  // Inner iterations are not really needed for incremental SfM because we are
  // *hopefully* already starting at a good local minima. Inner iterations are
  // disabled because they slow down BA a lot.
  bundle_adjustment_options_.use_inner_iterations = false;

  // If desired, select good tracks to optimize for BA. This dramatically
  // reduces the number of parameters in bundle adjustment, and does a decent
  // job of filtering tracks with outliers that may slow down the nonlinear
  // optimization.
  std::unordered_set<TrackId> tracks_to_optimize;
  if (options_.subsample_tracks_for_bundle_adjustment &&
      SelectGoodTracksForBundleAdjustment(
          *reconstruction_,
          options_.track_subset_selection_long_track_length_threshold,
          options_.track_selection_image_grid_cell_size_pixels,
          options_.min_num_optimized_tracks_per_view,
          &tracks_to_optimize)) {
    SetTracksInViewsToUnestimated(reconstructed_views_,
                                  tracks_to_optimize,
                                  reconstruction_);
  } else {
    GetEstimatedTracksFromReconstruction(*reconstruction_, &tracks_to_optimize);
  }
  LOG(INFO) << "Selected " << tracks_to_optimize.size()
            << " tracks to optimize.";

  std::unordered_set<ViewId> views_to_optimize;
  GetEstimatedViewsFromReconstruction(*reconstruction_,
                                      &views_to_optimize);
  const auto& ba_summary =
      BundleAdjustPartialReconstruction(bundle_adjustment_options_,
                                        views_to_optimize,
                                        tracks_to_optimize,
                                        reconstruction_);
  num_optimized_views_ = reconstructed_views_.size();

  const auto& track_ids = reconstruction_->TrackIds();
  const std::unordered_set<TrackId> all_tracks(track_ids.begin(),
                                               track_ids.end());
  RemoveOutlierTracks(all_tracks, options_.max_reprojection_error_in_pixels);

  return ba_summary.success;
}

bool HybridReconstructionEstimator::PartialBundleAdjustment() {
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

  // Do not optimize the camera orientations since they are considered to be
  // optimized.
  bundle_adjustment_options_.constant_camera_orientation = true;

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

  // If desired, select good tracks to optimize for BA. This dramatically
  // reduces the number of parameters in bundle adjustment, and does a decent
  // job of filtering tracks with outliers that may slow down the nonlinear
  // optimization.
  std::unordered_set<TrackId> tracks_to_optimize;
  if (options_.subsample_tracks_for_bundle_adjustment &&
      SelectGoodTracksForBundleAdjustment(
          *reconstruction_,
          views_to_optimize,
          options_.track_subset_selection_long_track_length_threshold,
          options_.track_selection_image_grid_cell_size_pixels,
          options_.min_num_optimized_tracks_per_view,
          &tracks_to_optimize)) {
    SetTracksInViewsToUnestimated(views_to_optimize,
                                  tracks_to_optimize,
                                  reconstruction_);
  } else {
    // If the track selection fails or is not desired, then add all tracks from
    // the views we wish to optimize.
    for (const ViewId view_to_optimize : views_to_optimize) {
      const View* view = reconstruction_->View(view_to_optimize);
      const auto& tracks_in_view = view->TrackIds();
      for (const TrackId track_in_view : tracks_in_view) {
        tracks_to_optimize.insert(track_in_view);
      }
    }
  }
  LOG(INFO) << "Selected " << tracks_to_optimize.size()
            << " tracks to optimize.";

  // TODO: We should write an interface to allow for track subset selection
  // during partial BA. This would provide obvious speedups, however, it may
  // actually be beneficial to optimize the entire local reconstruction and only
  // perform track subset selection for full BA. More testing should be done.
  ba_summary = BundleAdjustPartialReconstruction(bundle_adjustment_options_,
                                                 views_to_optimize,
                                                 tracks_to_optimize,
                                                 reconstruction_);

  RemoveOutlierTracks(tracks_to_optimize,
                      options_.max_reprojection_error_in_pixels);
  return ba_summary.success;
}

void HybridReconstructionEstimator::RemoveOutlierTracks(
    const std::unordered_set<TrackId>& tracks_to_check,
    const double max_reprojection_error_in_pixels) {
  // Remove the outlier points based on the reprojection error and how
  // well-constrained the 3D points are.
  int num_points_removed = SetOutlierTracksToUnestimated(
      tracks_to_check,
      max_reprojection_error_in_pixels,
      options_.min_triangulation_angle_degrees,
      reconstruction_);
  LOG(INFO) << num_points_removed << " outlier points were removed.";
}

void HybridReconstructionEstimator::SetUnderconstrainedAsUnestimated() {
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
          !ContainsKey(unlocalized_views_, view_id)) {
        unlocalized_views_.insert(view_id);

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
