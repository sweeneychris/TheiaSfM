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

#include "theia/sfm/reconstruction_estimator_utils.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/bundle_adjustment/optimize_relative_position_with_known_rotation.h"
#include "theia/sfm/camera/reprojection_error.h"
#include "theia/sfm/global_pose_estimation/nonlinear_position_estimator.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"

namespace theia {
namespace {

// Accumulate all two view feature matches between the input views. The features
// are normalized according to the camera intrinsics.
void GetNormalizedFeatureCorrespondences(
    const View& view1, const View& view2,
    std::vector<FeatureCorrespondence>* matches) {
  const Camera& camera1 = view1.Camera();
  const Camera& camera2 = view2.Camera();
  const std::vector<TrackId>& tracks = view1.TrackIds();
  for (const TrackId track_id : tracks) {
    const Feature* feature2 = view2.GetFeature(track_id);
    // If view 2 does not contain the current track then it cannot be a
    // correspondence.
    if (feature2 == nullptr) {
      continue;
    }

    FeatureCorrespondence match;
    const Feature* feature1 = view1.GetFeature(track_id);
    match.feature1 = *feature1;
    match.feature2 = *feature2;

    // Normalize for camera intrinsics.
    match.feature1 =
        camera1.PixelToNormalizedCoordinates(match.feature1).hnormalized();
    match.feature2 =
        camera2.PixelToNormalizedCoordinates(match.feature2).hnormalized();
    matches->emplace_back(match);
  }
}

}  // namespace

double ComputeResolutionScaledThreshold(const double threshold_pixels,
                                        const int image_width,
                                        const int image_height) {
  static const double kDefaultImageDimension = 1024.0;

  if (image_width == 0 && image_height == 0) {
    return threshold_pixels;
  }

  const int max_image_dimension = std::max(image_width, image_height);
  return threshold_pixels * static_cast<double>(max_image_dimension) /
         kDefaultImageDimension;
}

// Sets the bundle adjustment options from the reconstruction estimator options.
BundleAdjustmentOptions SetBundleAdjustmentOptions(
    const ReconstructionEstimatorOptions& options, const int num_views) {
  static const int kMinViewsForSparseSchur = 150;

  BundleAdjustmentOptions ba_options;
  ba_options.num_threads = options.num_threads;
  ba_options.loss_function_type = options.bundle_adjustment_loss_function_type;
  ba_options.robust_loss_width = options.bundle_adjustment_robust_loss_width;
  ba_options.use_inner_iterations = true;
  ba_options.intrinsics_to_optimize = options.intrinsics_to_optimize;

  if (num_views >= options.min_cameras_for_iterative_solver) {
    ba_options.linear_solver_type = ceres::ITERATIVE_SCHUR;
    // NOTE: this is an arbitrary scaling that was found to work well. It may
    // need to change depending on the application.
    ba_options.max_num_iterations *= 1.5;
  } else if (num_views >= kMinViewsForSparseSchur) {
    ba_options.linear_solver_type = ceres::SPARSE_SCHUR;
  } else {
    ba_options.linear_solver_type = ceres::DENSE_SCHUR;
  }
  ba_options.verbose = true;
  return ba_options;
}

// Sets the ransac parameters from the reconstruction estimator options.
RansacParameters SetRansacParameters(
    const ReconstructionEstimatorOptions& options) {
  RansacParameters ransac_params;
  ransac_params.rng = options.rng;
  ransac_params.failure_probability = 1.0 - options.ransac_confidence;
  ransac_params.min_iterations = options.ransac_min_iterations;
  ransac_params.max_iterations = options.ransac_max_iterations;
  ransac_params.use_mle = options.ransac_use_mle;
  return ransac_params;
}

// Collects the relative rotations for each view pair into a simple map.
std::unordered_map<ViewIdPair, Eigen::Vector3d> RelativeRotationsFromViewGraph(
    const ViewGraph& view_graph) {
  std::unordered_map<ViewIdPair, Eigen::Vector3d> relative_rotations;

  const auto& view_pairs = view_graph.GetAllEdges();
  for (const auto& view_pair : view_pairs) {
    relative_rotations[view_pair.first] = view_pair.second.rotation_2;
  }
  return relative_rotations;
}

// Each view that has a rotation and position estimated has the Camera pose set.
void SetReconstructionFromEstimatedPoses(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const std::unordered_map<ViewId, Eigen::Vector3d>& positions,
    Reconstruction* reconstruction) {
  for (const auto& position : positions) {
    View* view = reconstruction->MutableView(position.first);
    if (view == nullptr) {
      LOG(WARNING) << "Trying to set the pose of View " << position.first
                   << " which does not exist in the reconstruction.";
      continue;
    }
    CHECK(!view->IsEstimated()) << "Cannot set the pose of a view that has "
                                   "already been estimated. View Id "
                                << position.first;

    const Eigen::Vector3d* orientation =
        FindOrNull(orientations, position.first);
    if (orientation == nullptr) {
      LOG(WARNING) << "Cannot add View " << position.first
                   << " to the reconstruction because it does nto contain an "
                      "orientation estimation.";
      continue;
    }

    view->MutableCamera()->SetPosition(position.second);
    view->MutableCamera()->SetOrientationFromAngleAxis(*orientation);
    view->SetEstimated(true);
  }
}

void CreateEstimatedSubreconstruction(
    const Reconstruction& input_reconstruction,
    Reconstruction* estimated_reconstruction) {
  *estimated_reconstruction = input_reconstruction;
  const auto view_ids = estimated_reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    const View* view = estimated_reconstruction->View(view_id);
    if (view == nullptr) {
      continue;
    }

    if (!view->IsEstimated()) {
      CHECK(estimated_reconstruction->RemoveView(view_id));
    }
  }

  const auto track_ids = estimated_reconstruction->TrackIds();
  for (const TrackId track_id : track_ids) {
    const Track* track = estimated_reconstruction->Track(track_id);
    if (track == nullptr) {
      continue;
    }

    if (!track->IsEstimated() || track->NumViews() < 2) {
      CHECK(estimated_reconstruction->RemoveTrack(track_id));
    }
  }
}

// Outputs the ViewId of all estimated views in the reconstruction.
void GetEstimatedViewsFromReconstruction(const Reconstruction& reconstruction,
                                         std::unordered_set<ViewId>* views) {
  CHECK_NOTNULL(views)->clear();
  for (const ViewId view_id : reconstruction.ViewIds()) {
    const View* view = reconstruction.View(view_id);
    if (view != nullptr && view->IsEstimated()) {
      views->insert(view_id);
    }
  }
}

// Outputs the TrackId of all estimated tracks in the reconstruction.
void GetEstimatedTracksFromReconstruction(const Reconstruction& reconstruction,
                                          std::unordered_set<TrackId>* tracks) {
  CHECK_NOTNULL(tracks)->clear();
  const auto& track_ids = reconstruction.TrackIds();
  for (const TrackId track_id : track_ids) {
    const Track* track = reconstruction.Track(track_id);
    if (track != nullptr && track->IsEstimated()) {
      tracks->insert(track_id);
    }
  }
}

void RefineRelativeTranslationsWithKnownRotations(
    const Reconstruction& reconstruction,
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const int num_threads,
    ViewGraph* view_graph) {
  CHECK_GE(num_threads, 1);
  const auto& view_pairs = view_graph->GetAllEdges();

  ThreadPool pool(num_threads);
  // Refine the translation estimation for each view pair.
  for (const auto& view_pair : view_pairs) {
    // Get all feature correspondences common to both views.
    std::vector<FeatureCorrespondence> matches;
    const View* view1 = reconstruction.View(view_pair.first.first);
    const View* view2 = reconstruction.View(view_pair.first.second);
    GetNormalizedFeatureCorrespondences(*view1, *view2, &matches);

    TwoViewInfo* info = view_graph->GetMutableEdge(view_pair.first.first,
                                                   view_pair.first.second);
    pool.Add(OptimizeRelativePositionWithKnownRotation,
             matches,
             FindOrDie(orientations, view_pair.first.first),
             FindOrDie(orientations, view_pair.first.second),
             &info->position_2);
  }
}

int SetUnderconstrainedTracksToUnestimated(Reconstruction* reconstruction) {
  static const int kMinNumViews = 2;
  int num_underconstrained_tracks = 0;
  // Set all underconstrained tracks to be unestimated.
  for (const TrackId track_id : reconstruction->TrackIds()) {
    Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
    if (!track->IsEstimated()) {
      continue;
    }

    // Count the number of estimated views observing it.
    int num_estimated_views = 0;
    for (const ViewId view_id : track->ViewIds()) {
      if (reconstruction->View(view_id)->IsEstimated()) {
        ++num_estimated_views;
      }
      if (num_estimated_views >= kMinNumViews) {
        break;
      }
    }

    if (num_estimated_views < kMinNumViews) {
      track->SetEstimated(false);
      ++num_underconstrained_tracks;
    }
  }

  return num_underconstrained_tracks;
}

int SetUnderconstrainedViewsToUnestimated(Reconstruction* reconstruction) {
  // Set all underconstrained views to be unestimated.
  static const int kMinNumTracks = 3;
  int num_underconstrained_views = 0;
  // Set all underconstrained views to be unestimated.
  for (const ViewId view_id : reconstruction->ViewIds()) {
    View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
    if (!view->IsEstimated()) {
      continue;
    }

    // Count the number of estimated views observing it.
    int num_estimated_tracks = 0;
    for (const TrackId track_id : view->TrackIds()) {
      if (reconstruction->Track(track_id)->IsEstimated()) {
        ++num_estimated_tracks;
      }
      if (num_estimated_tracks >= kMinNumTracks) {
        break;
      }
    }

    if (num_estimated_tracks < kMinNumTracks) {
      view->SetEstimated(false);
      ++num_underconstrained_views;
    }
  }

  return num_underconstrained_views;
}

int NumEstimatedViews(const Reconstruction& reconstruction) {
  int num_estimated_views = 0;
  for (const ViewId view_id : reconstruction.ViewIds()) {
    const View* view = reconstruction.View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    ++num_estimated_views;
  }
  return num_estimated_views;
}

int NumEstimatedTracks(const Reconstruction& reconstruction) {
  int num_estimated_tracks = 0;
  for (const TrackId track_id : reconstruction.TrackIds()) {
    const Track* track = reconstruction.Track(track_id);
    if (track == nullptr || !track->IsEstimated()) {
      continue;
    }
    ++num_estimated_tracks;
  }
  return num_estimated_tracks;
}

}  // namespace theia
