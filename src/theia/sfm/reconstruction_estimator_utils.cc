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
#include "theia/sfm/pose/estimate_positions_nonlinear.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"

namespace theia {
namespace {

// Accumulate all two view feature matches between the input views. The features
// are normalized according to the camera intrinsics.
void GetFeatureCorrespondences(const View& view1, const View& view2,
                               std::vector<FeatureCorrespondence>* matches) {
  Eigen::Matrix3d calibration1, calibration2;
  view1.Camera().GetCalibrationMatrix(&calibration1);
  view2.Camera().GetCalibrationMatrix(&calibration2);
  Eigen::Matrix3d inv_calibration1, inv_calibration2;
  bool view1_invertible, view2_invertible;
  double determinant;
  calibration1.computeInverseAndDetWithCheck(inv_calibration1, determinant,
                                             view1_invertible);
  calibration2.computeInverseAndDetWithCheck(inv_calibration2, determinant,
                                             view2_invertible);
  if (!view1_invertible || !view2_invertible) {
    LOG(FATAL) << "Calibration matrices are ill formed. Cannot optimize "
                  "epipolar constraints.";
    return;
  }

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
        (inv_calibration1 * match.feature1.homogeneous()).eval().hnormalized();
    match.feature2 =
        (inv_calibration2 * match.feature2.homogeneous()).eval().hnormalized();
    matches->emplace_back(match);
  }
}

// Computes the relative rotation from two absolute rotations.
Eigen::Vector3d RelativeRotationFromTwoRotations(
    const Eigen::Vector3d& rotation1, const Eigen::Vector3d& rotation2) {
  Eigen::Matrix3d rotation1_mat, rotation2_mat;
  ceres::AngleAxisToRotationMatrix(
      rotation1.data(), ceres::ColumnMajorAdapter3x3(rotation1_mat.data()));
  ceres::AngleAxisToRotationMatrix(
      rotation2.data(), ceres::ColumnMajorAdapter3x3(rotation2_mat.data()));

  const Eigen::Matrix3d relative_rotation_mat =
      rotation2_mat * rotation1_mat.transpose();
  Eigen::Vector3d relative_rotation;
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(relative_rotation_mat.data()),
      relative_rotation.data());

  return relative_rotation;
}

}  // namespace

using Eigen::Vector3d;

// Sets the bundle adjustment optiosn from the reconstruction estimator options.
BundleAdjustmentOptions SetBundleAdjustmentOptions(
    const ReconstructionEstimatorOptions& options, const int num_views) {
  static const int kMinViewsForSparseSchur = 150;

  BundleAdjustmentOptions ba_options;
  ba_options.num_threads = options.num_threads;
  ba_options.use_inner_iterations = true;
  ba_options.constant_camera_intrinsics = options.constant_camera_intrinsics;

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
  ransac_params.failure_probability = 1.0 - options.ransac_confidence;
  ransac_params.min_iterations = options.ransac_min_iterations;
  ransac_params.max_iterations = options.ransac_max_iterations;
  ransac_params.use_mle = options.ransac_use_mle;
  return ransac_params;
}

// Initializes the focal length of each view. If EXIF data is available then the
// focal length is set using the EXIF value, otherwise it is set to
// 1.2 * max(image_width, image_height). This value is shown to be a decent
// initialization in the VisualSfM software.
void InitializeFocalLengthsFromImageSize(Reconstruction* reconstruction) {
  for (const ViewId view_id : reconstruction->ViewIds()) {
    View* view = reconstruction->MutableView(view_id);

    // Set it to the EXIF value if the focal length is known, otherwise set it
    // to the estimate based on the image size.
    if (view->CameraIntrinsicsPrior().focal_length.is_set) {
      view->MutableCamera()->SetFocalLength(
          view->CameraIntrinsicsPrior().focal_length.value);
    } else {
      view->MutableCamera()->SetFocalLength(
          1.2 * static_cast<double>(std::max(view->Camera().ImageWidth(),
                                             view->Camera().ImageHeight())));
    }
  }
}

void InitializeFocalLengthsFromMedian(const ViewGraph& view_graph,
                                      Reconstruction* reconstruction) {
  CHECK_GT(view_graph.NumEdges(), 0);
  std::unordered_map<ViewId, std::vector<double> > focal_lengths;
  focal_lengths.reserve(view_graph.NumViews());

  // Collect all focal lengths.
  const auto& edges = view_graph.GetAllEdges();
  for (const auto& edge : edges) {
    focal_lengths[edge.first.first].push_back(edge.second.focal_length_1);
    focal_lengths[edge.first.second].push_back(edge.second.focal_length_2);
  }

  // For each view, find the median focal length.
  for (auto& focal_length_vals : focal_lengths) {
    const int median_index = focal_length_vals.second.size() / 2;
    std::nth_element(focal_length_vals.second.begin(),
                     focal_length_vals.second.begin() + median_index,
                     focal_length_vals.second.end());
    View* view =
        CHECK_NOTNULL(reconstruction->MutableView(focal_length_vals.first));
    view->MutableCamera()->SetFocalLength(
        focal_length_vals.second[median_index]);
  }
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

    const Vector3d* orientation = FindOrNull(orientations, position.first);
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
    GetFeatureCorrespondences(*view1, *view2, &matches);

    const Eigen::Vector3d relative_rotation = RelativeRotationFromTwoRotations(
        FindOrDie(orientations, view_pair.first.first),
        FindOrDie(orientations, view_pair.first.second));

    TwoViewInfo* info = view_graph->GetMutableEdge(view_pair.first.first,
                                                   view_pair.first.second);
    pool.Add(OptimizeRelativePositionWithKnownRotation,
             matches,
             relative_rotation,
             &info->position_2);
  }
}

int RemoveOutlierFeatures(const double max_inlier_reprojection_error,
                           Reconstruction* reconstruction) {
  const double max_sq_reprojection_error =
      max_inlier_reprojection_error * max_inlier_reprojection_error;
  int num_features_removed = 0;
  int num_features = 0;
  // Compute reprojection error and remove all outlier features.
  for (const ViewId view_id : reconstruction->ViewIds()) {
    View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
    if (!view->IsEstimated()) {
      continue;
    }
    const Camera& camera = view->Camera();

    // Check the reprojection error of every track observed by the view.
    for (const TrackId track_id : view->TrackIds()) {
      Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
      const Feature* feature =
          reconstruction->View(view_id)->GetFeature(track_id);
      if (!track->IsEstimated()) {
        continue;
      }

      ++num_features;
      // Reproject the observations.
      Eigen::Vector2d projection;
      const double depth = camera.ProjectPoint(track->Point(), &projection);
      // Remove the feature if the reprojection error is too large or is behind
      // the camera.
      if (depth < 0 ||
          (projection - *feature).squaredNorm() > max_sq_reprojection_error) {
        view->RemoveFeature(track_id);
        track->RemoveView(view_id);
        ++num_features_removed;
      }
    }
  }
  LOG(INFO) << "Num features = " << num_features;
  return num_features_removed;
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

}  // namespace theia
