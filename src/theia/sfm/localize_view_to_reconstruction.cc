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

#include "theia/sfm/localize_view_to_reconstruction.h"

#include <glog/logging.h>
#include <vector>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/estimators/estimate_absolute_pose_with_known_orientation.h"
#include "theia/sfm/estimators/estimate_calibrated_absolute_pose.h"
#include "theia/sfm/estimators/estimate_uncalibrated_absolute_pose.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/types.h"
#include "theia/solvers/sample_consensus_estimator.h"

namespace theia {
namespace {

bool DoesViewHaveKnownIntrinsics(const Reconstruction& reconstruction,
                                 const ViewId view_id) {
  const View* view = reconstruction.View(view_id);
  // Return true if the EXIF data provides a focal length guess
  if (view->CameraIntrinsicsPrior().focal_length.is_set) {
    return true;
  }

  // If the camera has shared intrinsics, return true if one of the shared
  // cameras has been estimated.
  const CameraIntrinsicsGroupId intrinsics_group_id =
      reconstruction.CameraIntrinsicsGroupIdFromViewId(view_id);
  const std::unordered_set<ViewId> views_in_intrinsics_group =
      reconstruction.GetViewsInCameraIntrinsicGroup(intrinsics_group_id);
  for (const ViewId shared_intrinsics_view_id : views_in_intrinsics_group) {
    const View* view = reconstruction.View(view_id);
    if (view->IsEstimated() && view_id != shared_intrinsics_view_id) {
      return true;
    }
  }
  return false;
}

void GetNormalized2D3DMatches(const Reconstruction& reconstruction,
                              const View& view,
                              std::vector<FeatureCorrespondence2D3D>* matches) {
  const Camera& camera = view.Camera();
  const auto& tracks_in_view = view.TrackIds();
  matches->reserve(tracks_in_view.size());
  for (const TrackId track_id : tracks_in_view) {
    const Track* track = reconstruction.Track(track_id);
    // We only use 3D points that have been estimated.
    if (!track->IsEstimated()) {
      continue;
    }

    FeatureCorrespondence2D3D correspondence;
    const Feature& feature = *view.GetFeature(track_id);
    // Simply shift the pixel to remove the effect of the principal point.
    correspondence.feature =
        feature -
        Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY());
    correspondence.world_point = track->Point().hnormalized();
    matches->emplace_back(correspondence);
  }
}

void GetIntrinsicsNormalized2D3DMatches(
    const Reconstruction& reconstruction,
    const View& view,
    std::vector<FeatureCorrespondence2D3D>* matches) {
  const Camera& camera = view.Camera();
  const auto& tracks_in_view = view.TrackIds();
  matches->reserve(tracks_in_view.size());
  for (const TrackId track_id : tracks_in_view) {
    const Track* track = reconstruction.Track(track_id);
    // We only use 3D points that have been estimated.
    if (!track->IsEstimated()) {
      continue;
    }

    FeatureCorrespondence2D3D correspondence;
    const Feature& feature = *view.GetFeature(track_id);
    // Remove the camera intrinsics from the feature.
    correspondence.feature =
        camera.PixelToNormalizedCoordinates(feature).hnormalized();
    correspondence.world_point = track->Point().hnormalized();
    matches->emplace_back(correspondence);
  }
}

bool EstimateCameraPose(const bool known_intrinsics,
                        const LocalizeViewToReconstructionOptions& options,
                        const Reconstruction& reconstruction,
                        View* view,
                        RansacSummary* summary) {
  Camera* camera = view->MutableCamera();

  // Gather all 2D-3D correspondences.
  std::vector<FeatureCorrespondence2D3D> matches;
  if (known_intrinsics) {
    GetIntrinsicsNormalized2D3DMatches(reconstruction, *view, &matches);
  } else {
    GetNormalized2D3DMatches(reconstruction, *view, &matches);
  }

  // Exit early if there are not enough putative matches.
  if (matches.size() < options.min_num_inliers) {
    VLOG(2) << "Not enough 2D-3D correspondences to localize view "
            << view->Name();
    return false;
  }

  // Set up the ransac parameters for absolute pose estimation.
  RansacParameters ransac_parameters = options.ransac_params;

  // Compute the reprojection error threshold scaled to account for the image
  // resolution.
  const double resolution_scaled_reprojection_error_threshold_pixels =
      ComputeResolutionScaledThreshold(
          options.reprojection_error_threshold_pixels,
          camera->ImageWidth(),
          camera->ImageHeight());

  // If we are assuming that the orientation is known then first try to use
  // the simplified camera positions solver.
  if (options.assume_known_orientation) {
    ransac_parameters.error_thresh =
        resolution_scaled_reprojection_error_threshold_pixels *
        resolution_scaled_reprojection_error_threshold_pixels /
        (camera->FocalLength() * camera->FocalLength());

    // Return true if position estimation is successful. Otherwise the method
    // will proceed to estimate the full pose.
    const Eigen::Vector3d camera_orientation =
        camera->GetOrientationAsAngleAxis();
    Eigen::Vector3d camera_position;
    if (EstimateAbsolutePoseWithKnownOrientation(ransac_parameters,
                                                 RansacType::RANSAC,
                                                 camera_orientation,
                                                 matches,
                                                 &camera_position,
                                                 summary) &&
        summary->inliers.size() > options.min_num_inliers) {
      camera->SetPosition(camera_position);
      return true;
    } else {
      return false;
    }
  }

  // If calibrated, estimate the pose with P3P.
  bool success = false;
  if (known_intrinsics) {
    ransac_parameters.error_thresh =
        resolution_scaled_reprojection_error_threshold_pixels *
        resolution_scaled_reprojection_error_threshold_pixels /
        (camera->FocalLength() * camera->FocalLength());
    CalibratedAbsolutePose pose;
    if (EstimateCalibratedAbsolutePose(
            ransac_parameters, RansacType::RANSAC, matches, &pose, summary)) {
      camera->SetOrientationFromRotationMatrix(pose.rotation);
      camera->SetPosition(pose.position);
      return true;
    }
  } else {
    // If the focal length is not known, estimate the focal length and pose
    // together.
    ransac_parameters.error_thresh =
        resolution_scaled_reprojection_error_threshold_pixels *
        resolution_scaled_reprojection_error_threshold_pixels;

    UncalibratedAbsolutePose pose;
    if (EstimateUncalibratedAbsolutePose(
            ransac_parameters, RansacType::RANSAC, matches, &pose, summary)) {
      camera->SetOrientationFromRotationMatrix(pose.rotation);
      camera->SetPosition(pose.position);
      camera->SetFocalLength(pose.focal_length);
      return true;
    }
  }

  return false;
}

}  // namespace

bool LocalizeViewToReconstruction(
    const ViewId view_to_localize,
    const LocalizeViewToReconstructionOptions options,
    Reconstruction* reconstruction,
    RansacSummary* summary) {
  CHECK_NOTNULL(reconstruction);
  CHECK_NOTNULL(summary);

  View* view = reconstruction->MutableView(view_to_localize);
  // We assume that the intrinsics are known if the orientation is known.
  const bool known_intrinsics =
      options.assume_known_orientation ||
      DoesViewHaveKnownIntrinsics(*reconstruction, view_to_localize);

  // If localization failed or did not produce a sufficient number of inliers
  // then return false.
  bool success = EstimateCameraPose(
      known_intrinsics, options, *reconstruction, view, summary);
  if (!success || summary->inliers.size() < options.min_num_inliers) {
    VLOG(2) << "Failed to localize view id " << view_to_localize
            << " with only " << summary->inliers.size() << " out of "
            << summary->num_input_data_points << " features as inliers.";
    return false;
  }

  // Bundle adjust the view if desired.
  view->SetEstimated(true);
  if (options.bundle_adjust_view) {
    const BundleAdjustmentSummary summary =
        BundleAdjustView(options.ba_options, view_to_localize, reconstruction);
    success = summary.success;
  }

  VLOG(2) << "Estimated the camera pose for view " << view_to_localize
          << " with " << summary->inliers.size() << " inliers out of "
          << summary->num_input_data_points << " 2D-3D matches.";
  return success;
}

}  // namespace theia
