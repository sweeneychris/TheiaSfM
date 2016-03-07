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
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/estimators/estimate_calibrated_absolute_pose.h"
#include "theia/sfm/estimators/estimate_uncalibrated_absolute_pose.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/solvers/sample_consensus_estimator.h"

namespace theia {
namespace {

// Normalize a feature by its camera intrinsics.
void NormalizeFeature(const Camera& camera,
                      const double focal_length,
                      Feature* feature) {
  feature->y() = (feature->y() - camera.PrincipalPointY()) /
                 (focal_length * camera.AspectRatio());
  feature->x() =
      (feature->x() - camera.Skew() * feature->y() - camera.PrincipalPointX()) /
      focal_length;
}

void GetNormalized2D3DCorrespondencesForView(
        const View& view,
        const Reconstruction& reconstruction,
        const double focal_length,
        std::vector<FeatureCorrespondence2D3D>* correspondences) {
  const auto& tracks_in_view = view.TrackIds();
  const Camera& camera = view.Camera();
  correspondences->reserve(tracks_in_view.size());

  for (const TrackId track_id : tracks_in_view) {
    const Track* track = reconstruction.Track(track_id);
    // We only use 3D points that have been estimated.
    if (!track->IsEstimated()) {
      continue;
    }

    FeatureCorrespondence2D3D correspondence;
    correspondence.feature = *view.GetFeature(track_id);
    // Normalize the feature by the intrinsics.
    NormalizeFeature(camera, focal_length, &correspondence.feature);

    correspondence.world_point = track->Point().hnormalized();
    correspondences->emplace_back(correspondence);
  }
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
  Camera* camera = view->MutableCamera();

  // Normalizing the pixels to remove the camera intrinsics requires dividing by
  // the focal length. In the case where the focal length is unknown we simply
  // skip this division by making the focal length equal to 1.
  static const double kUnknownFocalLength = 1.0;
  const bool known_focal_length =
      view->CameraIntrinsicsPrior().focal_length.is_set;
  const double focal_length =
      known_focal_length ? camera->FocalLength() : kUnknownFocalLength;

  // Gather all 2D-3D correspondences.
  std::vector<FeatureCorrespondence2D3D> matches;
  GetNormalized2D3DCorrespondencesForView(*view,
                                          *reconstruction,
                                          focal_length,
                                          &matches);
  if (matches.size() < options.min_num_inliers) {
    VLOG(2) << "Not enough 2D-3D correspondences to localize view "
            << view_to_localize;
    return false;
  }

  // Set up the ransac parameters for absolute pose estimation.
  bool success = false;
  // NOTE: The value of focal_length depends on if the focal length is known or
  // not so this threshold will scale appropriately.
  RansacParameters ransac_parameters = options.ransac_params;
  ransac_parameters.error_thresh = options.reprojection_error_threshold_pixels *
                                   options.reprojection_error_threshold_pixels /
                                   (focal_length * focal_length);

  // If calibrated, estimate the pose with P3P.
  if (known_focal_length) {
    CalibratedAbsolutePose pose;
    success = EstimateCalibratedAbsolutePose(
        ransac_parameters, RansacType::RANSAC, matches, &pose, summary);
    camera->SetOrientationFromRotationMatrix(pose.rotation);
    camera->SetPosition(pose.position);
  } else {
    // If the focal length is not known, estimate the focal length and pose
    // together.
    UncalibratedAbsolutePose pose;
    success = EstimateUncalibratedAbsolutePose(
        ransac_parameters, RansacType::RANSAC, matches, &pose, summary);
    camera->SetOrientationFromRotationMatrix(pose.rotation);
    camera->SetPosition(pose.position);
    camera->SetFocalLength(pose.focal_length);
  }

  // If localization failed or did not produce a sufficient number of inliers
  // then return false.
  if (!success || summary->inliers.size() < options.min_num_inliers) {
    VLOG(2) << "Failed to localize view id " << view_to_localize
            << " with only " << summary->inliers.size() << " out of "
            << matches.size() << " features as inliers.";
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
          << matches.size() << " 2D-3D matches.";
  return success;
}

}  // namespace theia
