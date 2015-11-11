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

#include "theia/sfm/verify_two_view_matches.h"

#include <glog/logging.h>
#include <vector>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/bundle_adjustment/bundle_adjust_two_views.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/twoview_info.h"

namespace theia {

namespace {

void SetCameraIntrinsicsFromPriors(const CameraIntrinsicsPrior& prior,
                                   Camera* camera) {
  if (prior.principal_point[0].is_set && prior.principal_point[1].is_set) {
    camera->SetPrincipalPoint(prior.principal_point[0].value,
                              prior.principal_point[1].value);
  } else {
    camera->SetPrincipalPoint(prior.image_width / 2.0,
                              prior.image_height / 2.0);
  }

  if (prior.aspect_ratio.is_set) {
    camera->SetAspectRatio(prior.aspect_ratio.value);
  }

  if (prior.skew.is_set) {
    camera->SetSkew(prior.skew.value);
  }
}

void SetupCameras(const CameraIntrinsicsPrior& intrinsics1,
                  const CameraIntrinsicsPrior& intrinsics2,
                  const TwoViewInfo& info,
                  Camera* camera1,
                  Camera* camera2) {
  camera1->SetFocalLength(info.focal_length_1);
  SetCameraIntrinsicsFromPriors(intrinsics1, camera1);

  camera2->SetOrientationFromAngleAxis(info.rotation_2);
  camera2->SetPosition(info.position_2);
  camera2->SetFocalLength(info.focal_length_2);
  SetCameraIntrinsicsFromPriors(intrinsics2, camera2);
}

// Returns false if the reprojection error of the triangulated point is greater
// than the max allowable reprojection error and true otherwise.
bool AcceptableReprojectionError(
    const Camera& camera,
    const Feature& feature,
    const Eigen::Vector4d& triangulated_point,
    const double sq_max_reprojection_error_pixels) {
  Eigen::Vector2d reprojection;
  if (camera.ProjectPoint(triangulated_point, &reprojection) < 0) {
    return false;
  }
  const double sq_reprojection_error = (feature - reprojection).squaredNorm();
  return sq_reprojection_error < sq_max_reprojection_error_pixels;
}

// Triangulate the correspondences, keeping only the triangulated points with an
// acceptable initial reprojection error.
void TriangulatePoints(
    const double sq_max_reprojection_error_pixels,
    const Camera& camera1,
    const Camera& camera2,
    const std::vector<FeatureCorrespondence>& correspondences,
    std::vector<int>* inliers,
    std::vector<Eigen::Vector4d>* tracks) {
  Matrix3x4d projection_matrix1, projection_matrix2;
  camera1.GetProjectionMatrix(&projection_matrix1);
  camera2.GetProjectionMatrix(&projection_matrix2);

  tracks->reserve(inliers->size());
  std::vector<int> triangulated_inliers;
  for (int i = 0; i < inliers->size(); i++) {
    const auto& correspondence = correspondences[inliers->at(i)];
    Eigen::Vector4d point3d;
    if (!Triangulate(projection_matrix1,
                     projection_matrix2,
                     correspondence.feature1,
                     correspondence.feature2,
                     &point3d)) {
      continue;
    }

    // Only consider triangulation a success if the initial triangulation has a
    // small enough reprojection error.
    if (AcceptableReprojectionError(camera1, correspondence.feature1, point3d,
                                    sq_max_reprojection_error_pixels) &&
        AcceptableReprojectionError(camera2, correspondence.feature2, point3d,
                                    sq_max_reprojection_error_pixels)) {
      tracks->emplace_back(point3d);
      triangulated_inliers.emplace_back(inliers->at(i));
    }
  }

  // Set the output inliers to the successfull triangualted points.
  inliers->swap(triangulated_inliers);

  VLOG(1) << "Num acceptable triangulations = " << inliers->size()
          << " out of " << triangulated_inliers.size() << " total matches.";
}

TwoViewBundleAdjustmentOptions SetTwoViewBundleAdjustmentOptions(
    const CameraIntrinsicsPrior& intrinsics1,
    const CameraIntrinsicsPrior& intrinsics2) {
  TwoViewBundleAdjustmentOptions two_view_ba_options;
  two_view_ba_options.ba_options.verbose = false;
  two_view_ba_options.ba_options.linear_solver_type = ceres::DENSE_SCHUR;
  two_view_ba_options.constant_camera1_intrinsics =
      intrinsics1.focal_length.is_set;
  two_view_ba_options.constant_camera2_intrinsics =
      intrinsics2.focal_length.is_set;
  two_view_ba_options.ba_options.use_inner_iterations = false;
  return two_view_ba_options;
}

bool BundleAdjustRelativePose(
    TwoViewBundleAdjustmentOptions& two_view_ba_options,
    const std::vector<FeatureCorrespondence>& correspondences,
    const std::vector<int>& inliers,
    Camera* camera1,
    Camera* camera2,
    std::vector<Eigen::Vector4d>* triangulated_points) {

  std::vector<FeatureCorrespondence> triangulated_features(inliers.size());
  for (int i = 0; i < inliers.size(); i++) {
    triangulated_features[i] = correspondences[inliers[i]];
  }

  BundleAdjustmentSummary summary =
      BundleAdjustTwoViews(two_view_ba_options, triangulated_features, camera1,
                           camera2, triangulated_points);

  return summary.success;
}

}  // namespace

bool VerifyTwoViewMatches(
    const VerifyTwoViewMatchesOptions& options,
    const CameraIntrinsicsPrior& intrinsics1,
    const CameraIntrinsicsPrior& intrinsics2,
    const std::vector<FeatureCorrespondence>& correspondences,
    TwoViewInfo* twoview_info,
    std::vector<int>* inlier_indices) {
  if (correspondences.size() < options.min_num_inlier_matches) {
    return false;
  }

  // Estimate the two view info. If we fail to estimate a two view info then do
  // not add this view pair to the verified matches.
  if (!EstimateTwoViewInfo(options.estimate_twoview_info_options,
                           intrinsics1,
                           intrinsics2,
                           correspondences,
                           twoview_info,
                           inlier_indices)) {
    return false;
  }

  // Bundle adjustment (optional).
  if (options.bundle_adjustment) {
    if (inlier_indices->size() < options.min_num_inlier_matches) {
      return false;
    }

    // Get Camera objects for triangulation and bundle adjustment.
    Camera camera1, camera2;
    SetupCameras(intrinsics1, intrinsics2, *twoview_info, &camera1, &camera2);

    std::vector<Eigen::Vector4d> tracks;
    TriangulatePoints(options.triangulation_sq_max_reprojection_error, camera1,
                      camera2, correspondences, inlier_indices, &tracks);

    if (inlier_indices->size() < options.min_num_inlier_matches) {
      return false;
    }

    // Bundle adjust the relative pose and points.
    TwoViewBundleAdjustmentOptions ba_options =
        SetTwoViewBundleAdjustmentOptions(intrinsics1, intrinsics2);
    if (!BundleAdjustRelativePose(ba_options, correspondences, *inlier_indices,
                                  &camera1, &camera2, &tracks)) {
      return false;
    }

    // Remove points with high reprojection errors.
    std::vector<int> inliers_after_ba;
    inliers_after_ba.reserve(inlier_indices->size());
    for (int i = 0; i < inlier_indices->size(); i++) {
      const auto& correspondence = correspondences[inlier_indices->at(i)];
      const Eigen::Vector4d& point3d = tracks[i];
      if (AcceptableReprojectionError(
              camera1, correspondence.feature1, point3d,
              options.final_sq_max_reprojection_error) &&
          AcceptableReprojectionError(
              camera2, correspondence.feature2, point3d,
              options.final_sq_max_reprojection_error)) {
        inliers_after_ba.emplace_back(inlier_indices->at(i));
      }
    }
    inlier_indices->swap(inliers_after_ba);
    VLOG(1) << inlier_indices->size() << " valid matches after BA out of "
            << inliers_after_ba.size() << " total matches.";

    // Update the relative pose.
    twoview_info->rotation_2 = camera2.GetOrientationAsAngleAxis();
    twoview_info->position_2 = camera2.GetPosition();
    twoview_info->position_2.normalize();
    twoview_info->focal_length_1 = camera1.FocalLength();
    twoview_info->focal_length_2 = camera2.FocalLength();
  }
  return inlier_indices->size() > options.min_num_inlier_matches;
}

}  // namespace theia
