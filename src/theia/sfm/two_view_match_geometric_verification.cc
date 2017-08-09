// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include "theia/sfm/two_view_match_geometric_verification.h"

#include <glog/logging.h>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/matching/guided_epipolar_matcher.h"
#include "theia/sfm/bundle_adjustment/bundle_adjust_two_views.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/estimators/estimate_homography.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/set_camera_intrinsics_from_priors.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/twoview_info.h"

namespace theia {

namespace {

void SetupCameras(const CameraIntrinsicsPrior& intrinsics1,
                  const CameraIntrinsicsPrior& intrinsics2,
                  const TwoViewInfo& info,
                  Camera* camera1,
                  Camera* camera2) {
  camera1->SetFromCameraIntrinsicsPriors(intrinsics1);
  camera1->SetFocalLength(info.focal_length_1);

  camera2->SetOrientationFromAngleAxis(info.rotation_2);
  camera2->SetPosition(info.position_2);
  camera2->SetFromCameraIntrinsicsPriors(intrinsics2);
  camera2->SetFocalLength(info.focal_length_2);
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

}  // namespace

TwoViewMatchGeometricVerification::TwoViewMatchGeometricVerification(
    const TwoViewMatchGeometricVerification::Options& options,
    const CameraIntrinsicsPrior& intrinsics1,
    const CameraIntrinsicsPrior& intrinsics2,
    const KeypointsAndDescriptors& features1,
    const KeypointsAndDescriptors& features2,
    const std::vector<IndexedFeatureMatch>& matches)
    : options_(options),
      intrinsics1_(intrinsics1),
      intrinsics2_(intrinsics2),
      features1_(features1),
      features2_(features2),
      matches_(matches) {}

void TwoViewMatchGeometricVerification::CreateCorrespondencesFromIndexedMatches(
    std::vector<FeatureCorrespondence>* correspondences) {
  CHECK_NOTNULL(correspondences)->clear();
  correspondences->reserve(matches_.size());
  for (int i = 0; i < matches_.size(); i++) {
    const Keypoint& keypoint1 = features1_.keypoints[matches_[i].feature1_ind];
    const Keypoint& keypoint2 = features2_.keypoints[matches_[i].feature2_ind];
    correspondences->emplace_back(Feature(keypoint1.x(), keypoint1.y()),
                                  Feature(keypoint2.x(), keypoint2.y()));
  }
}

bool TwoViewMatchGeometricVerification::VerifyMatches(
    std::vector<FeatureCorrespondence>* verified_matches,
    TwoViewInfo* twoview_info) {
  if (matches_.size() < options_.min_num_inlier_matches) {
    return false;
  }

  std::vector<FeatureCorrespondence> correspondences;
  CreateCorrespondencesFromIndexedMatches(&correspondences);

  // Estimate a homography (before the matches_ container is modified).
  twoview_info->num_homography_inliers = CountHomographyInliers();

  // Estimate 2-view geometry from feature matches.
  std::vector<int> inlier_indices;
  if (!EstimateTwoViewInfo(options_.estimate_twoview_info_options,
                           intrinsics1_,
                           intrinsics2_,
                           correspondences,
                           twoview_info,
                           &inlier_indices)) {
    return false;
  }
  VLOG(2) << inlier_indices.size()
          << " matches passed initial geometric verification out of "
          << matches_.size() << " putative matches.";

  if (inlier_indices.size() < options_.min_num_inlier_matches) {
    return false;
  }

  // Update the current set of matches.
  std::vector<IndexedFeatureMatch> new_matches;
  new_matches.reserve(inlier_indices.size());
  for (int i = 0; i < inlier_indices.size(); ++i) {
    new_matches.emplace_back(matches_[inlier_indices[i]]);
  }
  matches_.swap(new_matches);

  // Get Camera objects for guided matching, triangulation, and bundle
  // adjustment.
  SetupCameras(intrinsics1_, intrinsics2_, *twoview_info, &camera1_, &camera2_);

  // Perform guided matching if desired.
  if (options_.guided_matching) {
    GuidedEpipolarMatcher::Options guided_matching_options;
    guided_matching_options.guided_matching_max_distance_pixels =
        options_.guided_matching_max_distance_pixels;
    guided_matching_options.lowes_ratio = options_.guided_matching_lowes_ratio;

    GuidedEpipolarMatcher guided_matcher(
        guided_matching_options, camera1_, camera2_, features1_, features2_);
    if (!guided_matcher.GetMatches(&matches_)) {
      return false;
    }
  }

  // Perform BA if desired.
  if (options_.bundle_adjustment &&
      matches_.size() > options_.min_num_inlier_matches) {
    if (!BundleAdjustRelativePose(twoview_info)) {
      return false;
    }
  }

  // Set the number of verified matches and the output verified_matches.
  CreateCorrespondencesFromIndexedMatches(verified_matches);
  twoview_info->num_verified_matches = verified_matches->size();
  return verified_matches->size() > options_.min_num_inlier_matches;
}

// Triangulates the points and updates the matches_
void TwoViewMatchGeometricVerification::TriangulatePoints(
    std::vector<Eigen::Vector4d>* triangulated_points) {
  CHECK_NOTNULL(triangulated_points)->reserve(matches_.size());
  const double triangulation_sq_max_reprojection_error_pixels =
      options_.triangulation_max_reprojection_error *
      options_.triangulation_max_reprojection_error;

  // Triangulate all points, throwing out the ones with bad initial reprojection
  // errors.
  const std::vector<Eigen::Vector3d> origins = {camera1_.GetPosition(),
                                                camera2_.GetPosition()};
  std::vector<IndexedFeatureMatch> triangulated_matches;
  triangulated_matches.reserve(matches_.size());
  int num_bad_triangulation_angles = 0;
  int num_failed_triangulations = 0;
  int num_bad_reprojection_errors = 0;
  for (int i = 0; i < matches_.size(); i++) {
    const Keypoint& keypoint1 = features1_.keypoints[matches_[i].feature1_ind];
    const Keypoint& keypoint2 = features2_.keypoints[matches_[i].feature2_ind];
    const Feature feature1(keypoint1.x(), keypoint1.y());
    const Feature feature2(keypoint2.x(), keypoint2.y());

    // Make sure that there is enough baseline between the point so that the
    // triangulation is well-constrained.
    std::vector<Eigen::Vector3d> ray_directions(2);
    ray_directions[0] = camera1_.PixelToUnitDepthRay(feature1).normalized();
    ray_directions[1] = camera2_.PixelToUnitDepthRay(feature2).normalized();
    if (!SufficientTriangulationAngle(
            ray_directions, options_.min_triangulation_angle_degrees)) {
      ++num_bad_triangulation_angles;
      continue;
    }

    Eigen::Vector4d point3d;
    if (!TriangulateMidpoint(origins, ray_directions, &point3d)) {
      ++num_failed_triangulations;
      continue;
    }

    // Only consider triangulation a success if the initial triangulation has a
    // small enough reprojection error.
    if (!AcceptableReprojectionError(
            camera1_,
            feature1,
            point3d,
            triangulation_sq_max_reprojection_error_pixels) ||
        !AcceptableReprojectionError(
            camera2_,
            feature2,
            point3d,
            triangulation_sq_max_reprojection_error_pixels)) {
      ++num_bad_reprojection_errors;
      continue;
    }

    triangulated_points->emplace_back(point3d);
    triangulated_matches.emplace_back(matches_[i]);
  }
  VLOG(2) << "Num acceptable triangulations = " << triangulated_matches.size()
          << " out of " << matches_.size() << " total matches. "
          << num_bad_triangulation_angles
          << " matches failed due to bad triangulation angles, "
          << num_bad_reprojection_errors
          << " matches failed due to inaccurate triangulation, and "
          << num_failed_triangulations
          << " failed to reliably triangulate and initial 3d point.";

  // Set the matches to be only the triangulated matches.
  matches_.swap(triangulated_matches);
}

bool TwoViewMatchGeometricVerification::BundleAdjustRelativePose(
    TwoViewInfo* twoview_info) {
  const double final_sq_max_reprojection_error_pixels =
      options_.final_max_reprojection_error *
      options_.final_max_reprojection_error;

  // Triangulate the points. This updates the matches_ container with only the
  // points that could be accurately triangulated.
  std::vector<Eigen::Vector4d> triangulated_points;
  TriangulatePoints(&triangulated_points);

  // Exit early if there are not enough inliers left.
  if (matches_.size() < options_.min_num_inlier_matches) {
    return false;
  }

  // Bundle adjust the relative pose and points.
  TwoViewBundleAdjustmentOptions two_view_ba_options;
  two_view_ba_options.ba_options.verbose = false;
  two_view_ba_options.ba_options.linear_solver_type = ceres::DENSE_SCHUR;
  two_view_ba_options.constant_camera1_intrinsics =
      intrinsics1_.focal_length.is_set;
  two_view_ba_options.constant_camera2_intrinsics =
      intrinsics2_.focal_length.is_set;
  two_view_ba_options.ba_options.use_inner_iterations = false;

  std::vector<FeatureCorrespondence> triangulated_correspondences;
  CreateCorrespondencesFromIndexedMatches(&triangulated_correspondences);
  BundleAdjustmentSummary summary =
      BundleAdjustTwoViews(two_view_ba_options,
                           triangulated_correspondences,
                           &camera1_,
                           &camera2_,
                           &triangulated_points);

  if (!summary.success) {
    return false;
  }

  // Remove points with high reprojection errors.
  std::vector<IndexedFeatureMatch> inliers_after_ba;
  inliers_after_ba.reserve(matches_.size());
  for (int i = 0; i < triangulated_correspondences.size(); i++) {
    const auto& correspondence = triangulated_correspondences[i];
    const Eigen::Vector4d& point3d = triangulated_points[i];
    if (AcceptableReprojectionError(camera1_,
                                    correspondence.feature1,
                                    point3d,
                                    final_sq_max_reprojection_error_pixels) &&
        AcceptableReprojectionError(camera2_,
                                    correspondence.feature2,
                                    point3d,
                                    final_sq_max_reprojection_error_pixels)) {
      inliers_after_ba.emplace_back(matches_[i]);
    }
  }
  VLOG(2) << inliers_after_ba.size() << " valid matches after BA out of "
          << matches_.size() << " triangulated matches.";
  matches_.swap(inliers_after_ba);

  // Update the relative pose.
  twoview_info->rotation_2 = camera2_.GetOrientationAsAngleAxis();
  twoview_info->position_2 = camera2_.GetPosition();
  twoview_info->position_2.normalize();
  twoview_info->focal_length_1 = camera1_.FocalLength();
  twoview_info->focal_length_2 = camera2_.FocalLength();

  return true;
}

// Compute a homography and return the number of inliers. This determines how
// well a plane fits the two view geometry.
int TwoViewMatchGeometricVerification::CountHomographyInliers() {
  const EstimateTwoViewInfoOptions& etvi_options =
      options_.estimate_twoview_info_options;
  RansacParameters homography_params;
  homography_params.rng = etvi_options.rng;

  // Compute the reprojection error threshold to account for the resolution of
  // the images.
  const double max_sampson_error_pixels1 =
      ComputeResolutionScaledThreshold(etvi_options.max_sampson_error_pixels,
                                       camera1_.ImageWidth(),
                                       camera1_.ImageHeight());
  const double max_sampson_error_pixels2 =
      ComputeResolutionScaledThreshold(etvi_options.max_sampson_error_pixels,
                                       camera2_.ImageWidth(),
                                       camera2_.ImageHeight());

  homography_params.error_thresh =
      max_sampson_error_pixels1 * max_sampson_error_pixels2;
  homography_params.max_iterations = etvi_options.max_ransac_iterations;
  homography_params.min_iterations = etvi_options.min_ransac_iterations;
  homography_params.use_mle = etvi_options.use_mle;
  homography_params.failure_probability =
      1.0 - etvi_options.expected_ransac_confidence;
  RansacSummary homography_summary;
  Eigen::Matrix3d unused_homography;
  std::vector<FeatureCorrespondence> correspondences;
  CreateCorrespondencesFromIndexedMatches(&correspondences);
  EstimateHomography(homography_params,
                     etvi_options.ransac_type,
                     correspondences,
                     &unused_homography,
                     &homography_summary);

  return homography_summary.inliers.size();
}

}  // namespace theia
