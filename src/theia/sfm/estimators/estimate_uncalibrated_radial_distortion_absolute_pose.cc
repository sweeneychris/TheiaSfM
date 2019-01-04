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


// This file was created by Steffen Urban (urbste@googlemail.com) or company address (steffen.urban@zeiss.com)
// January 2019

#include "theia/sfm/estimators/estimate_uncalibrated_radial_distortion_absolute_pose.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <vector>

#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/pose/four_point_focal_length_radial_distortion.h"
#include "theia/sfm/types.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

// An estimator for computing the uncalibrated absolute pose from 4 feature
// correspondences. The feature correspondences should be normalized such that
// the principal point is at (0, 0). This can be achieved by multiplying image points by
// normalized image points = inv([1 0 image_width/2; 0 1 image_height/2; 0 0 1]) * image_points;
class UncalibratedRadialAbsolutePoseEstimator
    : public Estimator<FeatureCorrespondence2D3D, UncalibratedRadialAbsolutePose> {
 public:
  UncalibratedRadialAbsolutePoseEstimator() {}

  // 4 correspondences are needed to determine a solution.
  double SampleSize() const { return 4; }

  // Estimates candidate absolute poses from correspondences.
  bool EstimateModel(
      const std::vector<FeatureCorrespondence2D3D>& correspondences,
      std::vector<UncalibratedRadialAbsolutePose>* absolute_poses_and_calib) const {
    Eigen::Matrix<double, 3, 4> features_hom;
    Eigen::Matrix4d world_points_hom;
    features_hom.setOnes();
    world_points_hom.setOnes();
    for (int i = 0; i < 4; ++i)
    {
        features_hom.col(i).topRows<2>() = correspondences[i].feature;
        world_points_hom.col(i).topRows<3>() = correspondences[i].world_point;
    }
    std::vector<Eigen::Matrix3d> rotations;
    std::vector<Eigen::Vector3d> translations;
    std::vector<double> focal_lengths;
    std::vector<double> radial_distortions;
    // estimate pose, focal length and radial distortion
    FourPointsPoseFocalLengthRadialDistortion(features_hom, world_points_hom,
                                              rotations,translations,
                                              radial_distortions, focal_lengths);
    for (int s = 0; s < rotations.size(); ++s)
    {
        UncalibratedRadialAbsolutePose solution;
        solution.focal_length = focal_lengths[s];
        solution.radial_distortion = radial_distortions[s];
        solution.translation = translations[s];
        solution.rotation = rotations[s];
        absolute_poses_and_calib->emplace_back(solution);
    }

    return rotations.size() > 0;
  }

  // The error for a correspondences given an absolute pose. This is the squared
  // reprojection error.
  double Error(const FeatureCorrespondence2D3D& correspondence,
               const UncalibratedRadialAbsolutePose& absolute_pose) const {
    // The reprojected point is computed as R * X + t where R is the camera
    // rotation, t is the translation, and X is the 3D point.
    Eigen::Vector3d point_in_cam =
            absolute_pose.rotation * correspondence.world_point +
            absolute_pose.translation;
    if (point_in_cam[2] < 0.0)
        return std::numeric_limits<double>::max();
    Eigen::Vector2d point_in_cam_persp_div(point_in_cam[0] / point_in_cam[2],
                                           point_in_cam[1] / point_in_cam[2]);

    point_in_cam_persp_div *= absolute_pose.focal_length;
    Eigen::Vector2d distorted_point;
    // see also division_undistortion_camera_model.h
    const double r_u_sq = point_in_cam_persp_div[0] * point_in_cam_persp_div[0] +
                          point_in_cam_persp_div[1] * point_in_cam_persp_div[1];

    const double denom = 2.0 * absolute_pose.radial_distortion * r_u_sq;
    const double inner_sqrt = 1.0 - 4.0 * absolute_pose.radial_distortion * r_u_sq;

    // If the denominator is nearly zero then we can evaluate the distorted
    // coordinates as k or r_u^2 goes to zero. Both evaluate to the identity.
    if (std::abs(denom) < std::numeric_limits<double>::epsilon() || inner_sqrt < 0.0) {
      distorted_point[0] = point_in_cam_persp_div[0];
      distorted_point[1] = point_in_cam_persp_div[1];
    } else {
      const double scale = (1.0 - std::sqrt(inner_sqrt)) / denom;
      distorted_point[0] = point_in_cam_persp_div[0] * scale;
      distorted_point[1] = point_in_cam_persp_div[1] * scale;
    }

    return (distorted_point - correspondence.feature).squaredNorm();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(UncalibratedRadialAbsolutePoseEstimator);
};

}  // namespace

bool EstimateUncalibratedRadialAbsolutePose(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence2D3D>& normalized_correspondences,
    UncalibratedRadialAbsolutePose* absolute_pose,
    RansacSummary* ransac_summary) {
  UncalibratedRadialAbsolutePoseEstimator absolute_pose_estimator;
  std::unique_ptr<SampleConsensusEstimator<UncalibratedRadialAbsolutePoseEstimator> >
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                absolute_pose_estimator);
  // Estimate the absolute pose.
  const bool success = ransac->Estimate(normalized_correspondences,
                                        absolute_pose,
                                        ransac_summary);
  return success;
}

}  // namespace theia
