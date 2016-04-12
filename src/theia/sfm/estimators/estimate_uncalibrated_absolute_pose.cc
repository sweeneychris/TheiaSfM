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

#include "theia/sfm/estimators/estimate_uncalibrated_absolute_pose.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <vector>

#include "theia/sfm/camera/projection_matrix_utils.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/estimators/feature_correspondence_2d_3d.h"
#include "theia/sfm/pose/four_point_focal_length.h"
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
// the principal point is at (0, 0).
class UncalibratedAbsolutePoseEstimator
    : public Estimator<FeatureCorrespondence2D3D, Matrix3x4d> {
 public:
  UncalibratedAbsolutePoseEstimator() {}

  // 3 correspondences are needed to determine the absolute pose.
  double SampleSize() const { return 4; }

  // Estimates candidate absolute poses from correspondences.
  bool EstimateModel(
      const std::vector<FeatureCorrespondence2D3D>& correspondences,
      std::vector<Matrix3x4d>* absolute_poses) const {
    const std::vector<Eigen::Vector2d> features = {
        correspondences[0].feature, correspondences[1].feature,
        correspondences[2].feature, correspondences[3].feature};
    const std::vector<Eigen::Vector3d> world_points = {
        correspondences[0].world_point, correspondences[1].world_point,
        correspondences[2].world_point, correspondences[3].world_point};

    const int num_solutions =
        FourPointPoseAndFocalLength(features, world_points, absolute_poses);
    return num_solutions > 0;
  }

  // The error for a correspondences given an absolute pose. This is the squared
  // reprojection error.
  double Error(const FeatureCorrespondence2D3D& correspondence,
               const Matrix3x4d& absolute_pose) const {
    // The reprojected point is computed as R * (X - c) where R is the camera
    // rotation, c is the position, and X is the 3D point.
    const Eigen::Vector2d reprojected_feature =
        (absolute_pose * correspondence.world_point.homogeneous())
            .eval()
            .hnormalized();
    return (reprojected_feature - correspondence.feature).squaredNorm();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(UncalibratedAbsolutePoseEstimator);
};

}  // namespace

bool EstimateUncalibratedAbsolutePose(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence2D3D>& normalized_correspondences,
    UncalibratedAbsolutePose* absolute_pose,
    RansacSummary* ransac_summary) {
  UncalibratedAbsolutePoseEstimator absolute_pose_estimator;
  std::unique_ptr<SampleConsensusEstimator<UncalibratedAbsolutePoseEstimator> >
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                absolute_pose_estimator);
  // Estimate the absolute pose.
  Matrix3x4d projection_matrix;
  const bool success = ransac->Estimate(normalized_correspondences,
                                        &projection_matrix,
                                        ransac_summary);

  // Recover the focal length and pose.
  Eigen::Matrix3d calibration_matrix;
  Eigen::Matrix3d rotation_matrix;
  DecomposeProjectionMatrix(projection_matrix,
                            &calibration_matrix,
                            &absolute_pose->rotation,
                            &absolute_pose->position);

  absolute_pose->focal_length =
      calibration_matrix(0, 0) / calibration_matrix(2, 2);

  return success;
}

}  // namespace theia
