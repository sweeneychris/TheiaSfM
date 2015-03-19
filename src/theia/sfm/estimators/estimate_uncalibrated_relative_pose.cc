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

#include "theia/sfm/estimators/estimate_uncalibrated_relative_pose.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <memory>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/pose/eight_point_fundamental_matrix.h"
#include "theia/sfm/pose/essential_matrix_utils.h"
#include "theia/sfm/pose/fundamental_matrix_util.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace {

// An estimator for computing the relative pose from 8 feature correspondences
// (via decomposition of the fundamental matrix).
//
// NOTE: Feature correspondences must be in pixel coordinates with the principal
// point removed i.e. principal point at (0, 0). This also assumes negligible
// skew (which is reasonable for most cameras).
class UncalibratedRelativePoseEstimator
    : public Estimator<FeatureCorrespondence, UncalibratedRelativePose> {
 public:
  UncalibratedRelativePoseEstimator() {}

  // 8 correspondences are needed to determine a fundamental matrix and thus a
  // relative pose.
  double SampleSize() const { return 8; }

  // Estimates candidate relative poses from correspondences.
  bool EstimateModel(
      const std::vector<FeatureCorrespondence>& centered_correspondences,
      std::vector<UncalibratedRelativePose>* relative_poses) const {
    std::vector<Eigen::Vector2d> image1_points, image2_points;
    for (int i = 0; i < 8; i++) {
      image1_points.emplace_back(centered_correspondences[i].feature1);
      image2_points.emplace_back(centered_correspondences[i].feature2);
    }

    UncalibratedRelativePose relative_pose;
    if (!NormalizedEightPointFundamentalMatrix(
            image1_points, image2_points, &relative_pose.fundamental_matrix)) {
      return false;
    }

    // Only consider fundamental matrices that we can decompose focal lengths
    // from.
    if (!FocalLengthsFromFundamentalMatrix(
            relative_pose.fundamental_matrix.data(),
            &relative_pose.focal_length1,
            &relative_pose.focal_length2)) {
      return false;
    }

    // TODO(cmsweeney): Should we check if the focal lengths are reasonable?

    // Compose the essential matrix from the fundamental matrix and focal
    // lengths.
    const Matrix3d essential_matrix =
        Eigen::DiagonalMatrix<double, 3>(relative_pose.focal_length2,
                                         relative_pose.focal_length2,
                                         1.0) *
        relative_pose.fundamental_matrix *
        Eigen::DiagonalMatrix<double, 3>(relative_pose.focal_length1,
                                         relative_pose.focal_length1,
                                         1.0);

    // Normalize the centered_correspondences.
    std::vector<FeatureCorrespondence> normalized_correspondences(
        centered_correspondences.size());
    for (int i = 0; i < centered_correspondences.size(); i++) {
      normalized_correspondences[i].feature1 =
          centered_correspondences[i].feature1 / relative_pose.focal_length1;
      normalized_correspondences[i].feature2 =
          centered_correspondences[i].feature2 / relative_pose.focal_length2;
    }

    GetBestPoseFromEssentialMatrix(essential_matrix,
                                   normalized_correspondences,
                                   &relative_pose.rotation,
                                   &relative_pose.position);
    relative_poses->emplace_back(relative_pose);
    return true;
  }

  // The error for a correspondences given a model. This is the squared sampson
  // error.
  double Error(const FeatureCorrespondence& centered_correspondence,
               const UncalibratedRelativePose& relative_pose) const {
    return SquaredSampsonDistance(relative_pose.fundamental_matrix,
                                  centered_correspondence.feature1,
                                  centered_correspondence.feature2);
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(UncalibratedRelativePoseEstimator);
};

}  // namespace

bool EstimateUncalibratedRelativePose(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence>& centered_correspondences,
    UncalibratedRelativePose* relative_pose,
    RansacSummary* ransac_summary) {
  UncalibratedRelativePoseEstimator relative_pose_estimator;
  std::unique_ptr<SampleConsensusEstimator<UncalibratedRelativePoseEstimator> >
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                relative_pose_estimator);

  // Estimate essential matrix.
  return ransac->Estimate(centered_correspondences,
                          relative_pose,
                          ransac_summary);
}

}  // namespace theia
