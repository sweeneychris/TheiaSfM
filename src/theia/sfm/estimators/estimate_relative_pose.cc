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

#include "theia/sfm/estimators/estimate_relative_pose.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <memory>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/pose/essential_matrix_utils.h"
#include "theia/sfm/pose/five_point_relative_pose.h"
#include "theia/sfm/pose/util.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

// An estimator for computing the relative pose from 5 feature
// correspondences. The feature correspondences should be normalized
// by the focal length with the principal point at (0, 0).
class RelativePoseEstimator
    : public Estimator<FeatureCorrespondence, RelativePose> {
 public:
  RelativePoseEstimator() {}

  // 5 correspondences are needed to determine an essential matrix and thus a
  // relative pose..
  double SampleSize() const { return 5; }

  // Estimates candidate relative poses from correspondences.
  bool EstimateModel(const std::vector<FeatureCorrespondence>& correspondences,
                     std::vector<RelativePose>* relative_poses) const {
    std::vector<Eigen::Vector2d> image1_points, image2_points;
    image1_points.reserve(correspondences.size());
    image2_points.reserve(correspondences.size());
    for (int i = 0; i < correspondences.size(); i++) {
      image1_points.emplace_back(correspondences[i].feature1);
      image2_points.emplace_back(correspondences[i].feature2);
    }

    std::vector<Matrix3d> essential_matrices;
    if (!FivePointRelativePose(image1_points,
                               image2_points,
                               &essential_matrices)) {
      return false;
    }

    relative_poses->reserve(essential_matrices.size() * 4);
    for (const Eigen::Matrix3d& essential_matrix : essential_matrices) {
      RelativePose relative_pose;
      relative_pose.essential_matrix = essential_matrix;

      // The best relative pose decomposition should have at least 4
      // triangulated points in front of the camera. This is because one point
      // may be at infinity.
      const int num_points_in_front_of_cameras = GetBestPoseFromEssentialMatrix(
          essential_matrix,
          correspondences,
          &relative_pose.rotation,
          &relative_pose.position);
      if (num_points_in_front_of_cameras >= 4) {
        relative_poses->push_back(relative_pose);
      }
    }
    return relative_poses->size() > 0;
  }

  // The error for a correspondences given a model. This is the squared sampson
  // error.
  double Error(const FeatureCorrespondence& correspondence,
               const RelativePose& relative_pose) const {
    if (IsTriangulatedPointInFrontOfCameras(correspondence,
                                            relative_pose.rotation,
                                            relative_pose.position)) {
      return SquaredSampsonDistance(relative_pose.essential_matrix,
                                    correspondence.feature1,
                                    correspondence.feature2);
    }
    return std::numeric_limits<double>::max();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(RelativePoseEstimator);
};

}  // namespace

bool EstimateRelativePose(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence>& normalized_correspondences,
    RelativePose* relative_pose,
    RansacSummary* ransac_summary) {
  RelativePoseEstimator relative_pose_estimator;
  std::unique_ptr<SampleConsensusEstimator<RelativePoseEstimator> > ransac =
      CreateAndInitializeRansacVariant(ransac_type,
                                       ransac_params,
                                       relative_pose_estimator);
  // Estimate the relative pose.
  return ransac->Estimate(normalized_correspondences,
                          relative_pose,
                          ransac_summary);
}

}  // namespace theia
