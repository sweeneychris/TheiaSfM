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

#include "theia/sfm/estimators/estimate_essential_matrix.h"

#include <Eigen/Core>
#include <vector>

#include "theia/alignment/alignment.h"
#include "theia/solvers/estimator.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/five_point_relative_pose.h"
#include "theia/sfm/pose/util.h"
#include "theia/util/util.h"

namespace theia {
namespace {

// An estimator for computing the essential matrix from 5 feature
// correspondences. The feature correspondences should be normalized
// by the focal length with the principal point at (0, 0).
class EssentialMatrixEstimator
    : public Estimator<FeatureCorrespondence, Eigen::Matrix3d> {
 public:
  EssentialMatrixEstimator() {}

  // 5 correspondences are needed to determine an essential matrix.
  double SampleSize() const { return 5; }

  // Estimates candidate essential matrices from correspondences.
  bool EstimateModel(const std::vector<FeatureCorrespondence>& correspondences,
                     std::vector<Eigen::Matrix3d>* essential_matrices) const {
    std::vector<Eigen::Vector2d> image1_points, image2_points;
    image1_points.reserve(correspondences.size());
    image2_points.reserve(correspondences.size());
    for (int i = 0; i < correspondences.size(); i++) {
      image1_points.emplace_back(correspondences[i].feature1);
      image2_points.emplace_back(correspondences[i].feature2);
    }

    return FivePointRelativePose(image1_points,
                                 image2_points,
                                 essential_matrices);
  }

  // The error for a correspondences given a model. This is the squared sampson
  // error.
  double Error(const FeatureCorrespondence& correspondence,
               const Eigen::Matrix3d& essential_matrix) const {
    return SquaredSampsonDistance(essential_matrix,
                                  correspondence.feature1,
                                  correspondence.feature2);
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(EssentialMatrixEstimator);
};

}  // namespace

bool EstimateEssentialMatrix(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence>& normalized_correspondences,
    Eigen::Matrix3d* essential_matrix,
    RansacSummary* ransac_summary) {
  EssentialMatrixEstimator essential_matrix_estimator;
  std::unique_ptr<SampleConsensusEstimator<EssentialMatrixEstimator> >
      ransac = CreateAndInitializeRansacVariant(ransac_type,
                                                ransac_params,
                                                essential_matrix_estimator);

  // Estimate essential matrix.
  return ransac->Estimate(normalized_correspondences,
                          essential_matrix,
                          ransac_summary);
}

}  // namespace theia
