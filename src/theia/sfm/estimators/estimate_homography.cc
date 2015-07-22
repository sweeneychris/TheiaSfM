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

#include "theia/sfm/estimators/estimate_homography.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <memory>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/pose/four_point_homography.h"
#include "theia/sfm/pose/util.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

// An estimator for computing a homography from 4 feature correspondences. The
// feature correspondences should be normalized by the focal length with the
// principal point at (0, 0).
class HomographyEstimator
    : public Estimator<FeatureCorrespondence, Eigen::Matrix3d> {
 public:
  HomographyEstimator() {}

  // 4 correspondences are needed to determine a homography.
  double SampleSize() const { return 4; }

  // Estimates candidate relative poses from correspondences.
  bool EstimateModel(const std::vector<FeatureCorrespondence>& correspondences,
                     std::vector<Eigen::Matrix3d>* homography) const {
    std::vector<Eigen::Vector2d> image1_points(4), image2_points(4);
    for (int i = 0; i < 4; i++) {
      image1_points[i] = correspondences[i].feature1;
      image2_points[i] = correspondences[i].feature2;
    }

    Eigen::Matrix3d homography_matrix;
    if (!FourPointHomography(image1_points,
                             image2_points,
                             &homography_matrix)) {
      return false;
    }

    homography->emplace_back(homography_matrix);
    return true;
  }

  // The error for a correspondences given a model. This is the asymmetric
  // distance that measures reprojection error in one image.
  double Error(const FeatureCorrespondence& correspondence,
               const Eigen::Matrix3d& homography) const {
    const Eigen::Vector3d reprojected_point =
        homography * correspondence.feature1.homogeneous();
    return (correspondence.feature2 - reprojected_point.hnormalized())
        .squaredNorm();
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(HomographyEstimator);
};

}  // namespace

bool EstimateHomography(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence>& correspondences,
    Eigen::Matrix3d* homography,
    RansacSummary* ransac_summary) {
  HomographyEstimator homography_estimator;
  std::unique_ptr<SampleConsensusEstimator<HomographyEstimator> > ransac =
      CreateAndInitializeRansacVariant(ransac_type,
                                       ransac_params,
                                       homography_estimator);
  // Estimate the homography.
  return ransac->Estimate(correspondences, homography, ransac_summary);
}

}  // namespace theia
