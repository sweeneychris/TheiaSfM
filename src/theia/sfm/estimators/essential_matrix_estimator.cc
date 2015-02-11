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

#include "theia/sfm/estimators/essential_matrix_estimator.h"

#include <Eigen/Core>
#include <vector>

#include "theia/solvers/estimator.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/five_point_relative_pose.h"
#include "theia/sfm/pose/util.h"

namespace theia {

bool EssentialMatrixEstimator::EstimateModel(
    const std::vector<FeatureCorrespondence>& correspondences,
    std::vector<Eigen::Matrix3d>* essential_matrices) const {
  Eigen::Vector2d image1_points[5], image2_points[5];
  for (int i = 0; i < 5; i++) {
    image1_points[i] = correspondences[i].feature1;
    image2_points[i] = correspondences[i].feature2;
  }

  return FivePointRelativePose(image1_points,
                               image2_points,
                               essential_matrices);
}

double EssentialMatrixEstimator::Error(
    const FeatureCorrespondence& correspondence,
    const Eigen::Matrix3d& essential_matrix) const {
  return SquaredSampsonDistance(essential_matrix,
                                correspondence.feature1,
                                correspondence.feature2);
}

}  // namespace theia
