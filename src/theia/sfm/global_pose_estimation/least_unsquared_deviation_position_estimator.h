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

#ifndef THEIA_SFM_GLOBAL_POSE_ESTIMATION_LEAST_UNSQUARED_DEVIATION_POSITION_ESTIMATOR_H_
#define THEIA_SFM_GLOBAL_POSE_ESTIMATION_LEAST_UNSQUARED_DEVIATION_POSITION_ESTIMATOR_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unordered_map>

#include "theia/util/hash.h"
#include "theia/util/util.h"
#include "theia/sfm/global_pose_estimation/position_estimator.h"
#include "theia/sfm/types.h"

namespace theia {

// Estimates the camera position of views given pairwise relative poses and the
// absolute orientations of cameras. Positions are estimated using a least
// unsquared deviations solver -- essentially an L1 solver that is wrapped in an
// Iteratively Reweighted Least Squares (IRLS) formulation. This method was
// proposed in "Robust Camera Location Estimation by Convex Programming" by
// Ozyesil and Singer (CVPR 2015). Please cite this paper when using this
// method.
class LeastUnsquaredDeviationPositionEstimator : public PositionEstimator {
 public:
  struct Options {
    // Options for ADMM QP solver.
    int max_num_iterations = 400;

    // Maximum number of reweighted iterations.
    int max_num_reweighted_iterations = 10;

    // A measurement for convergence criterion.
    double convergence_criterion = 1e-4;
  };

  LeastUnsquaredDeviationPositionEstimator(
      const LeastUnsquaredDeviationPositionEstimator::Options& options);

  // Returns true if the optimization was a success, false if there was a
  // failure.
  bool EstimatePositions(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientation,
      std::unordered_map<ViewId, Eigen::Vector3d>* positions);

 private:
  void InitializeIndexMapping(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientations);

  // Creates camera to camera constraints from relative translations.
  void SetupConstraintMatrix(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientations);

  // Computes the weight of the error terms for the IRLS system.
  void UpdateConstraintWeights();

  const LeastUnsquaredDeviationPositionEstimator::Options options_;

  std::unordered_map<ViewIdPair, int> view_id_pair_to_index_;
  std::unordered_map<ViewId, int> view_id_to_index_;
  static const int kConstantViewIndex = -3;

  Eigen::ArrayXd weights_;
  Eigen::SparseMatrix<double> constraint_matrix_;
  Eigen::VectorXd solution_;

  friend class EstimatePositionsLeastUnsquaredDeviationTest;

  DISALLOW_COPY_AND_ASSIGN(LeastUnsquaredDeviationPositionEstimator);
};

}  // namespace theia

#endif  // THEIA_SFM_GLOBAL_POSE_ESTIMATION_LEAST_UNSQUARED_DEVIATION_POSITION_ESTIMATOR_H_
