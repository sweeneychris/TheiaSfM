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

#ifndef THEIA_SFM_GLOBAL_POSE_ESTIMATION_ROBUST_ROTATION_ESTIMATOR_H_
#define THEIA_SFM_GLOBAL_POSE_ESTIMATION_ROBUST_ROTATION_ESTIMATOR_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unordered_map>

#include "theia/sfm/global_pose_estimation/rotation_estimator.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/util/hash.h"

namespace theia {

// Computes the global rotations given relative rotations and an initial guess
// for the global orientations. The robust algorithm of "Efficient and Large
// Scale Rotation Averaging" by Chatterjee and Govindu (ICCV 2013) is used to
// obtain accurate solutions that are robust to outliers.
//
// The general strategy of this algorithm is to minimize the relative rotation
// error (using the difference between relative rotations and the corresponding
// global rotations) with L1 minimization first, then a reweighted least
// squares. The L1 minimization is relatively slow, but provides excellent
// robustness to outliers. Then the L2 minimization (which is much faster) can
// refine the solution to be very accurate.
class RobustRotationEstimator : public RotationEstimator {
 public:
  struct Options {
    // Maximum number of times to run L1 minimization. L1 is very slow (compared
    // to L2), but is very robust to outliers. Typically only a few iterations
    // are needed in order for the solution to reside within the cone of
    // convergence for L2 solving.
    int max_num_l1_iterations = 5;

    // The number of iterative reweighted least squares iterations to perform.
    int max_num_irls_iterations = 100;
  };

  explicit RobustRotationEstimator(const Options& options)
      : options_(options) {}

  // Estimates the global orientations of all views based on an initial
  // guess. Returns true on successful estimation and false otherwise.
  bool EstimateRotations(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations);

 protected:
  // Sets up the sparse linear system such that dR_ij = dR_j - dR_i. This is the
  // first-order approximation of the angle-axis rotations. This should only be
  // called once.
  void SetupLinearSystem();

  // Computes the relative rotation error based on the current global
  // orientation estimates.
  void ComputeRotationError();

  // Performs the L1 robust loss minimization.
  bool SolveL1Regression();

  // Performs the iteratively reweighted least squares.
  bool SolveIRLS();

  // Updates the global rotations based on the current rotation change.
  void UpdateGlobalRotations();

  // We keep one of the rotations as constant to remove the ambiguity of the
  // linear system.
  static const int kConstantRotationIndex = -1;

  // Options for the solver.
  const Options options_;

  // The pairwise relative rotations used to compute the global rotations.
  const std::unordered_map<ViewIdPair, TwoViewInfo>* view_pairs_;

  // The global orientation estimates for each camera.
  std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations_;

  // The sparse matrix used to maintain the linear system. This is matrix A in
  // Ax = b.
  Eigen::SparseMatrix<double> sparse_matrix_;

  // Map of ViewIds to the corresponding positions of the view's orientation in
  // the linear system.
  std::unordered_map<ViewId, int> view_id_to_index_;

  // x in the linear system Ax = b.
  Eigen::VectorXd rotation_change_;

  // b in the linear system Ax = b.
  Eigen::VectorXd relative_rotation_error_;
};

}  // namespace theia

#endif  // THEIA_SFM_GLOBAL_POSE_ESTIMATION_ROBUST_ROTATION_ESTIMATOR_H_
