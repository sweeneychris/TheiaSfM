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

#include "theia/sfm/global_pose_estimation/least_unsquared_deviation_position_estimator.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <limits>
#include <unordered_map>
#include <vector>

#include "theia/math/qp_solver.h"
#include "theia/sfm/global_pose_estimation/pairwise_translation_and_scale_error.h"
#include "theia/sfm/types.h"
#include "theia/util/map_util.h"
#include "theia/util/util.h"

namespace theia {
namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

Vector3d GetRotatedTranslation(const Vector3d& rotation_angle_axis,
                               const Vector3d& translation) {
  Matrix3d rotation;
  ceres::AngleAxisToRotationMatrix(
      rotation_angle_axis.data(),
      ceres::ColumnMajorAdapter3x3(rotation.data()));
  return rotation.transpose() * translation;
}

}  // namespace

LeastUnsquaredDeviationPositionEstimator::
    LeastUnsquaredDeviationPositionEstimator(
        const LeastUnsquaredDeviationPositionEstimator::Options& options)
    : options_(options) {
  CHECK_GT(options_.max_num_iterations, 0);
  CHECK_GT(options_.max_num_reweighted_iterations, 0);
}

bool LeastUnsquaredDeviationPositionEstimator::EstimatePositions(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  CHECK_NOTNULL(positions)->clear();

  InitializeIndexMapping(view_pairs, orientations);

  // Set up the linear system.
  SetupConstraintMatrix(view_pairs, orientations);
  solution_.setZero(constraint_matrix_.cols());
  weights_.setConstant(constraint_matrix_.rows(), 1.0);

  // Set the lower bounds for the QP. The positions should be unbounded while
  // the scales havea lower bound of 1.0.
  Eigen::VectorXd lower_bound(constraint_matrix_.cols());
  lower_bound.head(3 * orientations.size())
      .setConstant(-std::numeric_limits<double>::infinity());
  lower_bound.tail(view_pairs.size()).setConstant(1.0);

  // Solve for the camera positions using an IRLS scheme. The values p and r are
  // constant at zero.
  Eigen::SparseMatrix<double> P(constraint_matrix_.rows(),
                                constraint_matrix_.rows());
  Eigen::VectorXd q(constraint_matrix_.cols());
  q.setZero();
  const double r = 0;
  QPSolver::Options qp_solver_options;
  qp_solver_options.max_num_iterations = 10;
  for (int i = 0; i < options_.max_num_reweighted_iterations; i++) {
    if (i > 0) {
      UpdateConstraintWeights();
    }

    // Compute P = A^t * W * A, the quadratic matrix term for our QP.
    P = constraint_matrix_.transpose() * weights_.matrix().asDiagonal() *
        constraint_matrix_;

    // Solve the quadratic program. Increase the number of possible iterations
    // each time.
    qp_solver_options.max_num_iterations = std::min(
        options_.max_num_iterations, 2 * qp_solver_options.max_num_iterations);
    QPSolver qp_solver(qp_solver_options, P, q, r);
    const Eigen::VectorXd prev_solution = solution_;
    qp_solver.SetLowerBound(lower_bound);
    if (!qp_solver.Solve(&solution_)) {
      LOG(WARNING) << "Could not solve the Quadratic Program for the least "
                      "unsquared deviations position solver.";
      return false;
    }
  }

  // Set the estimated positions.
  for (const auto& view_id_index : view_id_to_index_) {
    const int index = view_id_index.second;
    const ViewId view_id = view_id_index.first;
    if (index == kConstantViewIndex) {
      (*positions)[view_id] = Eigen::Vector3d::Zero();
    } else {
      (*positions)[view_id] = solution_.segment<3>(index);
    }
  }

  return true;
}

void LeastUnsquaredDeviationPositionEstimator::InitializeIndexMapping(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations) {
  std::unordered_set<ViewId> views;
  for (const auto& view_pair : view_pairs) {
    if (ContainsKey(orientations, view_pair.first.first)) {
      views.insert(view_pair.first.first);
    }
    if (ContainsKey(orientations, view_pair.first.second)) {
      views.insert(view_pair.first.second);
    }
  }

  // Create a mapping from the view id to the index of the linear system.
  int index = kConstantViewIndex;
  view_id_to_index_.reserve(orientations.size());
  for (const ViewId view_id : views) {
    view_id_to_index_[view_id] = index;
    index += 3;
  }

  // Create a mapping from the view id pair to the index of the linear system.
  view_id_pair_to_index_.reserve(view_pairs.size());
  for (const auto& view_pair : view_pairs) {
    view_id_pair_to_index_[view_pair.first] = index;
    ++index;
  }
}

void LeastUnsquaredDeviationPositionEstimator::SetupConstraintMatrix(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations) {
  constraint_matrix_.resize(3 * view_pairs.size(),
                            3 * orientations.size() + view_pairs.size());

  // Add the camera to camera constraints.
  std::vector<Eigen::Triplet<double> > triplet_list;
  triplet_list.reserve(9 * view_pairs.size());
  int row = 0;
  for (const auto& view_pair : view_pairs) {
    const ViewIdPair view_id_pair = view_pair.first;

    const int view1_index = FindOrDie(view_id_to_index_, view_id_pair.first);
    const int view2_index = FindOrDie(view_id_to_index_, view_id_pair.second);
    const int scale_index =
        FindOrDieNoPrint(view_id_pair_to_index_, view_id_pair);

    // Rotate the relative translation so that it is aligned to the global
    // orientation frame.
    const Vector3d translation_direction =
        GetRotatedTranslation(FindOrDie(orientations, view_id_pair.first),
                              view_pair.second.position_2);

    // Add the constraint for view 1 in the minimization:
    //   position2 - position1 - scale_1_2 * translation_direction.
    if (view1_index != kConstantViewIndex) {
      triplet_list.emplace_back(row + 0, view1_index + 0, -1.0);
      triplet_list.emplace_back(row + 1, view1_index + 1, -1.0);
      triplet_list.emplace_back(row + 2, view1_index + 2, -1.0);
    }

    // Add the constraint for view 2 in the minimization:
    //   position2 - position1 - scale_1_2 * translation_direction.
    if (view2_index != kConstantViewIndex) {
      triplet_list.emplace_back(row + 0, view2_index + 0, 1.0);
      triplet_list.emplace_back(row + 1, view2_index + 1, 1.0);
      triplet_list.emplace_back(row + 2, view2_index + 2, 1.0);
    }

    // Add the constraint for scale in the minimization:
    //   position2 - position1 - scale_1_2 * translation_direction.
    triplet_list.emplace_back(row + 0, scale_index, -translation_direction[0]);
    triplet_list.emplace_back(row + 1, scale_index, -translation_direction[1]);
    triplet_list.emplace_back(row + 2, scale_index, -translation_direction[2]);

    row += 3;
  }

  constraint_matrix_.setFromTriplets(triplet_list.begin(), triplet_list.end());

  VLOG(2) << view_pairs.size() << " camera to camera constraints were added "
                                    "to the position estimation problem.";
}

// Compute the error:
//     err_i_j = || c_j - c_i - scale_i_j * t_i_j ||^2
//
// For each pairwise constraint, set w_i_j = (err_i_j + delta)^(-1/2) to
// reweight the QP solver for robustness.
void LeastUnsquaredDeviationPositionEstimator::UpdateConstraintWeights() {
  static const double delta = 1e-12;
  // Compute the errors with a simple matrix multiplication.
  const Eigen::VectorXd errors = constraint_matrix_ * solution_;
  weights_ = (errors.array().square() + delta).sqrt().inverse();
}

}  // namespace theia
