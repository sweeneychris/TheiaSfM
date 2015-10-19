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

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Core>
#include <unordered_map>

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
  CHECK_GT(options_.num_threads, 0);
}

bool LeastUnsquaredDeviationPositionEstimator::EstimatePositions(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  CHECK_NOTNULL(positions);
  view_pairs_ = &view_pairs;

  // Initialize scale and weight terms.
  scales_.reserve(view_pairs_->size());
  weights_.reserve(view_pairs_->size());
  for (const auto& view_pair : *view_pairs_) {
    scales_[view_pair.first] = 10.0;
    weights_[view_pair.first] = 1.0;
  }

  if (options_.initialize_random_positions) {
    // Initialize positions to be random.
    InitializeRandomPositions(orientations, positions);
  }

  ceres::Solver::Summary summary;
  for (int i = 0; i < options_.max_num_reweighted_iterations; i++) {
    problem_.reset(new ceres::Problem());

    // Add the (weighted) constraints to the problem.
    AddCameraToCameraConstraints(orientations, positions);

    // Set one camera to be at the origin to remove the ambiguity of the origin.
    positions->begin()->second.setZero();
    problem_->SetParameterBlockConstant(positions->begin()->second.data());

    // Set the solver options.
    solver_options_.num_threads = options_.num_threads;
    solver_options_.num_linear_solver_threads = options_.num_threads;
    solver_options_.max_num_iterations = options_.max_num_iterations;
    solver_options_.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;

    ceres::Solve(solver_options_, problem_.get(), &summary);
    LOG(INFO) << "Iteration " << i + 1
              << " of the Iterative Reweighted Least Squares\n"
              << summary.FullReport();

    const double convergence_criterion =
        (summary.initial_cost - summary.final_cost) / summary.initial_cost;
    if (convergence_criterion < options_.convergence_criterion) {
      break;
    }

    ComputeWeights(orientations, *positions);
  }

  return summary.IsSolutionUsable();
}

void LeastUnsquaredDeviationPositionEstimator::InitializeRandomPositions(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  std::unordered_set<ViewId> constrained_positions;
  constrained_positions.reserve(orientations.size());
  for (const auto& view_pair : *view_pairs_) {
    constrained_positions.insert(view_pair.first.first);
    constrained_positions.insert(view_pair.first.second);
  }

  positions->reserve(orientations.size());
  for (const auto& orientation : orientations) {
    if (ContainsKey(constrained_positions, orientation.first)) {
      (*positions)[orientation.first] = 100.0 * Vector3d::Random();
    }
  }
}

void LeastUnsquaredDeviationPositionEstimator::AddCameraToCameraConstraints(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  // Add the camera to camera constraints.
  for (const auto& view_pair : *view_pairs_) {
    const ViewId view_id1 = view_pair.first.first;
    const ViewId view_id2 = view_pair.first.second;
    Vector3d* position1 = FindOrNull(*positions, view_id1);
    Vector3d* position2 = FindOrNull(*positions, view_id2);

    // Do not add this view pair if one or both of the positions do not exist.
    if (position1 == nullptr || position2 == nullptr) {
      continue;
    }

    // Rotate the relative translation so that it is aligned to the global
    // orientation frame.
    const Vector3d translation_direction = GetRotatedTranslation(
        FindOrDie(orientations, view_id1), view_pair.second.position_2);

    ceres::CostFunction* cost_function =
        PairwiseTranslationAndScaleError::Create(translation_direction);

    double* scale = FindOrNull(scales_, view_pair.first);
    const double weight = FindOrDieNoPrint(weights_, view_pair.first);
    problem_->AddResidualBlock(
        cost_function,
        new ceres::ScaledLoss(NULL, weight, ceres::TAKE_OWNERSHIP),
        position1->data(),
        position2->data(),
        scale);
    problem_->SetParameterLowerBound(scale, 0, 1.0);
  }

  VLOG(2) << problem_->NumResidualBlocks() << " camera to camera constraints "
                                              "were added to the position "
                                              "estimation problem.";
}

void LeastUnsquaredDeviationPositionEstimator::ComputeWeights(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    const std::unordered_map<ViewId, Vector3d>& positions) {
  const double gamma = 1e-8;
  // Add the camera to camera constraints.
  for (const auto& view_pair : *view_pairs_) {
    const ViewId view_id1 = view_pair.first.first;
    const ViewId view_id2 = view_pair.first.second;
    const Vector3d* position1 = FindOrNull(positions, view_id1);
    const Vector3d* position2 = FindOrNull(positions, view_id2);

    // Do not add this view pair if one or both of the positions do not exist.
    if (position1 == nullptr || position2 == nullptr) {
      continue;
    }

    // Rotate the relative translation so that it is aligned to the global
    // orientation frame.
    const Vector3d translation_direction = GetRotatedTranslation(
        FindOrDie(orientations, view_id1), view_pair.second.position_2);
    const double scale = FindOrDieNoPrint(scales_, view_pair.first);
    const double error =
        (*position2 - *position1 - scale * translation_direction).squaredNorm();
    weights_[view_pair.first] = std::sqrt(1.0 / (error + gamma));
  }
}

}  // namespace theia
