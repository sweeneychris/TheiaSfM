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

#include "theia/sfm/bundle_adjustment/optimize_relative_position_with_known_rotation.h"

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Core>
#include <vector>

#include "theia/math/util.h"
#include "theia/sfm/bundle_adjustment/orthogonal_vector_error.h"
#include "theia/sfm/bundle_adjustment/unit_norm_three_vector_parameterization.h"
#include "theia/matching/feature_correspondence.h"

namespace theia {

// The epipolar constraint x1' * [t]_x * R * x0 = 0 may be rewritten as:
//   -(x1.cross(R * x0))' * t = 0, or alternatively v' * t = 0 where
// v = -x1.cross(R * x0). We can use this parameterization to optimize for the
// relative position based on feature correspondences.
bool OptimizeRelativePositionWithKnownRotation(
    const std::vector<FeatureCorrespondence>& correspondences,
    const Eigen::Vector3d& relative_rotation,
    Eigen::Vector3d* relative_position) {
  CHECK_NOTNULL(relative_position);

  // Set problem options.
  ceres::Problem::Options problem_options;
  // Only use one loss function for the entire problem.
  std::unique_ptr<ceres::LossFunction> loss_function(new ceres::HuberLoss(0.1));
  problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  ceres::Problem problem(problem_options);

  ceres::Solver::Options solver_options;
  solver_options.linear_solver_type = ceres::DENSE_QR;
  solver_options.logging_type = ceres::SILENT;

  // Add position as a parameter with unit norm.
  // Add the position as a parameter block, ensuring that the norm is 1.
  const int kPositionSize = 3;
  ceres::LocalParameterization* position_parameterization =
      new ceres::AutoDiffLocalParameterization<
          UnitNormThreeVectorParameterization, 3, 3>;
  problem.AddParameterBlock(relative_position->data(),
                            kPositionSize,
                            position_parameterization);

  // Add correspondences to the optimization problem.
  Eigen::Matrix3d rotation_matrix;
  ceres::AngleAxisToRotationMatrix(
      relative_rotation.data(),
      ceres::ColumnMajorAdapter3x3(rotation_matrix.data()));

  for (const FeatureCorrespondence& match : correspondences) {
    const Eigen::Vector3d feature1 = match.feature1.homogeneous();
    const Eigen::Vector3d feature2 = match.feature2.homogeneous();
    const Eigen::Vector3d rotated_feature1 = rotation_matrix * feature1;

    const Eigen::Vector3d orthogonal_vector =
        feature2.cross(rotated_feature1).transpose() * rotation_matrix;
    problem.AddResidualBlock(
        OrthogonalVectorError::Create(orthogonal_vector),
        loss_function.get(),
        relative_position->data());
  }

  // Solve.
  ceres::Solver::Summary summary;
  ceres::Solve(solver_options, &problem, &summary);

  VLOG(2) << summary.FullReport();

  return summary.termination_type != ceres::FAILURE;
}

}  // namespace theia
