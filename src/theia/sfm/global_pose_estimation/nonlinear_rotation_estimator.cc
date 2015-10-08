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

#include "theia/sfm/global_pose_estimation/nonlinear_rotation_estimator.h"

#include <ceres/ceres.h>
#include <Eigen/Core>
#include <memory>
#include <unordered_map>

#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/sfm/global_pose_estimation/pairwise_rotation_error.h"
#include "theia/sfm/types.h"

namespace theia {

bool NonlinearRotationEstimator::EstimateRotations(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations) {
  CHECK_NOTNULL(global_orientations);
  if (global_orientations->size() == 0) {
    LOG(INFO) << "Skipping nonlinear rotation optimization because no "
                 "initialization was provivded.";
    return false;
  }
  if (view_pairs.size() == 0) {
    LOG(INFO) << "Skipping nonlinear rotation optimization because no "
                 "relative rotation constraints were provivded.";
    return false;
  }

  // Set up the problem and loss function.
  std::unique_ptr<ceres::Problem> problem(new ceres::Problem());
  ceres::LossFunction* loss_function =
      new ceres::SoftLOneLoss(robust_loss_width_);

  for (const auto& view_pair : view_pairs) {
    const ViewIdPair& view_id_pair = view_pair.first;
    Eigen::Vector3d* rotation1 =
        FindOrNull(*global_orientations, view_id_pair.first);
    Eigen::Vector3d* rotation2 =
        FindOrNull(*global_orientations, view_id_pair.second);

    // Do not add the relative rotation constaint if it requires an orientation
    // that we do not have an initialization for.
    if (rotation1 == nullptr || rotation2 == nullptr) {
      continue;
    }

    ceres::CostFunction* cost_function =
        PairwiseRotationError::Create(view_pair.second.rotation_2, 1.0);
    problem->AddResidualBlock(cost_function,
                              loss_function,
                              rotation1->data(),
                              rotation2->data());
  }

  // The problem should be relatively sparse so sparse cholesky is a good
  // choice.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.max_num_iterations = 200;

  ceres::Solver::Summary summary;
  ceres::Solve(options, problem.get(), &summary);
  VLOG(1) << summary.FullReport();
  return true;
}

}  // namespace theia
