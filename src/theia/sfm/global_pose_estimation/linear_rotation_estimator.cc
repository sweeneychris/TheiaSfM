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

#include "theia/sfm/global_pose_estimation/linear_rotation_estimator.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <unordered_map>

#include "theia/sfm/pose/util.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"

namespace theia {
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

typedef Eigen::Triplet<double> TripletEntry;

namespace {

void Fill3x3SparseMatrix(const Eigen::Matrix3d& mat,
                         const int row,
                         const int col,
                         std::vector<TripletEntry>* triplet_list) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      // Only add the value if it is not zero.
      if (mat(i, j) != 0.0) {
        triplet_list->emplace_back(row + i, col + j, mat(i, j));
      }
    }
  }
}

}  // namespace

// Set up linear system R_j - R_{i,j} * R_i = 0. This is a (3 * m) x (3 * n)
// system where m is the number of relative rotations and n is the number of
// views. Instead of the standard Ax = 0 formulation, we instead assume that one
// of the rotations is the identity and hold it constant. Using a bit of math,
// this ends up being a Cy = Z equation were C is all columns of A except the
// ones corresonding to the constant rotation. Z then, is the negative of the
// columns of A that were held constant. The constant rotation is indicated by
// an index of -1.
void LinearRotationEstimator::SetupSparseLinearSystem(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs) {
  static const int kRotationMatrixDimSize = 3;

  // Reserve size of matrices.
  lhs_.resize(kRotationMatrixDimSize * view_pairs.size(),
              kRotationMatrixDimSize * (view_id_map_.size() - 1));

  rhs_.resize(kRotationMatrixDimSize * view_pairs.size(),
              kRotationMatrixDimSize);
  rhs_.setZero();

  // Iterate through the relative rotation constraints and add them to the
  // sparse matrix.
  std::vector<TripletEntry> lhs_triplet_list, rhs_triplet_list;
  lhs_triplet_list.reserve(12 * view_pairs.size());
  int current_relative_rotation_index = 0;
  for (const auto& view_pair : view_pairs) {
    const int view1_index =
        FindOrDie(view_id_map_, view_pair.first.first);
    const int view2_index =
        FindOrDie(view_id_map_, view_pair.first.second);

    // Convert the relative rotation to a rotation matrix and store it.
    Eigen::Matrix3d relative_rotation;
    ceres::AngleAxisToRotationMatrix(
        view_pair.second.rotation_2.data(),
        ceres::ColumnMajorAdapter3x3(relative_rotation.data()));

    Eigen::Matrix3d identity_rotation = Eigen::Matrix3d::Identity();

    // Weight the terms by the inlier count, if desired.
    if (weight_terms_by_inliers_) {
      relative_rotation *= view_pair.second.num_verified_matches;
      identity_rotation *= view_pair.second.num_verified_matches;
    }

    // Set R_j term to identity.
    if (view2_index == -1) {
      rhs_.block<3, 3>(kRotationMatrixDimSize * current_relative_rotation_index,
                       0) = -identity_rotation;
    } else {
      Fill3x3SparseMatrix(
          identity_rotation,
          kRotationMatrixDimSize * current_relative_rotation_index,
          kRotationMatrixDimSize * view2_index,
          &lhs_triplet_list);
    }

    // Set R_i term to -R_{ij}.
    if (view1_index == -1) {
      rhs_.block<3, 3>(kRotationMatrixDimSize * current_relative_rotation_index,
                       0) = relative_rotation;
    } else {
      Fill3x3SparseMatrix(
          -relative_rotation,
          kRotationMatrixDimSize * current_relative_rotation_index,
          kRotationMatrixDimSize * view1_index,
          &lhs_triplet_list);
    }

    ++current_relative_rotation_index;
  }

  lhs_.setFromTriplets(lhs_triplet_list.begin(), lhs_triplet_list.end());
}

void LinearRotationEstimator::IndexInputViews(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs) {
  // Set one index to -1 because we keep it constant during the minimization.
  int current_index = -1;
  for (const auto& view_pair : view_pairs) {
    // Set the index of the view if it has not already been set.
    if (!ContainsKey(view_id_map_, view_pair.first.first)) {
      view_id_map_[view_pair.first.first] = current_index;
      ++current_index;
    }
    if (!ContainsKey(view_id_map_, view_pair.first.second)) {
      view_id_map_[view_pair.first.second] = current_index;
      ++current_index;
    }
  }

  return;
}

// Estimates the global orientations of all views based on an initial
// guess. Returns true on successful estimation and false otherwise.
bool LinearRotationEstimator::EstimateRotations(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations) {
  CHECK_GT(view_pairs.size(), 0);
  CHECK_NOTNULL(global_orientations);

  IndexInputViews(view_pairs);

  // Setup the sparse linear system.
  SetupSparseLinearSystem(view_pairs);

  // Setup sparse QR solver.
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
      sparse_qr_solver(lhs_);

  // Check that the input to SparseQR was valid and could be succesffully
  // factorized.
  if (sparse_qr_solver.info() != Eigen::Success) {
    VLOG(2) << "The sparse matrix of relative rotation constraints could not "
               "be factorized.";
    return false;
  }

  // Solve the linear least squares system.
  const Eigen::MatrixXd solution = sparse_qr_solver.solve(rhs_);

  // Check that the system could be solved successfully.
  if (sparse_qr_solver.info() != Eigen::Success) {
    VLOG(2) << "The sparse linear system could not be solved.";
    return false;
  }

  // Project all solutions into a valid SO3 rotation space. The linear system
  // above makes no constraint on the space of the solutions, so the final
  // solutions are not guaranteed to be valid rotations (e.g., det(R) may not be
  // +1).
  global_orientations->reserve(view_id_map_.size());
  for (const auto& view_id_map : view_id_map_) {
    // If this is the view held constant, set it to the identity.
    if (view_id_map.second == -1) {
      global_orientations->emplace(view_id_map.first, Eigen::Vector3d::Zero());
      continue;
    }

    const Matrix3d non_so3_rotation =
        solution.block<3, 3>(3 * view_id_map.second, 0);
    const Matrix3d rotation = ProjectToRotationMatrix(non_so3_rotation);

    // Convert to angle axis.
    Eigen::Vector3d rotation_aa;
    ceres::RotationMatrixToAngleAxis(
        ceres::ColumnMajorAdapter3x3(rotation.data()), rotation_aa.data());
    global_orientations->emplace(view_id_map.first, rotation_aa);
  }

  return true;
}

}  // namespace theia
