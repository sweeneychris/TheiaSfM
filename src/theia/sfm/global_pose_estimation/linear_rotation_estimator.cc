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
#include <vector>

#include "spectra/include/SymEigsSolver.h"

#include "theia/math/matrix/linear_operator.h"
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
      triplet_list->emplace_back(row + i, col + j, mat(i, j));
    }
  }
}

}  // namespace

// Estimates the global orientations of all views based on an initial
// guess. Returns true on successful estimation and false otherwise.
bool LinearRotationEstimator::EstimateRotations(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations) {
  for (const auto& view_pair : view_pairs) {
    AddRelativeRotationConstraint(view_pair.first, view_pair.second.rotation_2);
  }
  return EstimateRotations(global_orientations);
}

// Add the relative rotation matrix constraint that minimizes the Frobenius norm
// of the matrices. That is, set up the linear system such that
// ||R_j - R_ij * R_i|| is minimized
void LinearRotationEstimator::AddRelativeRotationConstraint(
    const ViewIdPair& view_id_pair, const Eigen::Vector3d& relative_rotation) {
  static const int kNumRotationMatrixDimensions = 3;

  // Add a new view-id to sparse matrix index mapping. This is a no-op if the
  // view has already been assigned to a matrix index.
  InsertIfNotPresent(&view_id_map_, view_id_pair.first, view_id_map_.size());
  InsertIfNotPresent(&view_id_map_, view_id_pair.second, view_id_map_.size());
  const int view1_index = view_id_map_[view_id_pair.first];
  const int view2_index = view_id_map_[view_id_pair.second];

  // Add the corresponding entries to A^t * A. Note that for each row, we may
  // represent the row as a block matrix: [B | C]. Thus, the corresponding
  // entries for A^t * A are:
  //
  //   [B | C]^t * [B | C] = [B^t * B | B^t * C]
  //                         [C^t * B | C^t * C]
  //
  // Note that for rotation matrices, B^t * B = C^t * B = I. Thus we have:
  //
  //   [B | C]^t * [B | C] = [   I    | B^t * C]
  //                         [C^t * B |    I   ]
  //
  // We only store the upper triangular portion of the matrix since it is
  // symmetric. In our case, one of B or C will be the identity, which further
  // simplifies the matrix entries.
  //
  // First, add the identity entries along the diagonal.
  for (int i = 0; i < 3; i++) {
    constraint_entries_.emplace_back(
        kNumRotationMatrixDimensions * view1_index + i,
        kNumRotationMatrixDimensions * view1_index + i,
        1.0);
    constraint_entries_.emplace_back(
        kNumRotationMatrixDimensions * view2_index + i,
        kNumRotationMatrixDimensions * view2_index + i,
        1.0);
  }

  // Add the 3x3 matrix B^t * C. This corresponds either to -R_ij or -R_ij^t
  // depending on the order of the view indices.
  Eigen::Matrix3d relative_rotation_matrix;
  ceres::AngleAxisToRotationMatrix(
      relative_rotation.data(),
      ceres::ColumnMajorAdapter3x3(relative_rotation_matrix.data()));

  // We seek to insert B^t * C into the linear system as noted above. If the
  // view1 index comes before the view2 index, then B = -R_ij and C =
  // I. Otherwise, B = I and C = -R_ij.
  if (view1_index < view2_index) {
    Fill3x3SparseMatrix(
        -relative_rotation_matrix.transpose(),
        kNumRotationMatrixDimensions * view1_index,
        kNumRotationMatrixDimensions * view2_index,
        &constraint_entries_);
  } else {
    Fill3x3SparseMatrix(
        -relative_rotation_matrix,
        kNumRotationMatrixDimensions * view2_index,
        kNumRotationMatrixDimensions * view1_index,
        &constraint_entries_);
  }
}

// Given the relative rotation constraints added with
// AddRelativeRotationConstraint, this method returns the robust estimation of
// global camera orientations. Like the method above, this requires an initial
// estimate of the global orientations.
bool LinearRotationEstimator::EstimateRotations(
    std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations) {
  CHECK_GT(constraint_entries_.size(), 0);
  CHECK_NOTNULL(global_orientations);
  static const int kNumRotationMatrixDimensions = 3;

  // Setup the sparse linear system.
  Eigen::SparseMatrix<double> constraint_matrix(
      view_id_map_.size() * kNumRotationMatrixDimensions,
      view_id_map_.size() * kNumRotationMatrixDimensions);
  constraint_matrix.setFromTriplets(constraint_entries_.begin(),
                                    constraint_entries_.end());

  // Compute the 3 eigenvectors corresponding to the smallest eigenvalues. These
  // orthogonal vectors will contain the solution rotation matrices.
  SparseSymShiftSolveLLT op(constraint_matrix);
  Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN,
                              SparseSymShiftSolveLLT> eigs(&op, 3, 6, 0.0);
  eigs.init();
  eigs.compute();

  // The solution appears in the first three eigenvectors.
  const Eigen::MatrixXd solution =
      eigs.eigenvectors().leftCols<kNumRotationMatrixDimensions>();

  // Project all solutions into a valid SO3 rotation space. The linear system
  // above makes no constraint on the space of the solutions, so the final
  // solutions are not guaranteed to be valid rotations (e.g., det(R) may not be
  // +1).
  global_orientations->reserve(view_id_map_.size());
  for (const auto& view_id_map : view_id_map_) {
    const Matrix3d non_so3_rotation =
        solution
            .block<kNumRotationMatrixDimensions, kNumRotationMatrixDimensions>(
                kNumRotationMatrixDimensions * view_id_map.second, 0);
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
