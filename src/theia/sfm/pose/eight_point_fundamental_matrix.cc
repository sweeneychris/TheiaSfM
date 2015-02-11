// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#include "theia/sfm/pose/eight_point_fundamental_matrix.h"

#include <glog/logging.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/LU>

#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::JacobiSVD;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector3d;

bool NormalizedEightPointFundamentalMatrix(
    const std::vector<Vector2d>& image_1_points,
    const std::vector<Vector2d>& image_2_points,
    Matrix3d* fundamental_matrix) {
  CHECK_EQ(image_1_points.size(), image_2_points.size());
  CHECK_GE(image_1_points.size(), 8);

  std::vector<Vector2d> norm_img1_points(image_1_points.size());
  std::vector<Vector2d> norm_img2_points(image_2_points.size());

  // Normalize the image points.
  Matrix3d img1_norm_mat, img2_norm_mat;
  NormalizeImagePoints(image_1_points, &norm_img1_points, &img1_norm_mat);
  NormalizeImagePoints(image_2_points, &norm_img2_points, &img2_norm_mat);

  // Build the constraint matrix based on x2' * F * x1 = 0.
  Matrix<double, Eigen::Dynamic, 9> constraint_matrix(image_1_points.size(), 9);
  for (int i = 0; i < image_1_points.size(); i++) {
    constraint_matrix.block<1, 3>(i, 0) = norm_img1_points[i].homogeneous();
    constraint_matrix.block<1, 3>(i, 0) *= norm_img2_points[i].x();
    constraint_matrix.block<1, 3>(i, 3) = norm_img1_points[i].homogeneous();
    constraint_matrix.block<1, 3>(i, 3) *= norm_img2_points[i].y();
    constraint_matrix.block<1, 3>(i, 6) = norm_img1_points[i].homogeneous();
  }

  // Solve the constraint equation for F from nullspace extraction.
  // An LU decomposition is efficient for the minimally constrained case.
  // Otherwise, use an SVD.
  Matrix<double, 9, 1> normalized_fvector;
  if (image_1_points.size() == 8) {
    const auto lu_decomposition = constraint_matrix.fullPivLu();
    if (lu_decomposition.dimensionOfKernel() != 1) {
      return false;
    }
    normalized_fvector = lu_decomposition.kernel();
  } else {
    JacobiSVD<Matrix<double, Eigen::Dynamic, 9> > cmatrix_svd(
       constraint_matrix, Eigen::ComputeFullV);
    normalized_fvector = cmatrix_svd.matrixV().col(8);
  }

  // NOTE: This is the transpose of a valid fundamental matrix! We implement a
  // "lazy" transpose and defer it to the SVD a few lines below.
  Eigen::Map<const Matrix3d> normalized_fmatrix(normalized_fvector.data());

  // Find the closest singular matrix to F under frobenius norm. We can compute
  // this matrix with SVD.
  JacobiSVD<Matrix3d> fmatrix_svd(normalized_fmatrix.transpose(),
                                  Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vector3d singular_values = fmatrix_svd.singularValues();
  singular_values[2] = 0.0;
  *fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                        fmatrix_svd.matrixV().transpose();

  // Correct for the point normalization.
  *fundamental_matrix =
      img2_norm_mat.transpose() * (*fundamental_matrix) * img1_norm_mat;

  return true;
}

}  // namespace theia
