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

#include "theia/sfm/pose/fundamental_matrix_util.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <glog/logging.h>

#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::DiagonalMatrix;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector3d;

// Given a fundmental matrix, decompose the fundmental matrix and recover focal
// lengths f1, f2 >0 such that diag([f2 f2 1]) F diag[f1 f1 1]) is a valid
// essential matrix. This assumes a principal point of (0, 0) for both cameras.
// Returns true on success, false otherwise.
bool FocalLengthsFromFundamentalMatrix(const double fmatrix[3 * 3],
                                       double* focal_length1,
                                       double* focal_length2) {
  Map<const Matrix3d> fundamental_matrix(fmatrix);

  // Compute the epipoles for each image.
  const Vector3d epipole1 = fundamental_matrix.jacobiSvd(Eigen::ComputeFullV)
                                .matrixV()
                                .rightCols<1>();
  const Vector3d epipole2 = fundamental_matrix.transpose()
                                .jacobiSvd(Eigen::ComputeFullV)
                                .matrixV()
                                .rightCols<1>();
  if (epipole1.x() == 0 || epipole2.x() == 0) {
    VLOG(3) << "Optical axes are collinear. Cannot recover the focal length.";
    return false;
  }

  // Find the rotation that takes epipole1 to (e_0, 0, e_2) and the
  // epipole2 to (e'_0, 0, e'_2). If we form a rotation matrix:
  // R = [ cos x  -sin x  0 ]
  //     [ sin x   cos x  0 ]
  //     [ 0       0      1 ]
  // then we can solve for the angle x such that R * e1 = (e_1, 0, e_3).
  // We can solve this simply be investigating the second row and noting that
  // e1(0) * sin x + e2 * cos x = 0.
  const double theta1 = atan2(-epipole1(1), epipole1(0));
  const double theta2 = atan2(-epipole2(1), epipole2(0));

  Matrix3d rotation1, rotation2;
  rotation1 <<
      cos(theta1), -sin(theta1), 0,
      sin(theta1), cos(theta1), 0,
      0, 0, 1;
  rotation2 <<
      cos(theta2), -sin(theta2), 0,
      sin(theta2), cos(theta2), 0,
      0, 0, 1;

  const Matrix3d rotated_fmatrix =
      rotation2 * fundamental_matrix * rotation1.transpose();
  // With the normalized epipoles, the fundamental matrix is now of the form:
  // F = [ e'_2   0    0   ] [ a b a ] [ e_2   0     0  ]
  //     [ 0      1    0   ] [ c d c ] [ 0     1     0  ]
  //     [ 0      0  -e'_1 ] [ a b a ] [ 0     0   -e_1 ]
  const Vector3d rotated_epipole1 = rotation1 * epipole1;
  const Vector3d rotated_epipole2 = rotation2 * epipole2;

  Matrix3d factorized_matrix =
      DiagonalMatrix<double, 3>(rotated_epipole2(2), 1, -rotated_epipole2(0))
          .inverse() *
      rotated_fmatrix *
      DiagonalMatrix<double, 3>(rotated_epipole1(2), 1, -rotated_epipole1(0))
          .inverse();

  // For convenience, as defined above.
  const double a = factorized_matrix(0, 0);
  const double b = factorized_matrix(0, 1);
  const double c = factorized_matrix(1, 0);
  const double d = factorized_matrix(1, 1);

  const double focal_length1_sq =
      (-a * c * rotated_epipole1(0) * rotated_epipole1(0)) /
      (a * c * rotated_epipole1(2) * rotated_epipole1(2) + b * d);
  const double focal_length2_sq =
      (-a * b * rotated_epipole2(0) * rotated_epipole2(0)) /
      (a * b * rotated_epipole2(2) * rotated_epipole2(2) + c * d);

  if (focal_length1_sq < 0 || focal_length2_sq < 0) {
    VLOG(3) << "Real focal length values could not be extracted.";
    return false;
  }

  *focal_length1 = sqrt(focal_length1_sq);
  *focal_length2 = sqrt(focal_length2_sq);
  return true;
}


void ProjectionMatricesFromFundamentalMatrix(const double fmatrix[3 * 3],
                                             double pmatrix1[3 * 4],
                                             double pmatrix2[3 * 4]) {
  Map<const Matrix3d> fmatrix_map(fmatrix);
  const Vector3d right_epipole = fmatrix_map.jacobiSvd(Eigen::ComputeFullV)
      .matrixV().rightCols<1>();

  Map<Matrix3d>(pmatrix1, 3, 3) = Matrix3d::Identity();
  Map<Vector3d>(pmatrix1 + 9, 3) = Vector3d::Zero();

  Map<Matrix3d>(pmatrix2, 3, 3) =
      CrossProductMatrix(right_epipole) * fmatrix_map.transpose();
  Map<Vector3d>(pmatrix2 + 9, 3) = right_epipole;
}

// Ported from Hartley and Zisserman:
// http://www.robots.ox.ac.uk/~vgg/hzbook/code/vgg_multiview/vgg_F_from_P.m
void FundamentalMatrixFromProjectionMatrices(const double pmatrix1[3 * 4],
                                             const double pmatrix2[3 * 4],
                                             double fmatrix[3 * 3]) {
  Map<const Matrix<double, 3, 4> > projection1(pmatrix1);
  Map<const Matrix<double, 3, 4> > projection2(pmatrix2);
  Map<Matrix3d> fundamental_matrix(fmatrix);

  const int index1[3] = {1, 2, 0};
  const int index2[3] = {2, 0, 1};
  Eigen::Matrix4d temp_mat;
  for (int r = 0; r < 3; r++) {
    temp_mat.row(2) = projection1.row(index1[r]);
    temp_mat.row(3) = projection1.row(index2[r]);
    for (int c = 0; c < 3; c++) {
      temp_mat.row(0) = projection2.row(index1[c]);
      temp_mat.row(1) = projection2.row(index2[c]);
      fundamental_matrix(r, c) = temp_mat.determinant();
    }
  }
}

void EssentialMatrixFromFundamentalMatrix(const double fmatrix[3 * 3],
                                          const double focal_length1,
                                          const double focal_length2,
                                          double ematrix[3 * 3]) {
  const Eigen::Map<const Matrix3d> fundamental_matrix(fmatrix);
  Eigen::Map<Matrix3d> essential_matrix(ematrix);
  essential_matrix =
      Eigen::DiagonalMatrix<double, 3>(focal_length2, focal_length2, 1.0) *
      fundamental_matrix *
      Eigen::DiagonalMatrix<double, 3>(focal_length1, focal_length1, 1.0);
}

}  // namespace theia
