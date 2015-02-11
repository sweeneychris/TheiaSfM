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

#include "theia/sfm/pose/five_point_relative_pose.h"

#include <Eigen/Dense>
#include <glog/logging.h>

#include <cmath>
#include <ctime>
#include <vector>

#include "theia/math/polynomial.h"
#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Matrix;
using Eigen::RowVector3d;
using Eigen::RowVector4d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::VectorXd;

typedef Matrix<double, 10, 10> Matrix10d;

namespace {

// Multiply two degree one polynomials of variables x, y, z.
// E.g. p1 = a[0]x + a[1]y + a[2]z + a[3]
// Output order: x^2 xy y^2 xz yz z^2 x y z 1 (GrevLex)
Matrix<double, 1, 10> MultiplyDegOnePoly(const RowVector4d& a,
                                         const RowVector4d& b) {
  Matrix<double, 1, 10> output;
  // x^2
  output(0) = a(0) * b(0);
  // xy
  output(1) = a(0) * b(1) + a(1) * b(0);
  // y^2
  output(2) = a(1) * b(1);
  // xz
  output(3) = a(0) * b(2) + a(2) * b(0);
  // yz
  output(4) = a(1) * b(2) + a(2) * b(1);
  // z^2
  output(5) = a(2) * b(2);
  // x
  output(6) = a(0) * b(3) + a(3) * b(0);
  // y
  output(7) = a(1) * b(3) + a(3) * b(1);
  // z
  output(8) = a(2) * b(3) + a(3) * b(2);
  // 1
  output(9) = a(3) * b(3);
  return output;
}

// Multiply a 2 deg poly (in x, y, z) and a one deg poly in GrevLex order.
// x^3 x^2y xy^2 y^3 x^2z xyz y^2z xz^2 yz^2 z^3 x^2 xy y^2 xz yz z^2 x y z 1
Matrix<double, 1, 20> MultiplyDegTwoDegOnePoly(const Matrix<double, 1, 10>& a,
                                               const RowVector4d& b) {
  Matrix<double, 1, 20> output;
  // x^3
  output(0) = a(0) * b(0);
  // x^2y
  output(1) = a(0) * b(1) + a(1) * b(0);
  // xy^2
  output(2) = a(1) * b(1) + a(2) * b(0);
  // y^3
  output(3) = a(2) * b(1);
  // x^2z
  output(4) = a(0) * b(2) + a(3) * b(0);
  // xyz
  output(5) = a(1) * b(2) + a(3) * b(1) + a(4) * b(0);
  // y^2z
  output(6) = a(2) * b(2) + a(4) * b(1);
  // xz^2
  output(7) = a(3) * b(2) + a(5) * b(0);
  // yz^2
  output(8) = a(4) * b(2) + a(5) * b(1);
  // z^3
  output(9) = a(5) * b(2);
  // x^2
  output(10) = a(0) * b(3) + a(6) * b(0);
  // xy
  output(11) = a(1) * b(3) + a(6) * b(1) + a(7) * b(0);
  // y^2
  output(12) = a(2) * b(3) + a(7) * b(1);
  // xz
  output(13) = a(3) * b(3) + a(6) * b(2) + a(8) * b(0);
  // yz
  output(14) = a(4) * b(3) + a(7) * b(2) + a(8) * b(1);
  // z^2
  output(15) = a(5) * b(3) + a(8) * b(2);
  // x
  output(16) = a(6) * b(3) + a(9) * b(0);
  // y
  output(17) = a(7) * b(3) + a(9) * b(1);
  // z
  output(18) = a(8) * b(3) + a(9) * b(2);
  // 1
  output(19) = a(9) * b(3);
  return output;
}

Matrix<double, 1, 20> GetDeterminantConstraint(
    const Matrix<double, 1, 4> null_space[3][3]) {
  // Singularity constraint.
  const Matrix<double, 1, 20> determinant =
      MultiplyDegTwoDegOnePoly(
          MultiplyDegOnePoly(null_space[0][1], null_space[1][2]) -
          MultiplyDegOnePoly(null_space[0][2], null_space[1][1]),
          null_space[2][0]) +
      MultiplyDegTwoDegOnePoly(
          MultiplyDegOnePoly(null_space[0][2], null_space[1][0]) -
          MultiplyDegOnePoly(null_space[0][0], null_space[1][2]),
          null_space[2][1]) +
      MultiplyDegTwoDegOnePoly(
          MultiplyDegOnePoly(null_space[0][0], null_space[1][1]) -
          MultiplyDegOnePoly(null_space[0][1], null_space[1][0]),
          null_space[2][2]);
  return determinant;
}

// Shorthand for multiplying the Essential matrix with its transpose.
Matrix<double, 1, 10> EETranspose(
    const Matrix<double, 1, 4> null_space[3][3], int i, int j) {
  return MultiplyDegOnePoly(null_space[i][0], null_space[j][0]) +
      MultiplyDegOnePoly(null_space[i][1], null_space[j][1]) +
      MultiplyDegOnePoly(null_space[i][2], null_space[j][2]);
}


// Builds the trace constraint: EEtE - 1/2 trace(EEt)E = 0
Matrix<double, 9, 20> GetTraceConstraint(
    const Matrix<double, 1, 4> null_space[3][3]) {
  Matrix<double, 9, 20> trace_constraint;

  // Comput EEt.
  Matrix<double, 1, 10> eet[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      eet[i][j] = 2 * EETranspose(null_space, i, j);
    }
  }

  // Compute the trace.
  const Matrix<double, 1, 10> trace = eet[0][0] + eet[1][1] + eet[2][2];

  // Multiply EEt with E.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      trace_constraint.row(3 * i + j) =
          MultiplyDegTwoDegOnePoly(eet[i][0], null_space[0][j]) +
          MultiplyDegTwoDegOnePoly(eet[i][1], null_space[1][j]) +
          MultiplyDegTwoDegOnePoly(eet[i][2], null_space[2][j]) -
          0.5 * MultiplyDegTwoDegOnePoly(trace, null_space[i][j]);
    }
  }

  return trace_constraint;
}

Matrix<double, 10, 20> BuildConstraintMatrix(
    const Matrix<double, 1, 4> null_space[3][3]) {
  Matrix<double, 10, 20> constraint_matrix;
  constraint_matrix.block<9, 20>(0, 0) = GetTraceConstraint(null_space);
  constraint_matrix.row(9) = GetDeterminantConstraint(null_space);
  return constraint_matrix;
}

}  // namespace

// Implementation of Nister from "An Efficient Solution to the Five-Point
// Relative Pose Problem"
bool FivePointRelativePose(const Vector2d image1_points[5],
                           const Vector2d image2_points[5],
                           std::vector<Matrix3d>* essential_matrices) {
  // Step 1. Create the 5x9 matrix containing epipolar constraints.
  //   Essential matrix is a linear combination of the 4 vectors spanning the
  //   null space of this matrix.
  Matrix<double, 5, 9> epipolar_constraint;
  for (int i = 0; i < 5; i++) {
    // Fill matrix with the epipolar constraint from q'_t*E*q = 0. Where q is
    // from the first image, and q' is from the second.
    epipolar_constraint.row(i) <<
        image2_points[i].x() * image1_points[i].x(),
        image2_points[i].y() * image1_points[i].x(),
        image1_points[i].x(),
        image2_points[i].x() * image1_points[i].y(),
        image2_points[i].y() * image1_points[i].y(),
        image1_points[i].y(),
        image2_points[i].x(),
        image2_points[i].y(),
        1.0;
  }

  const Eigen::FullPivLU<Matrix<double, 5, 9> > lu(epipolar_constraint);
  if (lu.dimensionOfKernel() != 4) {
    return false;
  }
  const Matrix<double, 9, 4>& null_space = lu.kernel();

  const Matrix<double, 1, 4> null_space_matrix[3][3] = {
    { null_space.row(0), null_space.row(3), null_space.row(6) },
    { null_space.row(1), null_space.row(4), null_space.row(7) },
    { null_space.row(2), null_space.row(5), null_space.row(8) }
  };

  // Step 2. Expansion of the epipolar constraints on the determinant and trace.
  const Matrix<double, 10, 20> constraint_matrix =
      BuildConstraintMatrix(null_space_matrix);

  // Step 3. Eliminate part of the matrix to isolate polynomials in z.
  Eigen::FullPivLU<Matrix10d> c_lu(constraint_matrix.block<10, 10>(0, 0));
  Matrix10d eliminated_matrix =
      c_lu.solve(constraint_matrix.block<10, 10>(0, 10));

  Matrix10d action_matrix = Matrix10d::Zero();
  action_matrix.block<3, 10>(0, 0) = eliminated_matrix.block<3, 10>(0, 0);
  action_matrix.row(3) = eliminated_matrix.row(4);
  action_matrix.row(4) = eliminated_matrix.row(5);
  action_matrix.row(5) = eliminated_matrix.row(7);
  action_matrix(6, 0) = -1.0;
  action_matrix(7, 1) = -1.0;
  action_matrix(8, 3) = -1.0;
  action_matrix(9, 6) = -1.0;

  Eigen::EigenSolver<Matrix10d> eigensolver(action_matrix);
  const auto& eigenvectors = eigensolver.eigenvectors();
  const auto& eigenvalues = eigensolver.eigenvalues();

  // Now that we have x, y, and z we need to substitute them back into the null
  // space to get a valid essential matrix solution.
  for (int i = 0; i < 10; i++) {
    // Only consider real solutions.
    if (eigenvalues(i).imag() != 0) {
      continue;
    }
    Matrix3d ematrix;
    Map<Matrix<double, 9, 1> >(ematrix.data()) =
        null_space * eigenvectors.col(i).tail<4>().real();
    essential_matrices->emplace_back(ematrix);
  }

  return essential_matrices->size() > 0;
}

}  // namespace theia
