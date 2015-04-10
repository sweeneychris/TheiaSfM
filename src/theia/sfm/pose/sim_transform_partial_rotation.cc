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

#include "theia/sfm/pose/sim_transform_partial_rotation.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <glog/logging.h>
#include <math.h>
#include <limits>

#include "theia/alignment/alignment.h"
#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::EigenSolver;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 5, 1> Vector5d;

namespace {

bool SolveQEP(const Matrix5d& M, const Matrix5d& C, const Matrix5d& K,
              std::vector<double>* eigenvalues,
              std::vector<Vector5d>* eigenvectors) {
  // Solves the quadratic eigenvalue problem:
  //
  //   Q(s)x = 0
  //
  // where:
  //
  //   Q(s) = s^2*M + s*C + K
  //
  // Returns true if the problem could be solved, false otherwise.
  //
  // This is converted to a generalized eigenvalue as described in:
  //   http://en.wikipedia.org/wiki/Quadratic_eigenvalue_problem
  //
  // The generalized eigenvalue problem is in this form:
  //
  // [ C  K ]z = s[ -M  0  ]z
  // [-I  0 ]     [ 0   -I ]
  //
  // With eigenvector z = [ sx ] and s is the eigenvalue.
  //                      [  x ]
  //
  // The eigenvector of the quadratic eigenvalue problem can be extracted from
  // z.
  //
  // Where z is the eigenvector, s is the eigenvalue.
  // This generalized eigenvalue can be converted to a standard eigenvalue
  // problem by multiplying by the inverse of the RHS matrix. The inverse of
  // the RHS matrix is particularly simple:
  //
  // [ -inv(M)  0 ]
  // [ 0       -I ]
  //
  // So the generalized eigenvalue problem reduces to a standard eigenvalue
  // problem on the constraint matrix:
  //
  // [ -inv(M)C  -inv(M)K ]z = sz
  // [ I         0        ]

  const Eigen::FullPivLU<Matrix5d> lu(M);
  if (!lu.isInvertible()) {
    VLOG(2) << "Matrix M is not invertible (rank " << lu.rank() << ")";
    return false;
  }

  // Negate inverse of M.
  const Matrix5d inv_M = -1.0 * lu.inverse();

  // Set up constraint matrix.
  Matrix<double, 10, 10> constraint = Matrix<double, 10, 10>::Zero();
  // Set upper-left to - inv(M) * C.
  constraint.block<5, 5>(0, 0) = inv_M * C;
  // Set upper-right to -inv(M) * K.
  constraint.block<5, 5>(0, 5) = inv_M * K;
  // Set lower-left to identity.
  constraint.block<5, 5>(5, 0) = Matrix5d::Identity();

  // Extract the left eigenvectors and values from the constraint matrix.
  EigenSolver<Matrix<double, 10, 10> > eig_solver(constraint);
  for (int i = 0; i < eig_solver.eigenvalues().size(); i++) {
    // Two solutions that always occur correspond to (s^2 + 1) so skip these
    // solutions.
    if (fabs(eig_solver.eigenvalues()[i].imag()) == 1 ||
        (i > 0 && eig_solver.eigenvalues()[i].real() ==
         eig_solver.eigenvalues()[i - 1].real())) {
      continue;
    }
    eigenvalues->push_back(eig_solver.eigenvalues()[i].real());
    eigenvectors->push_back(
        Vector5d(eig_solver.eigenvectors().col(i).tail<5>().real()));
  }
  return true;
}

}  // namespace

void SimTransformPartialRotation(const Vector3d& axis,
                                 const Vector3d image_one_ray_directions[5],
                                 const Vector3d image_one_ray_origins[5],
                                 const Vector3d image_two_ray_directions[5],
                                 const Vector3d image_two_ray_origins[5],
                                 std::vector<Quaterniond>* soln_rotations,
                                 std::vector<Vector3d>* soln_translations,
                                 std::vector<double>* soln_scales) {
  CHECK_DOUBLE_EQ(axis.squaredNorm(), 1.0);
  CHECK_NOTNULL(soln_rotations)->clear();
  CHECK_NOTNULL(soln_translations)->clear();

  // The generalized epipolar constraint between two sets of rays in different
  // coordinate systems (two generalized cameras).
  //
  // The rays are converted to plucker coordinates:
  //   L = [ q  ]
  //       [ p  ]
  //
  // Assuming the first coordinate system has zero rotation and translation,
  // the generalized epipolar constraint is:
  //
  //   q2' * [t]x R * q1 + q2' * R * p1 + q1' * R * p2 = 0
  //
  // v' means transpose v, t is unknown translation, [t]x is the matrix
  // cross product form.
  //
  // Re-arranging to isolate t:
  //
  //   -(q2 x (R * q1))' * t + q2' * R * p1 + q1' * R' * p2 = 0
  //
  // In vector form:
  //
  //  [-(q2 x (R * q1))', q2' * R * p1 + q1' * R' * p2 ] * [ t ] = 0
  //                                                       [ 1 ]
  //
  // If we consider that L2 may not have the same scale as L1, then we must
  // additionally reconcile a scale s. Our equation then becomes:
  //
  //  [-(q2 x (R * q1))',  q1' * R' * p2, q2' * R * p1 ] * [ t ] = 0
  //                                                       [ s ]
  //                                                       [ 1 ]
  //
  // Each correspondence gives another constraint and the constraints can
  // be stacked to create a constraint matrix.
  //
  // The rotation matrix R can be parameterized, up to scale factor, as:
  //
  //   R ~ 2 * (v * v' + alpha * [v]x) + (alpha^2 - 1)I
  //
  //   I = Identity matrix.
  //
  // where is v is the known (unit length) axis and alpha is related to the
  // unknown angle of rotation.
  //
  // The epipolar constraint holds if the rotation matrix is scaled, so the
  // rotation matrix parameterization above can be used. Each row of the
  // constraint matrix is thus a quadratic function in alpha, so the constraint
  // can be written as:
  //
  //   [ M * alpha^2 + C * alpha + k ] * [ t ] = 0
  //                                     [ s ]
  //                                     [ 1 ]
  //
  // This is standard quadratic eigenvalue problem (QEP) and is solved using
  // standard methods.

  // Creates the matrices for the QEP problem.
  Matrix5d M;
  Matrix5d C;
  Matrix5d K;

  // Based on the rotation parameterization above,
  //   R = alpha^2 * rot_alpha_sq + alpha * rot_alpha + rot_constant
  const Matrix3d rot_alpha_sq = Matrix3d::Identity();
  const Matrix3d rot_alpha = 2.0 * CrossProductMatrix(axis);
  const Matrix3d rot_constant =
      2 * axis * axis.transpose() - Matrix3d::Identity();

  for (int i = 0; i < 5; ++i) {
    const Vector3d& f1(image_one_ray_directions[i]);
    const Vector3d& f2(image_two_ray_directions[i]);
    const Matrix3d pos1 = CrossProductMatrix(image_one_ray_origins[i]);
    const Matrix3d pos2 = CrossProductMatrix(image_two_ray_origins[i]);

    M.row(i).head(3) = f1.cross(rot_alpha_sq * f2);
    M.row(i)[3] = -f1.transpose() * (rot_alpha_sq * pos2) * f2;
    M.row(i)[4] = f1.transpose() * (pos1 * rot_alpha_sq) * f2;

    C.row(i).head(3) = f1.cross(rot_alpha * f2);
    C.row(i)[3] = -f1.transpose() * (rot_alpha * pos2) * f2;
    C.row(i)[4] = f1.transpose() * (pos1 * rot_alpha) * f2;

    K.row(i).head(3) = f1.cross(rot_constant * f2);
    K.row(i)[3] = -f1.transpose() * (rot_constant * pos2) * f2;
    K.row(i)[4] = f1.transpose() * (pos1 * rot_constant) * f2;
  }

  std::vector<double> eigenvalues;
  std::vector<Vector5d> eigenvectors;
  static const double kWTolerance = 1e-12;
  if (SolveQEP(M, C, K, &eigenvalues, &eigenvectors)) {
    // Extracts the translations and rotations from the eigenvalues and
    // eigenvectors of the QEP problem.
    for (int i = 0; i < eigenvalues.size(); ++i) {
      if (fabs(eigenvectors[i][4]) < kWTolerance ||
          eigenvectors[i][3] / eigenvectors[i][4] < 0) {
      continue;
    }

      Quaterniond quat(eigenvalues[i], axis[0], axis[1], axis[2]);
      quat.normalize();
      soln_rotations->push_back(quat);
      const Vector3d rotated_translation =
          Vector3d(eigenvectors[i][0] / eigenvectors[i][4],
                   eigenvectors[i][1] / eigenvectors[i][4],
                   eigenvectors[i][2] / eigenvectors[i][4]);
      soln_translations->push_back(rotated_translation);
      soln_scales->push_back(eigenvectors[i][3] / eigenvectors[i][4]);
    }
  } else {
    // When there is zero rotation the vector part of the quaternion disappears
    // and it becomes ([0], 1) (where [0] is the zero vector) and the SolveQEP
    // method cannot be used.
    // However from the equations for M, C, and K above we can see that the
    // C and K matrices contain the axis, which is 0 assuming zero rotation,
    // so to solve for the translation we can directly extract the null space
    // of M.
    const Eigen::FullPivLU<Matrix5d> lu(M);
    if (lu.dimensionOfKernel() != 1) {
      LOG(WARNING) << "The input is a degenerate configuration resulting in a "
                      "rank deficient matrix!";
      return;
    }

    const Vector5d kernel = lu.kernel();
    if (fabs(kernel[4]) > kWTolerance) {
      soln_rotations->push_back(Quaterniond::Identity());
      soln_translations->push_back(kernel.head<3>() / kernel[4]);
      soln_scales->push_back(kernel[3] / kernel[4]);
    }
  }
}

}  // namespace theia
