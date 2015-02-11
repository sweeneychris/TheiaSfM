// Copyright (C) 2014 The Regents of the University of California (Regents)
// and Google, Inc. All rights reserved.
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
//     * Neither the name of The Regents or University of California, Google,
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu), John Flynn (jflynn@google.com)

#include "theia/sfm/pose/three_point_relative_pose_partial_rotation.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <math.h>

#include <limits>

namespace theia {

using Eigen::AngleAxisd;
using Eigen::EigenSolver;
using Eigen::JacobiSVD;
using Eigen::Map;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Quaterniond;
using Eigen::Vector3d;

namespace {

bool SolveQEP(const Matrix3d& M, const Matrix3d& C, const Matrix3d& K,
              std::vector<double>* eigenvalues,
              std::vector<Vector3d>* eigenvectors) {
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

  Matrix3d inv_M;
  static const double kDeterminantThreshold = 1e-12;
  bool invert_success;
  // Check that determinant of M is larger than threshold. This threshold only
  // seems to be reached when there is no rotation. TODO(jflynn): verify.
  M.computeInverseWithCheck(inv_M, invert_success, kDeterminantThreshold);
  if (!invert_success) {
    return false;
  }

  // Negate inverse of M.
  inv_M = -1.0 * inv_M;

  // Set up constraint matrix.
  Matrix<double, 6, 6> constraint = Matrix<double, 6, 6>::Zero();
  // Set upper-left to - inv(M) * C.
  constraint.block<3, 3>(0, 0) = inv_M * C;
  // Set upper-right to -inv(M) * K.
  constraint.block<3, 3>(0, 3) = inv_M * K;
  // Set lower-left to identity.
  constraint.block<3, 3>(3, 0) = Matrix3d::Identity();
  // Extract the left eigenvectors and values from the constraint matrix.
  EigenSolver<Matrix<double, 6, 6> > eig_solver(constraint);
  const double kImagEigenValueTolerance = 1e-12;

  for (int i = 0; i < eig_solver.eigenvalues().size(); i++) {
    // Ignore roots corresponding to s^2 + 1.
    if (fabs(eig_solver.eigenvalues()[i].imag() - 1) <
            kImagEigenValueTolerance ||
        fabs(eig_solver.eigenvalues()[i].imag() + 1) <
            kImagEigenValueTolerance) {
      continue;
    }
    // Only consider the real eigenvalues and corresponding eigenvectors.
    eigenvalues->push_back(eig_solver.eigenvalues()[i].real());
    eigenvectors->push_back(
        Vector3d(eig_solver.eigenvectors().col(i).tail<3>().real()));
  }
  return true;
}

}  // namespace

void ThreePointRelativePosePartialRotation(
    const Vector3d& axis,
    const Vector3d image_1_rays[3],
    const Vector3d image_2_rays[3],
    std::vector<Quaterniond>* soln_rotations,
    std::vector<Vector3d>* soln_translations) {
  CHECK_DOUBLE_EQ(axis.squaredNorm(), 1.0);
  CHECK_NOTNULL(soln_rotations)->clear();
  CHECK_NOTNULL(soln_translations)->clear();

  // Each correspondence gives another constraint and the constraints can
  // be stacked to create a constraint matrix.
  //
  // The rotation matrix R can be parameterized, up to scale factor, as:
  //
  //   R ~ 2 * (v * v' + s[v]x) + (s^2 - 1)I
  //
  //   I = Identity matrix.
  //
  // where is v is the known (unit length) axis and s is related to the
  // unknown angle of rotation.
  //
  // The epipolar constraint holds if the rotation matrix is scaled, so the
  // rotation matrix parameterization above can be used. Each row of the
  // constraint matrix is thus a function of s^2 and s, so the constraint can be
  // written as:
  //
  //   [ M * s^2 + C *s + k ] * [ t ] = 0
  //                            [ 1 ]
  //
  // This is standard quadratic eigenvalue problem (QEP) and is solved using
  // standard methods.

  // Creates the matrices for the QEP problem.
  Matrix3d M;
  Matrix3d C;
  Matrix3d K;
  for (int i = 0; i < 3; ++i) {
    const Vector3d& q1(image_1_rays[i]);
    const Vector3d& q2(image_2_rays[i]);
    M.row(i) = q2.cross(q1);
    C.row(i) = 2.0 * q2.cross(axis.cross(q1));
    K.row(i) = 2.0 * q1.dot(axis) * q2.cross(axis) - q2.cross(q1);
  }

  std::vector<double> eigenvalues;
  std::vector<Vector3d> eigenvectors;
  if (SolveQEP(M, C, K, &eigenvalues, &eigenvectors)) {
    // Extracts the translations and rotations from the eigenvalues and
    // eigenvectors of the QEP problem.
    for (int i = 0; i < eigenvalues.size(); ++i) {
      Quaterniond quat(eigenvalues[i], axis[0], axis[1], axis[2]);
      quat.normalize();

      soln_rotations->push_back(quat);
      soln_translations->push_back(eigenvectors[i]);
      soln_rotations->push_back(quat);
      soln_translations->push_back(-eigenvectors[i]);
    }
  } else {
    // When there is zero rotation the vector part of the quaternion disappears
    // and it becomes ([0], 1) (where [0] is the zero vector) and the SolveQEP
    // method cannot be used.
    // However from the equations for M, C, and K above we can see that the
    // C and K matrices contain the axis, which is 0 assuming zero rotation,
    // so to solve for the translation we can directly extract the null space
    // of M.
    // Alternatively this can be derived directly from the epipolar
    // constraint, after substituting identity for the rotation matrix.
    eigenvectors.clear();

    JacobiSVD<Matrix3d> svd = M.jacobiSvd(Eigen::ComputeFullV);
    const Vector3d eigenvector(svd.matrixV().col(2));

    soln_rotations->push_back(Quaterniond(AngleAxisd(0.0, axis)));
    soln_translations->push_back(eigenvector);
    soln_rotations->push_back(Quaterniond(AngleAxisd(0.0, axis)));
    soln_translations->push_back(-eigenvector);
  }
}

}  // namespace theia
