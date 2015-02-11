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

#include "theia/sfm/pose/four_point_relative_pose_partial_rotation.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <glog/logging.h>
#include <math.h>
#include <limits>

#include "theia/alignment/alignment.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::EigenSolver;
using Eigen::JacobiSVD;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::Vector4d;

namespace {

bool SolveQEP(const Matrix4d& M, const Matrix4d& C, const Matrix4d& K,
              std::vector<double>* eigenvalues,
              std::vector<Vector4d>* eigenvectors) {
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

  Matrix4d inv_M;
  static const double kDeterminantThreshold = 1e-7;
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
  Matrix<double, 8, 8> constraint = Matrix<double, 8, 8>::Zero();
  // Set upper-left to - inv(M) * C.
  constraint.block<4, 4>(0, 0) = inv_M * C;
  // Set upper-right to -inv(M) * K.
  constraint.block<4, 4>(0, 4) = inv_M * K;
  // Set lower-left to identity.
  constraint.block<4, 4>(4, 0) = Matrix4d::Identity();

  // Extract the left eigenvectors and values from the constraint matrix.
  static const double kImagEigenValueTolerance = 1e-9;
  EigenSolver<Matrix<double, 8, 8> > eig_solver(constraint);

  for (int i = 0; i < 8; i++) {
    // Only consider the real eigenvalues and corresponding eigenvectors.
    if (fabs(eig_solver.eigenvalues()[i].imag()) < kImagEigenValueTolerance) {
      eigenvalues->push_back(eig_solver.eigenvalues()[i].real());
      eigenvectors->push_back(
          Vector4d(eig_solver.eigenvectors().col(i).tail(4).real()));
    }
  }
  return true;
}

}  // namespace

void FourPointRelativePosePartialRotation(
    const Vector3d& axis,
    const Vector3d image_one_ray_directions[4],
    const Vector3d image_one_ray_origins[4],
    const Vector3d image_two_ray_directions[4],
    const Vector3d image_two_ray_origins[4],
    std::vector<Quaterniond>* soln_rotations,
    std::vector<Vector3d>* soln_translations) {
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
  // the constraint is:
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
  Matrix4d M;
  Matrix4d C;
  Matrix4d K;
  for (int i = 0; i < 4; ++i) {
    const Vector3d& q1(image_one_ray_directions[i]);
    const Vector3d p1(image_one_ray_origins[i].cross(q1));
    const Vector3d& q2(image_two_ray_directions[i]);
    const Vector3d p2(image_two_ray_origins[i].cross(q2));

    M.row(i).head(3) = -q2.cross(q1);
    M.row(i)[3] = q2.dot(p1) + q1.dot(p2);

    C.row(i).head(3) = -2.0 * q2.cross(axis.cross(q1));
    C.row(i)[3] = -2.0 * (q1.dot(axis.cross(p2)) - q2.dot(axis.cross(p1)));

    K.row(i).head(3) = -(2.0 * q1.dot(axis) * q2.cross(axis) - q2.cross(q1));
    K.row(i)[3] =
        -q2.dot(p1) - q1.dot(p2) +
        2.0 * (q2.dot(axis) * p1.dot(axis) + q1.dot(axis) * p2.dot(axis));
  }

  std::vector<double> eigenvalues;
  std::vector<Vector4d> eigenvectors;
  static const double kWTolerance = 1e-7;
  if (SolveQEP(M, C, K, &eigenvalues, &eigenvectors)) {
    // Extracts the translations and rotations from the eigenvalues and
    // eigenvectors of the QEP problem.
    for (int i = 0; i < eigenvalues.size(); ++i) {
      if (fabs(eigenvectors[i].w()) > kWTolerance) {
        Quaterniond quat(eigenvalues[i], axis[0], axis[1], axis[2]);
        quat.normalize();
        soln_rotations->push_back(quat);
        soln_translations->push_back(
            Vector3d(eigenvectors[i].x() / eigenvectors[i].w(),
                     eigenvectors[i].y() / eigenvectors[i].w(),
                     eigenvectors[i].z() / eigenvectors[i].w()));
      }
    }
  } else {
    // When there is zero rotation the vector part of the quaternion disappears
    // and it becomes ([0], 1) (where [0] is the zero vector) and the SolveQEP
    // method cannot be used.
    // However from the equations for M, C, and K above we can see that the
    // C and K matrices contain the axis, which is 0 assuming zero rotation,
    // so to solve for the translation we can directly extract the null space
    // of M.
    // Alternatively this can be derived directly from the generalized epipolar
    // constraint, after substituting identity for the rotation matrix.
    eigenvectors.clear();

    JacobiSVD<Matrix4d> svd = M.jacobiSvd(Eigen::ComputeFullV);
    const Vector4d eigenvector(svd.matrixV().col(3));

    if (fabs(eigenvector[3]) > kWTolerance) {
      soln_rotations->push_back(Quaterniond(AngleAxisd(0.0, axis)));
      soln_translations->push_back(Vector3d(eigenvector[0] / eigenvector[3],
                                            eigenvector[1] / eigenvector[3],
                                            eigenvector[2] / eigenvector[3]));
    }
  }
}

}  // namespace theia
