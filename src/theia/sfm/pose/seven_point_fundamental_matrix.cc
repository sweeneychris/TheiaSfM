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

#include "theia/sfm/pose/seven_point_fundamental_matrix.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <glog/logging.h>
#include <vector>

#include "theia/math/closed_form_polynomial_solver.h"

namespace theia {

using Eigen::Matrix;

namespace {

// Sets up the constraint y^t * F * x = 0 such that M * F_v = 0 where M is a 7x9
// matrix and F_v is the vector containing the entries of F.
const Matrix<double, 7, 9> SetupEpipolarConstraint(
    const std::vector<Eigen::Vector2d>& image1_points,
    const std::vector<Eigen::Vector2d>& image2_points) {
  Matrix<double, 7, 9> epipolar_constraint;
  for (int i = 0; i < 7; i++) {
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

  return epipolar_constraint;
}

// Multiplys two degree 1 polynomials:
//   (alpha * x1 + x2) * (alpha * y1 + y2)
//     = alpha^2 * x1 * y1 + alpha * (x1 * y2 + x2 * y1) + x2 * y2.
Eigen::Vector3d MultiplyDegOnePoly(const double x1,
                                   const double x2,
                                   const double y1,
                                   const double y2) {
  return Eigen::Vector3d(x1 * y1, x1 * y2 + x2 * y1, y2 * y2);
}

// Multiply a degree 1 and degree 2 polynomial:
//   (alpha * x1 + x2) * (alpha^2 * y1 + alpha * y2 + y3)
//     =   alpha^3 * x1 * y1
//       + alpha^2 (x1 * y2 + x2 * y1)
//       + alpha * (x1 * y3 + x2 * y2)
//       + x2 * y3
Eigen::Vector4d MultiplyDegTwoDegOnePoly(const double x1,
                                         const double x2,
                                         const Eigen::Vector3d& deg1_poly) {
  return Eigen::Vector4d(x1 * deg1_poly(0),
                         x1 * deg1_poly(1) + x2 * deg1_poly(0),
                         x1 * deg1_poly(2) + x2 * deg1_poly(1),
                         x2 * deg1_poly(2));
}

}  // namespace

bool SevenPointFundamentalMatrix(
    const std::vector<Eigen::Vector2d>& image1_points,
    const std::vector<Eigen::Vector2d>& image2_points,
    std::vector<Eigen::Matrix3d>* fundamental_matrices) {
  CHECK_EQ(image1_points.size(), 7);
  CHECK_EQ(image2_points.size(), 7);
  CHECK_NOTNULL(fundamental_matrices)->clear();

  const Matrix<double, 7, 9>& epipolar_constraint =
      SetupEpipolarConstraint(image1_points, image2_points);

  const Eigen::FullPivLU<Matrix<double, 7, 9> > lu(epipolar_constraint);
  if (lu.dimensionOfKernel() != 2) {
    return false;
  }

  // Represent F in terms of its null space such that F = x * F1' + (1 - x) * F2
  // where F1 and F2 are vectors in the null space of F. Note that this can also
  // be parameterized such that:
  //   F = x * F1' + (1 - x) * F2 = x * (F1' - F2) + F2 = x * F1 + F2.
  const Matrix<double, 9, 2>& null_space = lu.kernel();
  const Eigen::Matrix<double, 9, 1> F1 =
      null_space.col(0) - null_space.col(1);
  const Eigen::Matrix<double, 9, 1>& F2 = null_space.col(1);

  // Enforce that det(F) = 0. This results in a cubic equation in the unknown x.
  const Eigen::Vector4d determinant_constraint =
      MultiplyDegTwoDegOnePoly(F1(0), F2(0),
                               MultiplyDegOnePoly(F1(4), F2(4), F1(8), F2(8)) -
                               MultiplyDegOnePoly(F1(5), F2(5), F1(7), F2(7))) -
      MultiplyDegTwoDegOnePoly(F1(1), F2(1),
                               MultiplyDegOnePoly(F1(3), F2(3), F1(8), F2(8)) -
                               MultiplyDegOnePoly(F1(5), F2(5), F1(6), F2(6))) +
      MultiplyDegTwoDegOnePoly(F1(2), F2(2),
                               MultiplyDegOnePoly(F1(3), F2(3), F1(7), F2(7)) -
                               MultiplyDegOnePoly(F1(4), F2(4), F1(6), F2(6)));

  // Solve the cubic equation for x.
  double roots[3];
  const int num_solutions = SolveCubicReals(determinant_constraint(0),
                                            determinant_constraint(1),
                                            determinant_constraint(2),
                                            determinant_constraint(3),
                                            roots);
  const Eigen::Map<const Eigen::Matrix3d> F1_map(F1.data());
  const Eigen::Map<const Eigen::Matrix3d> F2_map(F2.data());
  for (int i = 0; i < num_solutions; i++) {
    // Compose the fundamental matrix solution from the null space and
    // determinant constraint: F = x * F1 + F2;
    fundamental_matrices->emplace_back(roots[i] * F1_map + F2_map);
  }
  return fundamental_matrices->size() > 0;
}

}  // namespace theia
