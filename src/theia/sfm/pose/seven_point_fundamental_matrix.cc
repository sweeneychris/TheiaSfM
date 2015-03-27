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
#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::Matrix;

namespace {

// Sets up the constraint y^t * F * x = 0 such that M * F_v = 0 where M is a 7x9
// matrix and F_v is the vector containing the entries of F.
Matrix<double, 7, 9> SetupEpipolarConstraint(
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

}  // namespace

bool SevenPointFundamentalMatrix(
    const std::vector<Eigen::Vector2d>& image1_points,
    const std::vector<Eigen::Vector2d>& image2_points,
    std::vector<Eigen::Matrix3d>* fundamental_matrices) {
  CHECK_EQ(image1_points.size(), 7);
  CHECK_EQ(image2_points.size(), 7);
  CHECK_NOTNULL(fundamental_matrices)->clear();

  std::vector<Eigen::Vector2d> norm_img1_points(image1_points.size());
  std::vector<Eigen::Vector2d> norm_img2_points(image2_points.size());

  // Normalize the image points.
  Eigen::Matrix3d img1_norm_mat, img2_norm_mat;
  NormalizeImagePoints(image1_points, &norm_img1_points, &img1_norm_mat);
  NormalizeImagePoints(image2_points, &norm_img2_points, &img2_norm_mat);

  const Matrix<double, 7, 9>& epipolar_constraint =
      SetupEpipolarConstraint(norm_img1_points, norm_img2_points);

  const Eigen::FullPivLU<Matrix<double, 7, 9> > lu(epipolar_constraint);
  if (lu.dimensionOfKernel() != 2) {
    return false;
  }

  // Represent F in terms of its null space such that F = x * F1' + (1 - x) * F2
  // where F1 and F2 are vectors in the null space of F. Note that this can also
  // be parameterized such that:
  //   F = x * F1' + (1 - x) * F2 = x * (F1' - F2) + F2 = x * F1 + F2.
  const Matrix<double, 9, 2>& null_space = lu.kernel();
  const Matrix<double, 9, 1> F1_vec = null_space.col(0) - null_space.col(1);
  const Eigen::Map<const Eigen::Matrix3d> F1(F1_vec.data());
  const Eigen::Map<const Eigen::Matrix3d> F2(null_space.col(1).data());

  // This is the cubic equation resulting from det(x * F1 + F2) = 0.
  const Eigen::Vector4d determinant_constraint(
      -(F2(1, 2) * F2(2, 1) - F2(1, 1) * F2(2, 2)) * F2(0, 0) +
          (F2(0, 2) * F2(2, 1) - F2(0, 1) * F2(2, 2)) * F2(1, 0) -
          (F2(0, 2) * F2(1, 1) - F2(0, 1) * F2(1, 2)) * F2(2, 0),
      -(F2(1, 2) * F2(2, 1) - F2(1, 1) * F2(2, 2)) * F1(0, 0) +
          (F2(0, 2) * F2(2, 1) - F2(0, 1) * F2(2, 2)) * F1(1, 0) -
          (F2(0, 2) * F2(1, 1) - F2(0, 1) * F2(1, 2)) * F1(2, 0) +
          (F1(2, 2) * F2(1, 1) - F1(2, 1) * F2(1, 2) - F1(1, 2) * F2(2, 1) +
           F1(1, 1) * F2(2, 2)) *
              F2(0, 0) -
          (F1(2, 2) * F2(0, 1) - F1(2, 1) * F2(0, 2) - F1(0, 2) * F2(2, 1) +
           F1(0, 1) * F2(2, 2)) *
              F2(1, 0) +
          (F1(1, 2) * F2(0, 1) - F1(1, 1) * F2(0, 2) - F1(0, 2) * F2(1, 1) +
           F1(0, 1) * F2(1, 2)) *
              F2(2, 0),
      (F1(2, 2) * F2(1, 1) - F1(2, 1) * F2(1, 2) - F1(1, 2) * F2(2, 1) +
       F1(1, 1) * F2(2, 2)) *
              F1(0, 0) -
          (F1(2, 2) * F2(0, 1) - F1(2, 1) * F2(0, 2) - F1(0, 2) * F2(2, 1) +
           F1(0, 1) * F2(2, 2)) *
              F1(1, 0) +
          (F1(1, 2) * F2(0, 1) - F1(1, 1) * F2(0, 2) - F1(0, 2) * F2(1, 1) +
           F1(0, 1) * F2(1, 2)) *
              F1(2, 0) -
          (F1(1, 2) * F1(2, 1) - F1(1, 1) * F1(2, 2)) * F2(0, 0) +
          (F1(0, 2) * F1(2, 1) - F1(0, 1) * F1(2, 2)) * F2(1, 0) -
          (F1(0, 2) * F1(1, 1) - F1(0, 1) * F1(1, 2)) * F2(2, 0),
      -(F1(1, 2) * F1(2, 1) - F1(1, 1) * F1(2, 2)) * F1(0, 0) +
          (F1(0, 2) * F1(2, 1) - F1(0, 1) * F1(2, 2)) * F1(1, 0) -
          (F1(0, 2) * F1(1, 1) - F1(0, 1) * F1(1, 2)) * F1(2, 0));

  // Solve the cubic equation for x.
  double roots[3];
  const int num_solutions = SolveCubicReals(determinant_constraint(3),
                                            determinant_constraint(2),
                                            determinant_constraint(1),
                                            determinant_constraint(0),
                                            roots);

  for (int i = 0; i < num_solutions; i++) {
    // Compose the fundamental matrix solution from the null space and
    // determinant constraint: F = x * F1 + F2;
    fundamental_matrices->emplace_back(img2_norm_mat.transpose() *
                                       (roots[i] * F1 + F2) * img1_norm_mat);
  }
  return fundamental_matrices->size() > 0;
}

}  // namespace theia
