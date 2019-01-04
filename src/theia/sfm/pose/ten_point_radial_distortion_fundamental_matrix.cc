// Copyright (C) 2019 The Regents of the University of California (Regents).
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

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// January 2019

#include "theia/sfm/pose/ten_point_radial_distortion_fundamental_matrix.h"
#include "theia/sfm/pose/ten_point_radial_distortion_fundamental_matrix_helper.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Jacobi>

namespace theia {

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix3d;
using Eigen::Vector4d;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Vector10d = Eigen::Matrix<double, 10, 1>;
using Matrix102d = Eigen::Matrix<double, 10, 2>;
using Matrix210d = Eigen::Matrix<double, 2, 10>;

bool TenPointRadialDistortionFundamentalMatrix(
    const std::vector<Eigen::Vector2d>& normalized_feature_points_left,
    const std::vector<Eigen::Vector2d>& normalized_feature_points_right,
    std::vector<RadialFundamentalMatrixResult>* results,
    double lmin, double lmax) {

  Matrix102d X;
  Matrix102d U;
  for (int i = 0; i < 10; ++i) {
    X.row(i) = normalized_feature_points_left[i];
    U.row(i) = normalized_feature_points_right[i];
  }

  Vector10d Z1;
  Vector10d Z2;
  Matrix<double, 10, 16> A;
  Z1.array() =
      X.col(0).array() * X.col(0).array() + X.col(1).array() * X.col(1).array();
  Z2.array() =
      U.col(0).array() * U.col(0).array() + U.col(1).array() * U.col(1).array();

  A.col(0).array() = X.col(0).array() * U.col(0).array();
  A.col(1).array() = X.col(0).array() * U.col(1).array();
  A.col(2).array() = X.col(1).array() * U.col(0).array();
  A.col(3).array() = X.col(1).array() * U.col(1).array();
  A.col(4).array() = U.col(0).array() * Z1.array();
  A.col(5).array() = U.col(0).array();
  A.col(6).array() = U.col(1).array() * Z1.array();
  A.col(7).array() = U.col(1).array();
  A.col(8).array() = X.col(0).array() * Z2.array();
  A.col(9).array() = X.col(0).array();
  A.col(10).array() = X.col(1).array() * Z2.array();
  A.col(11).array() = X.col(1).array();
  A.col(12).array() = Z1.array() * Z2.array();
  A.col(13).array() = Z1.array();
  A.col(14).array() = Z2.array();
  A.col(15).fill(1.0);

  const Matrix<double, 10, 6> Mr =
      A.block<10, 10>(0, 0).lu().solve(A.block<10, 6>(0, 10));

  Matrix<double, 29, 1> params;

  params << Mr(5, 0), Mr(5, 1), -Mr(4, 0), -Mr(4, 1), Mr(5, 2), Mr(5, 3),
      Mr(5, 4) - Mr(4, 2), Mr(5, 5) - Mr(4, 3), -Mr(4, 4), -Mr(4, 5), Mr(7, 0),
      Mr(7, 1), -Mr(6, 0), -Mr(6, 1), Mr(7, 2), Mr(7, 3), Mr(7, 4) - Mr(6, 2),
      Mr(7, 5) - Mr(6, 3), -Mr(6, 4), -Mr(6, 5), Mr(9, 0), Mr(9, 1) - Mr(8, 0),
      -Mr(8, 1), Mr(9, 2), Mr(9, 4), Mr(9, 3) - Mr(8, 2), -Mr(8, 3),
      Mr(9, 5) - Mr(8, 4), -Mr(8, 5);

  Matrix210d Ls;
  int nsols = TenPointRadialDistortionFundamentalMatrixHelper(params, Ls);

  if (nsols > 0) {
    Vector4d m1;
    Vector6d m2;
    Vector6d m3;
    Vector10d b;

    b << Mr(5, 0), Mr(5, 1), -Mr(4, 0), -Mr(4, 1), Mr(5, 2), Mr(5, 3),
        Mr(5, 4) - Mr(4, 2), Mr(5, 5) - Mr(4, 3), -Mr(4, 4), -Mr(4, 5);

    for (int i = 0; i < nsols; i++) {
      const double l1 = Ls(0, i);
      const double l2 = Ls(1, i);

      if (l1 < lmin || l2 < lmin || l1 > lmax || l2 > lmax) {
        continue;
      }

      const double l1l1 = l1 * l1;
      const double l1l2 = l1 * l2;
      double f23;

      m1 << l1l2, l1, l2, 1;
      m2 << l1l2 * l1, l1l1, l1l2, l1, l2, 1;
      f23 = -b.block<6, 1>(4, 0).dot(m2) / b.block<4, 1>(0, 0).dot(m1);
      m3 << l2 * f23, f23, l1l2, l1, l2, 1;

      RadialFundamentalMatrixResult res;
      res.l1 = l1;
      res.l2 = l2;
      res.F <<
          m3.dot(-Mr.row(0)), m3.dot(-Mr.row(1)), m3.dot(-Mr.row(9)),
          m3.dot(-Mr.row(2)), m3.dot(-Mr.row(3)), f23,
          m3.dot(-Mr.row(5)), m3.dot(-Mr.row(7)), 1.0;
      results->push_back(res);
    }

    return true;
  } else
    return false;
}
}  // theia namespace
