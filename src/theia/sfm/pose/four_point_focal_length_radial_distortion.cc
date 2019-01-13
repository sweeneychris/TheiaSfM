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
// December 2018

#include <glog/logging.h>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <random>

#include "theia/sfm/pose/four_point_focal_length_radial_distortion.h"
#include "theia/sfm/pose/four_point_focal_length_radial_distortion_helper.h"
#include "theia/util/random.h"

namespace theia {

double sgn(double val) { return (0.0 < val) - (val < 0.0); }

using Matrix34d = Eigen::Matrix<double, 3, 4>;
using Matrix24d = Eigen::Matrix<double, 2, 4>;
using Matrix42d = Eigen::Matrix<double, 4, 2>;
using Vector8d  = Eigen::Matrix<double, 8, 1>;
using Vector5d  = Eigen::Matrix<double, 5, 1>;
using Eigen::Map;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix4d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Matrix;

bool FourPointsPoseFocalLengthRadialDistortion(
    const std::vector<Eigen::Vector2d>& feature_vectors,
    const std::vector<Eigen::Vector3d>& world_points,
    std::vector<Eigen::Matrix3d>* rotations,
    std::vector<Eigen::Vector3d>* translations,
    std::vector<double>* radial_distortions,
    std::vector<double>* focal_lengths) {
  // check that input size of features and world points is 4
  CHECK_GE(feature_vectors.size(), 4);
  CHECK_EQ(feature_vectors.size(), world_points.size());

  Vector4d d;
  Matrix34d world_points_;
  Matrix34d u;  // image points. Will be normalized.
  for (int i = 0; i < 4; ++i) {
    d[i] = feature_vectors[i].squaredNorm();
    world_points_.col(i) = world_points[i];
    u.col(i) = feature_vectors[i].homogeneous();
  }

  const Vector3d t0 = world_points_.rowwise().mean();
  Matrix34d t0_mat;
  t0_mat << t0, t0, t0, t0;

  Matrix4d U;
  U.topRows<3>() = world_points_ - t0_mat;
  U.bottomRows<1>() << 1.0, 1.0, 1.0, 1.0;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      U.topRows<3>(), Eigen::ComputeFullV | Eigen::ComputeFullU);
  Matrix3d R0 = svd.matrixU();

  if (sgn(R0.determinant()) < 0.0) {
    R0.col(0) *= -1;
  }
  R0.transposeInPlace();

  U.topRows<3>() = R0 * U.topRows<3>();
  const double scale =
      U.topRows<3>().array().pow(2).colwise().sum().sqrt().mean();

  U.topRows<3>() /= scale;

  // rescale image points
  const double f0 = u.topRows<2>().array().pow(2).colwise().sum().sqrt().mean();
  u.topRows<2>() /= f0;

  const double k0 = d.array().mean();
  d /= k0;

  Matrix<double, 5, 8> M;
  M.fill(0.0);
  M.row(0).leftCols<4>() = U.col(0);
  M.row(1).rightCols<4>() = U.col(0);
  for (int k = 1; k < 4; ++k) {
    M.row(k + 1).leftCols<4>() = u.row(1).col(k) * U.col(k).transpose();
    M.row(k + 1).rightCols<4>() = -u.row(0).col(k) * U.col(k).transpose();
  }

  Vector5d b;
  b.fill(0.0);
  b.topRows<2>() = u.col(0).topRows<2>();
  Eigen::HouseholderQR<Eigen::MatrixXd> qr;
  qr.compute(M.transpose());

  Matrix<double, 8, 8> Q = qr.householderQ();
  Matrix<double, 8, 5> R = qr.matrixQR().triangularView<Eigen::Upper>();

  Matrix<double, 8, 4> N;
  N.fill(0.0);
  N.leftCols<3>() = Q.rightCols<3>();

  // this random rotation is supposed to make the solver more stable
  static RandomNumberGenerator random_number_gen(42);
  Vector3d rot_vec(random_number_gen.RandDouble(-0.5, 0.5),
                   random_number_gen.RandDouble(-0.5, 0.5),
                   random_number_gen.RandDouble(-0.5, 0.5));
  Eigen::AngleAxisd random_rot(rot_vec.norm(), rot_vec);
  N.leftCols<3>() *= random_rot.toRotationMatrix();
  Matrix<double, 8, 1> x0 =
      Q.leftCols<5>() * (R.topRows<5>().transpose().fullPivLu().solve(b));
  N.rightCols<1>() = x0;

  Matrix<double, 6, 3> C;
  C.fill(0.0);
  Matrix34d UN1 = U.rightCols<3>().transpose() * N.topRows<4>();
  Matrix34d UN2 = U.rightCols<3>().transpose() * N.bottomRows<4>();
  Matrix<double, 6, 9> B;
  B.fill(0.0);

  B.topLeftCorner<3, 3>() = UN1.leftCols<3>();
  B.bottomLeftCorner<3, 3>() = UN2.leftCols<3>();
  B.topRightCorner<3, 1>() = UN1.rightCols<1>();
  B.bottomRightCorner<3, 1>() = UN2.rightCols<1>();

  B.block<3, 4>(0, 3) =
      d.bottomRows<3>().transpose().replicate(4, 1).transpose().cwiseProduct(
          UN1);
  B.block<3, 4>(3, 3) =
      d.bottomRows<3>().transpose().replicate(4, 1).transpose().cwiseProduct(
          UN2);

  B.col(7).topRows<3>() = -u.row(0).rightCols<3>().transpose().cwiseProduct(
      U.row(2).rightCols<3>().transpose());
  B.col(7).bottomRows<3>() = -u.row(1).rightCols<3>().transpose().cwiseProduct(
      U.row(2).rightCols<3>().transpose());

  // fill these guys
  Matrix3d Utmp;
  Utmp.row(0) = U.row(0).rightCols<3>();
  Utmp.row(1) = U.row(1).rightCols<3>();
  Utmp.row(2) = U.row(3).rightCols<3>();

  Matrix3d u1temp;
  u1temp.row(0) = u.row(0).rightCols<3>();
  u1temp.row(1) = u.row(0).rightCols<3>();
  u1temp.row(2) = u.row(0).rightCols<3>();
  Matrix3d u2temp;
  u2temp.row(0) = u.row(1).rightCols<3>();
  u2temp.row(1) = u.row(1).rightCols<3>();
  u2temp.row(2) = u.row(1).rightCols<3>();

  C.block<3, 3>(0, 0) = Utmp.transpose().cwiseProduct(u1temp.transpose());
  C.block<3, 3>(3, 0) = Utmp.transpose().cwiseProduct(u2temp.transpose());

  Matrix<double, 3, 9> D = C.colPivHouseholderQr().solve(B);
  Map<Eigen::RowVectorXd> N_(N.data(), N.size());
  Map<Eigen::RowVectorXd> D_(D.data(), D.size());

  Matrix<double, 64, 1> data;
  data(0) = 0.0;  // used to keep matlab indices, just an index offset of 1
  data.block<32, 1>(1, 0) = N_;
  data.block<27, 1>(33, 0) = D_;
  data(60, 0) = d(0, 0);
  data.bottomRows<3>() = U.topRows<3>().col(0);

  std::vector<Vector5d> valid_solutions;
  FourPointsPoseFocalLengthRadialDistortionSolver(data, &valid_solutions);

  rotations->resize(valid_solutions.size());
  translations->resize(valid_solutions.size());
  radial_distortions->resize(valid_solutions.size());
  focal_lengths->resize(valid_solutions.size());

  for (int i = 0; i < valid_solutions.size(); ++i) {
    const double k = valid_solutions[i][3];
    const double P33 = valid_solutions[i][4];
    Eigen::Vector4d alpha(valid_solutions[i][0], valid_solutions[i][1],
                          valid_solutions[i][2], 1.0);

    Vector8d P12_ = N * alpha;
    Map<Matrix42d> P12(P12_.data(), 4, 2);

    Matrix<double, 9, 1> tmp;
    tmp(0, 0) = alpha[0];
    tmp(1, 0) = alpha[1];
    tmp(2, 0) = alpha[2];
    tmp(3, 0) = k * alpha(0);
    tmp(4, 0) = k * alpha(1);
    tmp(5, 0) = k * alpha(2);
    tmp(6, 0) = k;
    tmp(7, 0) = P33;
    tmp(8, 0) = 1.0;

    const Vector3d P3_124 = D * tmp;

    Vector4d P3;
    P3(0) = P3_124(0);
    P3(1) = P3_124(1);
    P3(2) = P33;
    P3(3) = P3_124(2);

    Matrix34d P;
    P.topRows<2>() = P12.transpose();
    P.bottomRows<1>() = P3;
    P /= P.bottomRows<1>().leftCols<3>().norm();
    const double f = P.topRows<1>().leftCols<3>().norm();
    Matrix3d K = Matrix3d::Identity();
    K(0, 0) = 1. / f;
    K(1, 1) = 1. / f;

    (*focal_lengths)[i] = f * f0;
    (*radial_distortions)[i] = k / k0;

    Matrix34d Rt = K * P;

    if (Rt.topLeftCorner<3, 3>().determinant() < 0.0) Rt *= -1.0;

    // scale radial distortion
    // radial_distortions[i] *= (focal_lengths[i]*focal_lengths[i]);

    Rt.col(3) = Rt.col(3) * scale - Rt.topLeftCorner<3, 3>() * R0 * t0;
    Rt.topLeftCorner<3, 3>() = Rt.topLeftCorner<3, 3>() * R0;

    (*rotations)[i] = Rt.topLeftCorner<3, 3>();
    (*translations)[i] = Rt.col(3);
  }

  return valid_solutions.size() > 0;
}
}
