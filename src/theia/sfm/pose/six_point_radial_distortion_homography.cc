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

#include "theia/sfm/pose/six_point_radial_distortion_homography.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Jacobi>

namespace theia {

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Array51d = Eigen::Array<double, 5, 1>;
using Matrix68d = Eigen::Matrix<double, 6, 8>;
using Matrix62d = Eigen::Matrix<double, 6, 2>;
using Matrix65d = Eigen::Matrix<double, 6, 5>;
using Eigen::Matrix3d;

bool IsNearZero(double val) {
  return (+val < (100.0 * std::numeric_limits<double>::epsilon())) &&
         (-val < (100.0 * std::numeric_limits<double>::epsilon()));
}

bool SixPointRadialDistortionHomography(
    const std::vector<Eigen::Vector2d>& normalized_feature_points_left,
    const std::vector<Eigen::Vector2d>& normalized_feature_points_right,
    std::vector<RadialHomographyResult>* results, const double lmin,
    const double lmax) {
  Matrix62d X;
  Matrix62d U;

  for (int i = 0; i < 6; ++i) {
    X.row(i) = normalized_feature_points_left[i];
    U.row(i) = normalized_feature_points_right[i];
  }

  Matrix68d M;
  Vector6d u2 = U.col(0).array().square() + U.col(1).array().square();

  M.col(0) = -X.col(1).array() * U.col(0).array();
  M.col(1) = -X.col(1).array() * U.col(1).array();
  M.col(2) = -X.col(1);
  M.col(3) = X.col(0).array() * U.col(0).array();
  M.col(4) = X.col(0).array() * U.col(1).array();
  M.col(5) = X.col(0);
  M.col(6) = -X.col(1).array() * u2.array();
  M.col(7) = X.col(0).array() * u2.array();

  Eigen::JacobiSVD<Matrix68d, Eigen::FullPivHouseholderQRPreconditioner> Svd1(
      M, Eigen::ComputeFullV);
  const Eigen::Matrix<double, 8, 8>& V1 = Svd1.matrixV();

  const double a = -V1(2, 6) * V1(7, 6) + V1(5, 6) * V1(6, 6);
  const double b = -V1(2, 6) * V1(7, 7) - V1(2, 7) * V1(7, 6) +
                   V1(5, 6) * V1(6, 7) + V1(5, 7) * V1(6, 6);
  const double c = -V1(2, 7) * V1(7, 7) + V1(5, 7) * V1(6, 7);
  const double d = b * b - 4.0 * a * c;

  int nsols = 0;
  Vector2d rs;

  if (IsNearZero(d)) {
    nsols = 1;
    rs(0) = (-b) / (2.0 * a);
  } else if (d > 0.0) {
    nsols = 2;
    double d2 = std::sqrt(d);
    rs(0) = (-b + d2) / (2.0 * a);
    rs(1) = (-b - d2) / (2.0 * a);
  } else {
    return false;
  }

  const Vector6d x2 = X.col(0).array().square() + X.col(1).array().square();
  Vector6d u3, r;
  Matrix<double, 8, 1> n;
  Matrix65d T;
  T.col(0) = -M.col(3);
  T.col(1) = -M.col(4);

  for (int i = 0; i < nsols; i++) {
    n = rs(i) * V1.col(6) + V1.col(7);
    const double l2 = n(6) / n(2);
    // skip this solution early if radial distortion is spurious
    if (l2 < lmin || l2 > lmax) {
      continue;
    }

    u3 = u3.Ones() + l2 * u2;
    r = n(0) * U.col(0) + n(1) * U.col(1) + n(2) * u3;

    T.col(2) = -X.col(0).array() * u3.array();
    T.col(3) = x2.array() * r.array();
    T.col(4) = r;

    Eigen::JacobiSVD<Matrix65d> Svd2(T, Eigen::ComputeFullV);
    Matrix<double, 5, 1> v2 = Svd2.matrixV().col(4);

    v2.head(4) /= v2(4);
    const double l1 = v2(3);
    // skip this solution early if radial distortion is spurious
    if (l1 < lmin || l1 > lmax) {
      continue;
    }

    RadialHomographyResult res;
    // fill homograhapy
    res.H << n(0), n(1), n(2), n(3), n(4), n(5), v2(0), v2(1), v2(2);
    // fill radial distortion values
    res.l1 = l1;
    res.l2 = l2;
    results->push_back(res);
  }

  return nsols > 0;
}

// Some helper functions
void DistortPoint(const Vector3d& point_in_camera, const double focal_length,
                  const double radial_distortion, Vector2d& distorted_point) {
  Vector2d point_in_cam_persp_div;
  point_in_cam_persp_div[0] =
      focal_length * point_in_camera[0] / point_in_camera[2];
  point_in_cam_persp_div[1] =
      focal_length * point_in_camera[1] / point_in_camera[2];

  // see also division_undistortion_camera_model.h
  const double r_u_sq = point_in_cam_persp_div.dot(point_in_cam_persp_div);

  const double denom = 2.0 * radial_distortion * r_u_sq;
  const double inner_sqrt = 1.0 - 4.0 * radial_distortion * r_u_sq;

  // If the denominator is nearly zero then we can evaluate the distorted
  // coordinates as k or r_u^2 goes to zero. Both evaluate to the identity.
  if (std::abs(denom) < std::numeric_limits<double>::epsilon() ||
      inner_sqrt < 0.0) {
    distorted_point = point_in_cam_persp_div;
  } else {
    const double scale = (1.0 - std::sqrt(inner_sqrt)) / denom;
    distorted_point = point_in_cam_persp_div * scale;
  }
}

void UndistortPoint(const Vector2d& distorted_point, const double focal_length,
                    const double radial_distortion,
                    Vector3d& undistorted_point) {
  double px = distorted_point[0];
  double py = distorted_point[1];
  // The squared radius of the distorted image point.
  const double r_d_sq = distorted_point.dot(distorted_point);

  const double undistortion = 1.0 / (1.0 + radial_distortion * r_d_sq);
  undistorted_point[0] = px * undistortion / focal_length;
  undistorted_point[1] = py * undistortion / focal_length;
  undistorted_point[2] = 1.0;
}

void ProjectCameraToCamera(const Matrix3d& H, const Vector3d& X, Vector3d* Y) {
  (*Y) = H * X;
  (*Y) /= (*Y)(2);
}

double CheckRadialSymmetricError(
    const RadialHomographyResult& radial_homography, const Vector2d& pt_left,
    const Vector2d& pt_right, const double focal_length1,
    const double focal_length2) {
  // todo: in the estimator this gets calculated for each sample every time
  const double l1_scaled =
      radial_homography.l1 / (focal_length1 * focal_length1);
  const double l2_scaled =
      radial_homography.l2 / (focal_length2 * focal_length2);

  Vector3d bearing_vector_left, bearing_vector_right;
  UndistortPoint(pt_left, focal_length1, l1_scaled, bearing_vector_left);
  UndistortPoint(pt_right, focal_length2, l2_scaled, bearing_vector_right);

  Vector3d ray2_in_1, ray1_in_2;
  ProjectCameraToCamera(radial_homography.H, bearing_vector_right, &ray2_in_1);
  ProjectCameraToCamera(radial_homography.H.inverse(), bearing_vector_left,
                        &ray1_in_2);

  Vector2d pt_left_from_right, pt_right_from_left;
  DistortPoint(ray2_in_1, focal_length1, l1_scaled, pt_left_from_right);
  DistortPoint(ray1_in_2, focal_length2, l2_scaled, pt_right_from_left);

  const Vector2d diff_left = pt_left - pt_left_from_right;
  const Vector2d diff_right = pt_right - pt_right_from_left;

  const double square_left = diff_left.dot(diff_left);
  const double square_right = diff_right.dot(diff_right);

  double sym_error = 0.5 * (square_left + square_right);

  return sym_error;
}
}
