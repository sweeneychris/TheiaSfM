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

#include "theia/sfm/triangulation/triangulation.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <glog/logging.h>
#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/pose/fundamental_matrix_util.h"
#include "theia/sfm/pose/util.h"

namespace theia {
namespace {

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

// Given either a fundamental or essential matrix and two corresponding images
// points such that ematrix * point2 produces a line in the first image,
// this method finds corrected image points such that
// corrected_point1^t * ematrix * corrected_point2 = 0.
void FindOptimalImagePoints(const Matrix3d& ematrix,
                            const Vector2d& point1,
                            const Vector2d& point2,
                            Vector2d* corrected_point1,
                            Vector2d* corrected_point2) {
  const Vector3d point1_homog = point1.homogeneous();
  const Vector3d point2_homog = point2.homogeneous();

  // A helper matrix to isolate certain coordinates.
  Matrix<double, 2, 3> s_matrix;
  s_matrix <<
      1, 0, 0,
      0, 1, 0;

  const Eigen::Matrix2d e_submatrix = ematrix.topLeftCorner<2, 2>();

  // The epipolar line from one image point in the other image.
  Vector2d epipolar_line1 = s_matrix * ematrix * point2_homog;
  Vector2d epipolar_line2 = s_matrix * ematrix.transpose() * point1_homog;

  const double a = epipolar_line1.transpose() * e_submatrix * epipolar_line2;
  const double b =
      (epipolar_line1.squaredNorm() + epipolar_line2.squaredNorm()) / 2.0;
  const double c = point1_homog.transpose() * ematrix * point2_homog;

  const double d = sqrt(b * b - a * c);

  double lambda = c / (b + d);
  epipolar_line1 -= e_submatrix * lambda * epipolar_line1;
  epipolar_line2 -= e_submatrix.transpose() * lambda * epipolar_line2;

  lambda *=
      (2.0 * d) / (epipolar_line1.squaredNorm() + epipolar_line2.squaredNorm());

  *corrected_point1 = (point1_homog - s_matrix.transpose() * lambda *
                                          epipolar_line1).hnormalized();
  *corrected_point2 = (point2_homog - s_matrix.transpose() * lambda *
                                          epipolar_line2).hnormalized();
}

}  // namespace

// Triangulates 2 posed views
bool Triangulate(const Matrix3x4d& pose1,
                 const Matrix3x4d& pose2,
                 const Vector2d& point1,
                 const Vector2d& point2,
                 Vector4d* triangulated_point) {
  Matrix3d fmatrix;
  FundamentalMatrixFromProjectionMatrices(pose1.data(),
                                          pose2.data(),
                                          fmatrix.data());

  Vector2d corrected_point1, corrected_point2;
  FindOptimalImagePoints(fmatrix, point1, point2,
                         &corrected_point1, &corrected_point2);

  // Now the two points are guaranteed to intersect. We can use the DLT method
  // since it is easy to construct.
  return TriangulateDLT(pose1 ,
                        pose2,
                        corrected_point1,
                        corrected_point2,
                        triangulated_point);
}

// Triangulates a 3D point by determining the closest point between the two
// rays. This method is known to be suboptimal in terms of reprojection error
// but it is extremely fast.
bool TriangulateMidpoint(const Vector3d& ray_origin1,
                         const Vector3d& ray_direction1,
                         const Vector3d& ray_origin2,
                         const Vector3d& ray_direction2,
                         Eigen::Vector4d* triangulated_point) {
  const double dir1_dot_dir2 = ray_direction1.dot(ray_direction2);
  const double dir1_dot_pos = ray_direction1.dot(ray_origin2 - ray_origin1);
  const double dir2_dot_pos = ray_direction2.dot(ray_origin2 - ray_origin1);
  const double scale_part =  1.0 - dir1_dot_dir2 * dir1_dot_dir2;

  const Vector3d scaled_dir1 =
      (dir1_dot_pos - dir1_dot_dir2 * dir2_dot_pos) * ray_direction1;
  const Vector3d scaled_dir2 =
      (dir1_dot_pos * dir1_dot_dir2 - dir2_dot_pos) * ray_direction2;

  // The point is at infinity if the scale division == 0.
  if (scale_part == 0) {
    triangulated_point->head<3>() =
        (ray_origin1 + scaled_dir1 + ray_origin2 + scaled_dir2) / 2.0;
    (*triangulated_point)[3] = 0;
    return true;
  }

  triangulated_point->head<3>() =
      (ray_origin1 + scaled_dir1 / scale_part +
       ray_origin2 + scaled_dir2 / scale_part) / 2.0;
  (*triangulated_point)[3] = 1;
  return true;
}

// Triangulates 2 posed views
bool TriangulateDLT(const Matrix3x4d& pose1,
                    const Matrix3x4d& pose2,
                    const Vector2d& point1,
                    const Vector2d& point2,
                    Vector4d* triangulated_point) {
  Matrix4d design_matrix;
  design_matrix.row(0) = point1[0] * pose1.row(2) - pose1.row(0);
  design_matrix.row(1) = point1[1] * pose1.row(2) - pose1.row(1);
  design_matrix.row(2) = point2[0] * pose2.row(2) - pose2.row(0);
  design_matrix.row(3) = point2[1] * pose2.row(2) - pose2.row(1);

  // Extract nullspace.
  *triangulated_point =
      design_matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();
  return true;
}

// Triangulates N views by computing SVD that minimizes the error.
bool TriangulateNViewSVD(const std::vector<Matrix3x4d>& poses,
                         const std::vector<Vector2d>& points,
                         Vector4d* triangulated_point) {
  CHECK_EQ(poses.size(), points.size());

  MatrixXd design_matrix(3 * points.size(), 4 + points.size());

  for (int i = 0; i < points.size(); i++) {
    design_matrix.block<3, 4>(3 * i, 0) = -poses[i].matrix();
    design_matrix.block<3, 1>(3 * i, 4 + i) = points[i].homogeneous();
  }

  *triangulated_point = design_matrix.jacobiSvd(Eigen::ComputeFullV)
                            .matrixV()
                            .rightCols<1>()
                            .head(4);
  return true;
}

bool TriangulateNView(const std::vector<Matrix3x4d>& poses,
                      const std::vector<Vector2d>& points,
                      Vector4d* triangulated_point) {
  CHECK_EQ(poses.size(), points.size());

  Matrix4d design_matrix = Matrix4d::Zero();
  for (int i = 0; i < points.size(); i++) {
    const Vector3d norm_point = points[i].homogeneous().normalized();
    const Eigen::Matrix<double, 3, 4> cost_term =
        poses[i].matrix() -
        norm_point * norm_point.transpose() * poses[i].matrix();
    design_matrix = design_matrix + cost_term.transpose() * cost_term;
  }

  Eigen::SelfAdjointEigenSolver<Matrix4d> eigen_solver(design_matrix);
  *triangulated_point = eigen_solver.eigenvectors().col(0);
  return true;
}

bool IsTriangulatedPointInFrontOfCameras(
    const FeatureCorrespondence& correspondence,
    const Matrix3d& rotation,
    const Vector3d& position) {
  const Vector3d dir1 = correspondence.feature1.homogeneous();
  const Vector3d dir2 =
      rotation.transpose() * correspondence.feature2.homogeneous();

  const double dir1_sq = dir1.squaredNorm();
  const double dir2_sq = dir2.squaredNorm();
  const double dir1_dir2 = dir1.dot(dir2);
  const double dir1_pos = dir1.dot(position);
  const double dir2_pos = dir2.dot(position);

  return (dir2_sq * dir1_pos - dir1_dir2 * dir2_pos > 0 &&
          dir1_dir2 * dir1_pos - dir1_sq * dir2_pos > 0);
}

}  // namespace theia
