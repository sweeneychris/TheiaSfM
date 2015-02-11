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

#include "theia/sfm/pose/four_point_focal_length.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <glog/logging.h>
#include <vector>

#include "theia/alignment/alignment.h"
#include "theia/sfm/pose/four_point_focal_length_helper.h"

namespace theia {
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace {
void GetRigidTransform(const Matrix<double, 3, 4>& points1,
                       const Matrix<double, 3, 4>& points2,
                       const bool left_handed_coordinates,
                       Eigen::Matrix3d* rotation,
                       Vector3d* translation) {
  // Move the centroid to th origin.
  const Vector3d mean_points1 = points1.rowwise().mean();
  const Vector3d mean_points2 = points2.rowwise().mean();

  const Matrix<double, 3, 4> points1_shifted = points1.colwise() - mean_points1;
  const Matrix<double, 3, 4> points2_shifted = points2.colwise() - mean_points2;

  // Normalize to unit size.
  const Matrix<double, 3, 4> points1_normalized =
      points1_shifted.colwise().normalized();
  const Matrix<double, 3, 4> points2_normalized =
      points2_shifted.colwise().normalized();

  // Compute the necessary rotation from the difference in points.
  Matrix3d rotation_diff = points2_normalized * points1_normalized.transpose();
  Eigen::JacobiSVD<Matrix3d> svd(rotation_diff,
                                 Eigen::ComputeFullU | Eigen::ComputeFullV);

  Matrix3d s = Matrix3d::Zero();
  s(0, 0) = svd.singularValues()(0) < 0 ? -1.0 : 1.0;
  s(1, 1) = svd.singularValues()(1) < 0 ? -1.0 : 1.0;
  const double sign =
      (svd.matrixU() * svd.matrixV().transpose()).determinant() < 0 ? -1.0
                                                                    : 1.0;

  if (left_handed_coordinates) {
    s(2, 2) = -sign;
  } else {
    s(2, 2) = sign;
  }
  *rotation = svd.matrixU() * s * svd.matrixV().transpose();
  *translation = - *rotation * mean_points1 + mean_points2;
}

}  // namespace

int FourPointPoseAndFocalLength(
    const std::vector<Eigen::Vector2d>& feature_vectors,
    const std::vector<Eigen::Vector3d>& world_points_vector,
    std::vector<Eigen::Matrix<double, 3, 4> >* projection_matrices) {
  Eigen::Map<const Matrix<double, 2, 4> > features(feature_vectors[0].data());
  Eigen::Map<const Matrix<double, 3, 4> > world_points(world_points_vector[0]
                                                           .data());

  // Normalize the points such that the mean = 0, variance = sqrt(2.0).
  const Vector3d mean_world_point = world_points.rowwise().mean();
  Eigen::Matrix<double, 3, 4> world_point_normalized =
      world_points.colwise() - mean_world_point;
  const double world_point_variance =
      world_point_normalized.colwise().norm().mean();
  world_point_normalized /= world_point_variance;

  // Scale 2D data so variance = sqrt(2.0).
  const double features_variance = features.colwise().norm().mean();
  Eigen::Matrix<double, 2, 4> features_normalized =
      features / features_variance;

  // Precompute monomials.
  const double glab = (world_point_normalized.col(0) -
                       world_point_normalized.col(1)).squaredNorm();
  const double glac = (world_point_normalized.col(0) -
                       world_point_normalized.col(2)).squaredNorm();
  const double glad = (world_point_normalized.col(0) -
                       world_point_normalized.col(3)).squaredNorm();
  const double glbc = (world_point_normalized.col(1) -
                       world_point_normalized.col(2)).squaredNorm();
  const double glbd = (world_point_normalized.col(1) -
                       world_point_normalized.col(3)).squaredNorm();
  const double glcd = (world_point_normalized.col(2) -
                       world_point_normalized.col(3)).squaredNorm();

  if (glab * glac * glad * glbc * glbd * glcd < 1e-15) {
    return -1;
  }

  // Call the helper function.
  std::vector<double> focal_length;
  std::vector<Vector3d> depths;

  FourPointFocalLengthHelper(glab, glac, glad, glbc, glbd, glcd,
                             features_normalized, &focal_length, &depths);

  if (focal_length.size() == 0) {
    return -1;
  }

  // Get the rotation and translation.
  for (int i = 0; i < focal_length.size(); i++) {
    // Create world points in camera coordinate system.
    Matrix<double, 3, 4> adjusted_world_points;
    adjusted_world_points.block<2, 4>(0, 0) = features_normalized;
    adjusted_world_points.row(2).setConstant(focal_length[i]);
    adjusted_world_points.col(1) *= depths[i].x();
    adjusted_world_points.col(2) *= depths[i].y();
    adjusted_world_points.col(3) *= depths[i].z();

    // Fix the scale.
    Matrix<double, 6, 1> d;
    d(0) = sqrt(glab / (adjusted_world_points.col(0) -
                        adjusted_world_points.col(1)).squaredNorm());
    d(1) = sqrt(glac / (adjusted_world_points.col(0) -
                        adjusted_world_points.col(2)).squaredNorm());
    d(2) = sqrt(glad / (adjusted_world_points.col(0) -
                        adjusted_world_points.col(3)).squaredNorm());
    d(3) = sqrt(glbc / (adjusted_world_points.col(1) -
                        adjusted_world_points.col(2)).squaredNorm());
    d(4) = sqrt(glbd / (adjusted_world_points.col(1) -
                        adjusted_world_points.col(3)).squaredNorm());
    d(5) = sqrt(glcd / (adjusted_world_points.col(2) -
                        adjusted_world_points.col(3)).squaredNorm());

    const double gta = d.mean();

    adjusted_world_points *= gta;

    // Get the transformation by aligning the points.
    Matrix3d rotation;
    Vector3d translation;
    GetRigidTransform(world_point_normalized, adjusted_world_points, false,
                      &rotation, &translation);

    translation =
        world_point_variance * translation - rotation * mean_world_point;

    focal_length[i] *= features_variance;

    Matrix<double, 3, 4> transformation_matrix;
    transformation_matrix.block<3, 3>(0, 0) = rotation;
    transformation_matrix.col(3) = translation;
    Matrix3d camera_matrix =
        Eigen::DiagonalMatrix<double, 3>(focal_length[i], focal_length[i], 1.0);
    projection_matrices->push_back(camera_matrix * transformation_matrix);
  }
  return projection_matrices->size();
}

}  // namespace theia
