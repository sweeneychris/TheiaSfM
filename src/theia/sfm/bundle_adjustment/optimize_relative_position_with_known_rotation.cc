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

#include "theia/sfm/bundle_adjustment/optimize_relative_position_with_known_rotation.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <algorithm>
#include <vector>

#include "theia/math/util.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/triangulation/triangulation.h"

namespace theia {
namespace {

// Creates the constraint matrix such that ||A * t|| is minimized, where A is
// R_i * f_i x R_j * f_j. Given known rotations, we can solve for the
// relative translation from this constraint matrix.
void CreateConstraintMatrix(
    const std::vector<FeatureCorrespondence>& correspondences,
    const Eigen::Vector3d& rotation1,
    const Eigen::Vector3d& rotation2,
    Eigen::MatrixXd* constraint_matrix) {
  constraint_matrix->resize(3, correspondences.size());

  Eigen::Matrix3d rotation_matrix1;
  ceres::AngleAxisToRotationMatrix(
      rotation1.data(), ceres::ColumnMajorAdapter3x3(rotation_matrix1.data()));
  Eigen::Matrix3d rotation_matrix2;
  ceres::AngleAxisToRotationMatrix(
      rotation2.data(), ceres::ColumnMajorAdapter3x3(rotation_matrix2.data()));

  for (int i = 0; i < correspondences.size(); i++) {
    const Eigen::Vector3d rotated_feature1 =
        rotation_matrix1.transpose() *
        correspondences[i].feature1.homogeneous();
    const Eigen::Vector3d rotated_feature2 =
        rotation_matrix2.transpose() *
        correspondences[i].feature2.homogeneous();

    constraint_matrix->col(i) =
        rotated_feature2.cross(rotated_feature1).transpose() *
        rotation_matrix1.transpose();
  }
}

// Determines if the majority of the points are in front of the cameras. This is
// useful for determining the sign of the relative position. Returns true if
// more than 50% of correspondences are in front of both cameras and false
// otherwise.
bool MajorityOfPointsInFrontOfCameras(
    const std::vector<FeatureCorrespondence>& correspondences,
    const Eigen::Vector3d& rotation1,
    const Eigen::Vector3d& rotation2,
    const Eigen::Vector3d& relative_position) {
  // Compose the relative rotation.
  Eigen::Matrix3d rotation_matrix1, rotation_matrix2;
  ceres::AngleAxisToRotationMatrix(
      rotation1.data(), ceres::ColumnMajorAdapter3x3(rotation_matrix1.data()));
  ceres::AngleAxisToRotationMatrix(
      rotation2.data(), ceres::ColumnMajorAdapter3x3(rotation_matrix2.data()));
  const Eigen::Matrix3d relative_rotation_matrix =
      rotation_matrix2 * rotation_matrix1.transpose();

  // Tests all points for cheirality.
  int num_points_in_front_of_cameras = 0;
  for (const FeatureCorrespondence& match : correspondences) {
    if (IsTriangulatedPointInFrontOfCameras(match,
                                            relative_rotation_matrix,
                                            relative_position)) {
      ++num_points_in_front_of_cameras;
    }
  }

  return num_points_in_front_of_cameras > (correspondences.size() / 2);
}

}  // namespace

// Given known camera rotations and feature correspondences, this method solves
// for the relative translation that optimizes the epipolar error
// f_i * E * f_j^t = 0.
bool OptimizeRelativePositionWithKnownRotation(
    const std::vector<FeatureCorrespondence>& correspondences,
    const Eigen::Vector3d& rotation1,
    const Eigen::Vector3d& rotation2,
    Eigen::Vector3d* relative_position) {
  CHECK_NOTNULL(relative_position);

  *relative_position = Eigen::Vector3d::Random().normalized();

  // Constants used for the IRLS solving.
  const double eps = 1e-5;
  const int kMaxIterations = 100;
  const int kMaxInnerIterations = 10;
  const double kMinWeight = 1e-7;

  // Create the constraint matrix from the known correspondences and rotations.
  Eigen::MatrixXd constraint_matrix;
  CreateConstraintMatrix(correspondences,
                         rotation1,
                         rotation2,
                         &constraint_matrix);

  // Initialize the weighting terms for each correspondence.
  Eigen::VectorXd weights(correspondences.size());
  weights.setConstant(1.0);

  // Solve for the relative positions using a robust IRLS.
  double cost = 0;
  int num_inner_iterations = 0;
  for (int i = 0;
       i < kMaxIterations && num_inner_iterations < kMaxInnerIterations;
       i++) {
    // Limit the minimum weight at kMinWeight.
    weights = (weights.array() < kMinWeight).select(kMinWeight, weights);

    // Apply the weights to the constraint matrix.
    const Eigen::Matrix3d lhs = constraint_matrix *
                                weights.asDiagonal().inverse() *
                                constraint_matrix.transpose();

    // Solve for the relative position which is the null vector of the weighted
    // constraints.
    const Eigen::Vector3d new_relative_position =
        lhs.jacobiSvd(Eigen::ComputeFullU).matrixU().rightCols<1>();

    // Update the weights based on the current errors.
    weights =
        (new_relative_position.transpose() * constraint_matrix).array().abs();

    // Compute the new cost.
    const double new_cost = weights.sum();

    // Check for convergence.
    const double delta = std::max(std::abs(cost - new_cost),
                                  1 - new_relative_position.squaredNorm());

    // If we have good convergence, attempt an inner iteration.
    if (delta <= eps) {
      ++num_inner_iterations;
    } else {
      num_inner_iterations = 0;
    }

    cost = new_cost;
    *relative_position = new_relative_position;
  }

  // The position solver above does not consider the sign of the relative
  // position. We can determine the sign by choosing the sign that puts the most
  // points in front of the camera.
  if (!MajorityOfPointsInFrontOfCameras(correspondences,
                                        rotation1,
                                        rotation2,
                                        *relative_position)) {
    *relative_position *= -1.0;
  }

  return true;
}

}  // namespace theia
