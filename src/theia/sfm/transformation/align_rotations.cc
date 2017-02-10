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

#include "theia/sfm/transformation/align_rotations.h"

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Core>
#include <glog/logging.h>
#include <vector>

namespace theia {

namespace {

// A cost function whose error is the difference in rotations after the current
// alignemnt is applied. That is,
//    error = unaligned_rotation * rotation_alignment - gt_rotation.
struct RotationAlignmentError {
  RotationAlignmentError(const Eigen::Vector3d& gt_rotation,
                         const Eigen::Vector3d& unaligned_rotation)
      : gt_rotation_(gt_rotation) {
    // Convert the unaligned rotation to rotation matrix.
    Eigen::Matrix3d unaligned_rotation_mat;
    ceres::AngleAxisToRotationMatrix(
        unaligned_rotation.data(),
        ceres::ColumnMajorAdapter3x3(unaligned_rotation_mat_.data()));
  }

  // Compute the alignment error of the two rotations after applying the
  // rotation transformation.
  template <typename T>
  bool operator()(const T* rotation, T* residuals) const {
    // Convert the rotation transformation to a matrix.
    Eigen::Matrix<T, 3, 3> rotation_mat;
    ceres::AngleAxisToRotationMatrix(
        rotation, ceres::ColumnMajorAdapter3x3(rotation_mat.data()));

    // Apply the rotation transformation.
    const Eigen::Matrix<T, 3, 3> aligned_rotation_mat =
        unaligned_rotation_mat_.cast<T>() * rotation_mat;

    // Convert back to angle axis.
    Eigen::Matrix<T, 3, 1> aligned_rotation;
    ceres::RotationMatrixToAngleAxis(
        ceres::ColumnMajorAdapter3x3(aligned_rotation_mat.data()),
        aligned_rotation.data());

    // Compute the error of the aligned rotation to the gt_rotation.
    residuals[0] = gt_rotation_[0] - aligned_rotation[0];
    residuals[1] = gt_rotation_[1] - aligned_rotation[1];
    residuals[2] = gt_rotation_[2] - aligned_rotation[2];

    return true;
  }

  static ceres::CostFunction* Create(
      const Eigen::Vector3d& gt_rotation,
      const Eigen::Vector3d& unaligned_rotation) {
    return new ceres::AutoDiffCostFunction<RotationAlignmentError, 3, 3>(
        new RotationAlignmentError(gt_rotation, unaligned_rotation));
  }

  Eigen::Vector3d gt_rotation_;
  Eigen::Matrix3d unaligned_rotation_mat_;
};

// Apply the rotation alignment to all rotations in the vector.
void ApplyRotationTransformation(const Eigen::Vector3d& rotation_alignment,
                                 std::vector<Eigen::Vector3d>* rotation) {
  Eigen::Matrix3d rotation_alignment_mat;
  ceres::AngleAxisToRotationMatrix(
      rotation_alignment.data(),
      ceres::ColumnMajorAdapter3x3(rotation_alignment_mat.data()));

  for (int i = 0; i < rotation->size(); i++) {
    // Convert the current rotation to a rotation matrix.
    Eigen::Matrix3d rotation_mat;
    ceres::AngleAxisToRotationMatrix(
        rotation->at(i).data(),
        ceres::ColumnMajorAdapter3x3(rotation_mat.data()));

    // Apply the rotation transformation.
    const Eigen::Matrix3d aligned_rotation =
        rotation_mat * rotation_alignment_mat;

    // Convert back to angle axis.
    ceres::RotationMatrixToAngleAxis(
        ceres::ColumnMajorAdapter3x3(aligned_rotation.data()),
        rotation->at(i).data());
  }
}

}  // namespace

// We solve a nonlinear least squares system to determine a rotation R that will
// align the rotation to the gt_rotation such that rotation * R = gt_rotation.
// This could potentially be set up as a linear system, however, that does not
// guarantee that R will remaind a valid rotation. Instead, we simply use a
// nonlinear system to ensure that R is a valid rotation.
void AlignRotations(const std::vector<Eigen::Vector3d>& gt_rotation,
                    std::vector<Eigen::Vector3d>* rotation) {
  CHECK_EQ(gt_rotation.size(), rotation->size());

  Eigen::Vector3d rotation_alignment = Eigen::Vector3d::Zero();

  // Set up the nonlinear system and adds all residuals.
  ceres::Problem problem;
  for (int i = 0; i < gt_rotation.size(); i++) {
    problem.AddResidualBlock(RotationAlignmentError::Create(gt_rotation[i],
                                                            rotation->at(i)),
                             NULL,
                             rotation_alignment.data());
  }

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.function_tolerance = 0.0;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  VLOG(2) << summary.FullReport();

  // Apply the solved rotation transformation to the rotations.
  ApplyRotationTransformation(rotation_alignment, rotation);
}

}  // namespace theia
