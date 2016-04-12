// Copyright (C) 2014 The Regents of the University of California (Regents).
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

#include "theia/sfm/camera/projection_matrix_utils.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <glog/logging.h>

#include "theia/math/matrix/rq_decomposition.h"
#include "theia/sfm/pose/util.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

void IntrinsicsToCalibrationMatrix(const double focal_length,
                                   const double skew,
                                   const double aspect_ratio,
                                   const double principal_point_x,
                                   const double principal_point_y,
                                   Matrix3d* calibration_matrix) {
  *calibration_matrix <<
      focal_length, skew, principal_point_x,
      0, focal_length * aspect_ratio, principal_point_y,
      0, 0, 1.0;
}

void CalibrationMatrixToIntrinsics(const Matrix3d& calibration_matrix,
                                   double* focal_length,
                                   double* skew,
                                   double* aspect_ratio,
                                   double* principal_point_x,
                                   double* principal_point_y) {
  CHECK_NE(calibration_matrix(2, 2), 0);
  *focal_length = calibration_matrix(0, 0) / calibration_matrix(2, 2);
  *skew = calibration_matrix(0, 1) / calibration_matrix(2, 2);
  *aspect_ratio = calibration_matrix(1, 1) / calibration_matrix(0, 0);
  *principal_point_x = calibration_matrix(0, 2) / calibration_matrix(2, 2);
  *principal_point_y = calibration_matrix(1, 2) / calibration_matrix(2, 2);
}

bool DecomposeProjectionMatrix(const Matrix3x4d pmatrix,
                               Matrix3d* calibration_matrix,
                               Matrix3d* rotation_matrix,
                               Vector3d* position) {
  RQDecomposition<Matrix3d> rq(pmatrix.block<3, 3>(0, 0));

  *rotation_matrix = ProjectToRotationMatrix(rq.matrixQ());

  const double k_det = rq.matrixR().determinant();
  if (k_det == 0) {
    return false;
  }

  Matrix3d& kmatrix = *calibration_matrix;
  if (k_det > 0) {
    kmatrix = rq.matrixR();
  } else {
    kmatrix = -rq.matrixR();
  }

  // Fix the matrix such that all internal parameters are greater than 0.
  for (int i = 0; i < 3; ++i) {
    if (kmatrix(i, i) < 0) {
      kmatrix.col(i) *= -1.0;
      rotation_matrix->row(i) *= -1.0;
    }
  }

  // Solve for t.
  const Vector3d t =
      kmatrix.triangularView<Eigen::Upper>().solve(pmatrix.col(3));

  // c = - R' * t, and flip the sign according to k_det;
  if (k_det > 0) {
    *position = - rotation_matrix->transpose() * t;
  } else {
    *position = rotation_matrix->transpose() * t;
  }

  return true;
}

bool ComposeProjectionMatrix(const Matrix3d& calibration_matrix,
                             const Matrix3d& rotation_matrix,
                             const Vector3d& position,
                             Matrix3x4d* pmatrix) {
  pmatrix->block<3, 3>(0, 0) = rotation_matrix;

  pmatrix->col(3) = - (pmatrix->block<3, 3>(0, 0) *  position);
  *pmatrix = calibration_matrix * (*pmatrix);
  return true;
}

}  // namespace theia
