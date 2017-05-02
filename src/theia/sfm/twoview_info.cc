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

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <glog/logging.h>
#include <algorithm>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/twoview_info.h"

namespace theia {

void SwapCameras(TwoViewInfo* twoview_info) {
  CHECK_NE(twoview_info->focal_length_1, 0.0);
  CHECK_NE(twoview_info->focal_length_2, 0.0);

  // Swap the focal lengths.
  std::swap(twoview_info->focal_length_1, twoview_info->focal_length_2);

  // Invert the translation.
  Eigen::Matrix3d rotation_mat;
  ceres::AngleAxisToRotationMatrix(
      twoview_info->rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(rotation_mat.data()));
  twoview_info->position_2 = -rotation_mat * twoview_info->position_2;

  // Invert the rotation.
  twoview_info->rotation_2 *= -1.0;
}

// Constructs the TwoViewInfo object from two cameras
void TwoViewInfoFromTwoCameras(const Camera& camera1,
                               const Camera& camera2,
                               TwoViewInfo* info) {
  CHECK_NOTNULL(info);

  // Fetch the "world-to-camera" rotation matrices for convenience.
  const Eigen::Matrix3d rotation1 = camera1.GetOrientationAsRotationMatrix();
  const Eigen::Matrix3d rotation2 = camera2.GetOrientationAsRotationMatrix();

  // Construct the two view info such that camera1 is the reference view.
  info->focal_length_1 = camera1.FocalLength();
  info->focal_length_2 = camera2.FocalLength();

  // The relative rotation of camera2 is: R_12 = R2 * R1^t. This is constructed
  // such that rotations map from the coordinate system of camera 1 into the
  // coordinate system of camera2.
  const Eigen::Matrix3d relative_rotation = rotation2 * rotation1.transpose();
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(relative_rotation.data()),
      info->rotation_2.data());

  // Compute the position of camera 2 in the coordinate system of camera 1 using
  // the standard projection equation:
  //    X' = R * (X - c)
  // which yields:
  //    c2' = R1 * (c2 - c1).
  info->position_2 =
      rotation1 * (camera2.GetPosition() - camera1.GetPosition());
  // Scale the relative position to be a unit-length vector.
  info->position_2.normalize();
}

}  // namespace theia
