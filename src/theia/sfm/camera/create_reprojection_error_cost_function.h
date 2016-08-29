// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_CAMERA_CREATE_REPROJECTION_ERROR_COST_FUNCTION_H_
#define THEIA_SFM_CAMERA_CREATE_REPROJECTION_ERROR_COST_FUNCTION_H_

#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/camera/fisheye_camera_model.h"
#include "theia/sfm/camera/fov_camera_model.h"
#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/camera/pinhole_radial_tangential_camera_model.h"
#include "theia/sfm/camera/reprojection_error.h"

namespace theia {
// Create the appropriate reprojection error cost function based on the camera
// intrinsics model that is passed in. The ReprojectionError struct is templated
// on the camera intrinsics model class and so it will automatically model the
// reprojection error appropriately.
inline ceres::CostFunction* CreateReprojectionErrorCostFunction(
    const CameraIntrinsicsModelType& camera_model_type,
    const Feature& feature) {
  static const int kResidualSize = 2;
  static const int kPointSize = 4;
  // Return the appropriate reprojection error cost function based on the camera
  // model type.
  switch (camera_model_type) {
    case CameraIntrinsicsModelType::PINHOLE:
      return new ceres::AutoDiffCostFunction<
          ReprojectionError<PinholeCameraModel>, kResidualSize,
          Camera::kExtrinsicsSize, PinholeCameraModel::kIntrinsicsSize,
          kPointSize>(new ReprojectionError<PinholeCameraModel>(feature));
      break;
    case CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL:
      return new ceres::AutoDiffCostFunction<
          ReprojectionError<PinholeRadialTangentialCameraModel>, kResidualSize,
          Camera::kExtrinsicsSize,
          PinholeRadialTangentialCameraModel::kIntrinsicsSize, kPointSize>(
          new ReprojectionError<PinholeRadialTangentialCameraModel>(feature));
      break;
    case CameraIntrinsicsModelType::FISHEYE:
      return new ceres::AutoDiffCostFunction<
          ReprojectionError<FisheyeCameraModel>, kResidualSize,
          Camera::kExtrinsicsSize, FisheyeCameraModel::kIntrinsicsSize,
          kPointSize>(new ReprojectionError<FisheyeCameraModel>(feature));
      break;
    case CameraIntrinsicsModelType::FOV:
      return new ceres::AutoDiffCostFunction<
          ReprojectionError<FOVCameraModel>, kResidualSize,
          Camera::kExtrinsicsSize, FOVCameraModel::kIntrinsicsSize,
          kPointSize>(new ReprojectionError<FOVCameraModel>(feature));
      break;
    default:
      LOG(FATAL) << "Invalid camera type. Please see camera_intrinsics_model.h "
                    "for a list of valid camera models.";
      break;
  }
}

}  // namespace theia

#endif  // THEIA_SFM_CAMERA_CREATE_REPROJECTION_ERROR_COST_FUNCTION_H_
