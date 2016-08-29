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

#include "theia/sfm/camera/camera_intrinsics_model.h"

#include <Eigen/Core>
#include <memory>

#include "theia/sfm/camera/fisheye_camera_model.h"
#include "theia/sfm/camera/fov_camera_model.h"
#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/camera/pinhole_radial_tangential_camera_model.h"

namespace theia {

// Creates a camera model object based on the model type.
std::unique_ptr<CameraIntrinsicsModel>
CameraIntrinsicsModel::Create(const CameraIntrinsicsModelType& camera_type) {
  std::unique_ptr<CameraIntrinsicsModel> camera_model;
  switch (camera_type) {
    case CameraIntrinsicsModelType::PINHOLE:
      camera_model.reset(new PinholeCameraModel());
      break;
    case CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL:
      camera_model.reset(new PinholeRadialTangentialCameraModel());
      break;
    case CameraIntrinsicsModelType::FISHEYE:
      camera_model.reset(new FisheyeCameraModel());
      break;
    case CameraIntrinsicsModelType::FOV:
      camera_model.reset(new FOVCameraModel());
      break;
    default:
      LOG(FATAL) << "Invalid Camera model chosen.";
      break;
  }

  return camera_model;
}

CameraIntrinsicsModel& CameraIntrinsicsModel::operator=(
    const CameraIntrinsicsModel& camera) {
  CHECK(this->Type() == camera.Type())
      << "Cannot assign camera intrinsics model of type "
      << static_cast<int>(camera.Type()) << " to a camera model of type "
      << static_cast<int>(this->Type())
      << " See camera_intrinsics_model.h for more information.";

  parameters_ = camera.parameters_;

  return *this;
}

// This macro is used to call static methods within derived classes. Rather than
// repeating this switch statment everywhere, we simply use this macro. Each use
// of this macro requries that the macro CAMERA_MODEL_CASE_BODY is defined to be
// the lines that are executed for all cases.
#ifndef CAMERA_MODEL_SWITCH_STATEMENT
#define CAMERA_MODEL_SWITCH_STATEMENT                                     \
  switch (this->Type()) {                                                 \
    CAMERA_MODEL_CASE(PINHOLE, PinholeCameraModel)                        \
    CAMERA_MODEL_CASE(PINHOLE_RADIAL_TANGENTIAL,                          \
                      PinholeRadialTangentialCameraModel)                 \
    CAMERA_MODEL_CASE(FISHEYE, FisheyeCameraModel)                        \
    CAMERA_MODEL_CASE(FOV, FOVCameraModel)                                \
    default:                                                              \
      LOG(FATAL)                                                          \
          << "Invalid camera type. Please see camera_intrinsics_model.h " \
             "for a list of valid camera models.";                        \
      break;                                                              \
  }
#endif  // CAMERA_MODEL_SWITCH_STATEMENT

#ifndef CAMERA_MODEL_CASE
#define CAMERA_MODEL_CASE(CameraModelType, CameraModel) \
  case CameraIntrinsicsModelType::CameraModelType:      \
    CAMERA_MODEL_CASE_BODY(CameraModel)                 \
    break;
#endif  // CAMERA_MODEL_CASE

Eigen::Vector2d CameraIntrinsicsModel::CameraToImageCoordinates(
    const Eigen::Vector3d& point) const {
  Eigen::Vector2d pixel;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel)                         \
  CameraModel::CameraToPixelCoordinates(parameters(), point.data(), \
                                        pixel.data());

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return pixel;
}

Eigen::Vector3d CameraIntrinsicsModel::ImageToCameraCoordinates(
    const Eigen::Vector2d& pixel) const {
  Eigen::Vector3d point;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel)                         \
  CameraModel::PixelToCameraCoordinates(parameters(), pixel.data(), \
                                        point.data());

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return point;
}

Eigen::Vector2d CameraIntrinsicsModel::DistortPoint(
    const Eigen::Vector2d& undistorted_point) const {
  Eigen::Vector2d distorted_point;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel)                         \
  CameraModel::DistortPoint(parameters(), undistorted_point.data(), \
                            distorted_point.data());

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return distorted_point;
}

Eigen::Vector2d CameraIntrinsicsModel::UndistortPoint(
    const Eigen::Vector2d& distorted_point) const {
  Eigen::Vector2d undistorted_point;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel)                         \
  CameraModel::UndistortPoint(parameters(), distorted_point.data(), \
                              undistorted_point.data());

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return undistorted_point;
}

void CameraIntrinsicsModel::SetFocalLength(const double focal_length) {
  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
  CHECK_GT(focal_length, 0.0)
      << "Invalid focal length value. Focal length must be greater than 0.0.";
#define CAMERA_MODEL_CASE_BODY(CameraModel) \
  SetParameter(CameraModel::FOCAL_LENGTH, focal_length);

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY
}

double CameraIntrinsicsModel::FocalLength() const {
  double focal_length = -1.0;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel) \
  focal_length = GetParameter(CameraModel::FOCAL_LENGTH);

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return focal_length;
}

void CameraIntrinsicsModel::SetPrincipalPoint(const double principal_point_x,
    const double principal_point_y) {
  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel)                        \
  SetParameter(CameraModel::PRINCIPAL_POINT_X, principal_point_x); \
  SetParameter(CameraModel::PRINCIPAL_POINT_Y, principal_point_y);

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY
}

double CameraIntrinsicsModel::PrincipalPointX() const {
  double principal_point_x = 0.0;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel) \
  principal_point_x = GetParameter(CameraModel::PRINCIPAL_POINT_X);

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return principal_point_x;
}

double CameraIntrinsicsModel::PrincipalPointY() const {
  double principal_point_y = 0.0;

  // Define the functions that we want to execute in every case of the switch
  // statement. CameraModel will be filled in with the appropriate derived
  // class.
#define CAMERA_MODEL_CASE_BODY(CameraModel) \
  principal_point_y = GetParameter(CameraModel::PRINCIPAL_POINT_Y);

  // Execute the switch statement.
  CAMERA_MODEL_SWITCH_STATEMENT

#undef CAMERA_MODEL_CASE_BODY

  return principal_point_y;
}

void CameraIntrinsicsModel::SetParameter(const int parameter_index,
                                         const double parameter_value) {
  DCHECK_GE(parameter_index, 0);
  DCHECK_LT(parameter_index, parameters_.size());
  parameters_[parameter_index] = parameter_value;
}

const double CameraIntrinsicsModel::GetParameter(
    const int parameter_index) const {
  DCHECK_GE(parameter_index, 0);
  DCHECK_LT(parameter_index, parameters_.size());
  return parameters_[parameter_index];
}

const double* CameraIntrinsicsModel::parameters() const {
  return parameters_.data();
}

double* CameraIntrinsicsModel::mutable_parameters() {
  return parameters_.data();
}

#undef CAMERA_MODEL_SWITCH_STATEMENT

}  // namespace theia
