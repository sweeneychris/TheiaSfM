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

#include "theia/sfm/camera/camera.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <algorithm>

#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/camera/projection_matrix_utils.h"
#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

Camera::Camera() {
  // Set rotation and position to zero (i.e. identity).
  Map<Matrix<double, 1, 6> >(mutable_extrinsics()).setZero();
  camera_intrinsics_.reset(new PinholeCameraModel());

  image_size_[0] = 0;
  image_size_[1] = 0;
}

Camera::Camera(const CameraIntrinsicsModelType& camera_type) {
  // Set rotation and position to zero (i.e. identity).
  Map<Matrix<double, 1, 6> >(mutable_extrinsics()).setZero();
  camera_intrinsics_ = CameraIntrinsicsModel::Create(camera_type);

  image_size_[0] = 0;
  image_size_[1] = 0;
}

Camera::Camera(const Camera& camera) {
  // Copy the extrinsics.
  std::copy(camera.parameters(),
            camera.parameters() + kExtrinsicsSize,
            camera_parameters_);

  camera_intrinsics_ =
      CameraIntrinsicsModel::Create(camera.GetCameraIntrinsicsModelType());
  const CameraIntrinsicsModel& other_camera_intrinsics =
      camera.CameraIntrinsics();
  std::copy(other_camera_intrinsics.parameters(),
            other_camera_intrinsics.parameters() +
                other_camera_intrinsics.NumParameters(),
            camera_intrinsics_->mutable_parameters());

  image_size_[0] = camera.image_size_[0];
  image_size_[1] = camera.image_size_[1];
}

Camera& Camera::operator=(const Camera& camera) {
  // Copy the extrinsics.
  std::copy(camera.parameters(),
            camera.parameters() + kExtrinsicsSize,
            camera_parameters_);

  camera_intrinsics_ =
      CameraIntrinsicsModel::Create(camera.GetCameraIntrinsicsModelType());
  const CameraIntrinsicsModel& other_camera_intrinsics =
      camera.CameraIntrinsics();
  std::copy(other_camera_intrinsics.parameters(),
            other_camera_intrinsics.parameters() +
                other_camera_intrinsics.NumParameters(),
            camera_intrinsics_->mutable_parameters());

  image_size_[0] = camera.image_size_[0];
  image_size_[1] = camera.image_size_[1];
  return *this;
}

void Camera::SetFromCameraIntrinsicsPriors(const CameraIntrinsicsPrior& prior) {
  camera_intrinsics_ =
      CameraIntrinsicsModel::Create(prior.camera_intrinsics_model_type);
  image_size_[0] = prior.image_width;
  image_size_[1] = prior.image_height;
  camera_intrinsics_->SetFromCameraIntrinsicsPriors(prior);
}

CameraIntrinsicsModelType Camera::GetCameraIntrinsicsModelType() const {
  return camera_intrinsics_->Type();
}

void Camera::SetCameraIntrinsicsModelType(
    const CameraIntrinsicsModelType& camera_model_type) {
  // Only reset and change the camera model if the camera intrinsics model type
  // has changed.
  if (GetCameraIntrinsicsModelType() != camera_model_type) {
    camera_intrinsics_ = CameraIntrinsicsModel::Create(camera_model_type);
  }
}

bool Camera::InitializeFromProjectionMatrix(
      const int image_width,
      const int image_height,
      const Matrix3x4d projection_matrix) {
  DCHECK_GT(image_width, 0);
  DCHECK_GT(image_height, 0);
  image_size_[0] = image_width;
  image_size_[1] = image_height;

  camera_intrinsics_.reset(new PinholeCameraModel());

  Vector3d orientation, position;
  Matrix3d calibration_matrix;
  DecomposeProjectionMatrix(projection_matrix,
                            &calibration_matrix,
                            &orientation,
                            &position);

  Map<Vector3d>(mutable_extrinsics() + ORIENTATION) = orientation;
  Map<Vector3d>(mutable_extrinsics() + POSITION) = position;

  if (calibration_matrix(0, 0) == 0 || calibration_matrix(1, 1) == 0) {
    LOG(INFO) << "Cannot set focal lengths to zero!";
    return false;
  }

  double* mutable_intrinsics = camera_intrinsics_->mutable_parameters();
  CalibrationMatrixToIntrinsics(
      calibration_matrix,
      mutable_intrinsics + PinholeCameraModel::FOCAL_LENGTH,
      mutable_intrinsics + PinholeCameraModel::SKEW,
      mutable_intrinsics + PinholeCameraModel::ASPECT_RATIO,
      mutable_intrinsics + PinholeCameraModel::PRINCIPAL_POINT_X,
      mutable_intrinsics + PinholeCameraModel::PRINCIPAL_POINT_Y);
  return true;
}

void Camera::GetProjectionMatrix(Matrix3x4d* pmatrix) const {
  Matrix3d calibration_matrix;
  camera_intrinsics_->GetCalibrationMatrix(&calibration_matrix);
  ComposeProjectionMatrix(calibration_matrix,
                          GetOrientationAsAngleAxis(),
                          GetPosition(),
                          pmatrix);
}

void Camera::GetCalibrationMatrix(Matrix3d* kmatrix) const {
  camera_intrinsics_->GetCalibrationMatrix(kmatrix);
}

double Camera::ProjectPoint(const Vector4d& point, Vector2d* pixel) const {
  Eigen::Vector3d adjusted_point = point.head<3>() - point[3] * GetPosition();
  Eigen::Vector3d rotated_point;
  ceres::AngleAxisRotatePoint(extrinsics() + ORIENTATION,
                              adjusted_point.data(),
                              rotated_point.data());
  *pixel = camera_intrinsics_->CameraToImageCoordinates(rotated_point);

  return rotated_point[2] / point[3];
}

Vector3d Camera::PixelToUnitDepthRay(const Vector2d& pixel) const {
  // Remove the effect of calibration.
  const Vector3d undistorted_point = PixelToNormalizedCoordinates(pixel);

  // Apply rotation.
  const Matrix3d& rotation = GetOrientationAsRotationMatrix();
  const Vector3d direction = rotation.transpose() * undistorted_point;
  return direction;
}

Vector3d Camera::PixelToNormalizedCoordinates(const Vector2d& pixel) const {
  return camera_intrinsics_->ImageToCameraCoordinates(pixel);
}

  // ----------------------- Getter and Setter methods ---------------------- //
void Camera::SetPosition(const Vector3d& position) {
  Map<Vector3d>(mutable_extrinsics() + POSITION) = position;
}

Vector3d Camera::GetPosition() const {
  return Map<const Vector3d>(extrinsics() + POSITION);
}

void Camera::SetOrientationFromRotationMatrix(const Matrix3d& rotation) {
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(rotation.data()),
      mutable_extrinsics() + ORIENTATION);
}

void Camera::SetOrientationFromAngleAxis(const Vector3d& angle_axis) {
  Map<Vector3d>(mutable_extrinsics() + ORIENTATION) = angle_axis;
}

Matrix3d Camera::GetOrientationAsRotationMatrix() const {
  Matrix3d rotation;
  ceres::AngleAxisToRotationMatrix(
      extrinsics() + ORIENTATION,
      ceres::ColumnMajorAdapter3x3(rotation.data()));
  return rotation;
}

Vector3d Camera::GetOrientationAsAngleAxis() const {
  return Map<const Vector3d>(extrinsics() + ORIENTATION);
}

void Camera::SetFocalLength(const double focal_length) {
  camera_intrinsics_->SetFocalLength(focal_length);
}

double Camera::FocalLength() const {
  return camera_intrinsics_->FocalLength();
}

void Camera::SetPrincipalPoint(const double principal_point_x,
                               const double principal_point_y) {
  camera_intrinsics_->SetPrincipalPoint(principal_point_x,
                                        principal_point_y);
}

double Camera::PrincipalPointX() const {
  return camera_intrinsics_->PrincipalPointX();
}

double Camera::PrincipalPointY() const {
  return camera_intrinsics_->PrincipalPointY();
}

void Camera::SetImageSize(const int image_width, const int image_height) {
  image_size_[0] = image_width;
  image_size_[1] = image_height;
}

}  // namespace theia
