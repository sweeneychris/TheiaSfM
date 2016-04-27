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

#include "theia/sfm/camera/projection_matrix_utils.h"
#include "theia/sfm/camera/project_point_to_image.h"
#include "theia/sfm/camera/radial_distortion.h"

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

  SetFocalLength(1.0);
  SetAspectRatio(1.0);
  SetSkew(0.0);
  SetPrincipalPoint(0.0, 0.0);
  SetRadialDistortion(0.0, 0.0);

  image_size_[0] = 0;
  image_size_[1] = 0;
}

bool Camera::InitializeFromProjectionMatrix(
      const int image_width,
      const int image_height,
      const Matrix3x4d projection_matrix) {
  DCHECK_GT(image_width, 0);
  DCHECK_GT(image_height, 0);
  image_size_[0] = image_width;
  image_size_[1] = image_height;

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

  CalibrationMatrixToIntrinsics(calibration_matrix,
                                mutable_intrinsics() + FOCAL_LENGTH,
                                mutable_intrinsics() + SKEW,
                                mutable_intrinsics() + ASPECT_RATIO,
                                mutable_intrinsics() + PRINCIPAL_POINT_X,
                                mutable_intrinsics() + PRINCIPAL_POINT_Y);
  return true;
}

void Camera::GetProjectionMatrix(Matrix3x4d* pmatrix) const {
  Matrix3d calibration_matrix;
  GetCalibrationMatrix(&calibration_matrix);
  ComposeProjectionMatrix(calibration_matrix,
                          GetOrientationAsAngleAxis(),
                          GetPosition(),
                          pmatrix);
}

void Camera::GetCalibrationMatrix(Matrix3d* kmatrix) const {
  IntrinsicsToCalibrationMatrix(FocalLength(),
                                Skew(),
                                AspectRatio(),
                                PrincipalPointX(),
                                PrincipalPointY(),
                                kmatrix);
}

double Camera::ProjectPoint(const Vector4d& point, Vector2d* pixel) const {
  return ProjectPointToImage(extrinsics(),
                             intrinsics(),
                             point.data(),
                             pixel->data());
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
  // First, undo the calibration.
  const double focal_length_y = FocalLength() * AspectRatio();
  const double y_normalized = (pixel[1] - PrincipalPointY()) / focal_length_y;
  const double x_normalized =
      (pixel[0] - PrincipalPointX() - y_normalized * Skew()) / FocalLength();

  // Undo radial distortion.
  const Vector2d normalized_point(x_normalized, y_normalized);
  Vector2d undistorted_pixel;
  RadialUndistortPoint(normalized_point,
                       RadialDistortion1(),
                       RadialDistortion2(),
                       &undistorted_pixel);
  const Vector3d undistorted_point = undistorted_pixel.homogeneous();
  return undistorted_point;
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
  mutable_intrinsics()[FOCAL_LENGTH] = focal_length;
}

double Camera::FocalLength() const {
  return intrinsics()[FOCAL_LENGTH];
}

void Camera::SetAspectRatio(const double aspect_ratio) {
  mutable_intrinsics()[ASPECT_RATIO] = aspect_ratio;
}
double Camera::AspectRatio() const {
  return intrinsics()[ASPECT_RATIO];
}

void Camera::SetSkew(const double skew) {
  mutable_intrinsics()[SKEW] = skew;
}

double Camera::Skew() const {
  return intrinsics()[SKEW];
}

void Camera::SetPrincipalPoint(const double principal_point_x,
                               const double principal_point_y) {
  mutable_intrinsics()[PRINCIPAL_POINT_X] = principal_point_x;
  mutable_intrinsics()[PRINCIPAL_POINT_Y] = principal_point_y;
}

double Camera::PrincipalPointX() const {
  return intrinsics()[PRINCIPAL_POINT_X];
}

double Camera::PrincipalPointY() const {
  return intrinsics()[PRINCIPAL_POINT_Y];
}

void Camera::SetRadialDistortion(const double radial_distortion_1,
                                 const double radial_distortion_2) {
  mutable_intrinsics()[RADIAL_DISTORTION_1] = radial_distortion_1;
  mutable_intrinsics()[RADIAL_DISTORTION_2] = radial_distortion_2;
}

double Camera::RadialDistortion1() const {
  return intrinsics()[RADIAL_DISTORTION_1];
}

double Camera::RadialDistortion2() const {
  return intrinsics()[RADIAL_DISTORTION_2];
}

void Camera::SetImageSize(const int image_width, const int image_height) {
  image_size_[0] = image_width;
  image_size_[1] = image_height;
}

}  // namespace theia
