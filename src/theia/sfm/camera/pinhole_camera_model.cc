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

#include "theia/sfm/camera/pinhole_camera_model.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/camera/projection_matrix_utils.h"
#include "theia/sfm/camera/radial_distortion.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

PinholeCameraModel::PinholeCameraModel() {
  SetFocalLength(1.0);
  SetAspectRatio(1.0);
  SetSkew(0.0);
  SetPrincipalPoint(0.0, 0.0);
  SetRadialDistortion(0.0, 0.0);
}

int PinholeCameraModel::NumParameters() const {return kIntrinsicsSize;}

void PinholeCameraModel::GetCalibrationMatrix(Matrix3d* kmatrix) const {
  IntrinsicsToCalibrationMatrix(FocalLength(),
                                Skew(),
                                AspectRatio(),
                                PrincipalPointX(),
                                PrincipalPointY(),
                                kmatrix);
}


// Returns the camera model type of the object.
CameraIntrinsicsModelType PinholeCameraModel::Type() const {
  return CameraIntrinsicsModelType::PINHOLE;
}

// Set the intrinsic camera parameters from the priors.
void PinholeCameraModel::SetFromCameraIntrinsicsPriors(
    const CameraIntrinsicsPrior& prior) {
  // Set the focal length.
  if (prior.focal_length.is_set) {
    SetFocalLength(prior.focal_length.value);
  } else {
    SetFocalLength(1.2 * static_cast<double>(std::max(
        prior.image_width, prior.image_height)));
  }

  // Set the principal point.
  if (prior.principal_point[0].is_set && prior.principal_point[1].is_set) {
    SetPrincipalPoint(prior.principal_point[0].value,
                      prior.principal_point[1].value);
  } else {
    SetPrincipalPoint(prior.image_width / 2.0, prior.image_height / 2.0);
  }

  // Set aspect ratio if available.
  if (prior.aspect_ratio.is_set) {
    SetAspectRatio(prior.aspect_ratio.value);
  }

  // Set skew if available.
  if (prior.skew.is_set) {
    SetSkew(prior.skew.value);
  }

  // Set radial distortion if available.
  if (prior.radial_distortion[0].is_set &&
      prior.radial_distortion[1].is_set) {
    SetRadialDistortion(prior.radial_distortion[0].value,
                        prior.radial_distortion[1].value);
  }
}

// Returns the indices of the parameters that will be optimized during bundle
// adjustment.
std::vector<int> PinholeCameraModel::GetSubsetFromOptimizeIntrinsicsType(
    const OptimizeIntrinsicsType& intrinsics_to_optimize) {
  std::vector<int> constant_intrinsics;
  if (intrinsics_to_optimize == OptimizeIntrinsicsType::ALL) {
    return constant_intrinsics;
  }

  if ((intrinsics_to_optimize &
      OptimizeIntrinsicsType::FOCAL_LENGTH) == OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(FOCAL_LENGTH);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::ASPECT_RATIO) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(ASPECT_RATIO);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::SKEW) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(SKEW);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::PRINCIPAL_POINTS) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(PRINCIPAL_POINT_X);
    constant_intrinsics.emplace_back(PRINCIPAL_POINT_Y);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::RADIAL_DISTORTION) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(RADIAL_DISTORTION_1);
    constant_intrinsics.emplace_back(RADIAL_DISTORTION_2);
  }
  return constant_intrinsics;
}

Eigen::Vector2d PinholeCameraModel::CameraToImageCoordinates(
    const Eigen::Vector3d& point) const {
  const Eigen::Vector2d normalized_pixel = point.hnormalized();
  Eigen::Vector2d distorted_pixel;
  RadialDistortPoint(normalized_pixel.x(),
                     normalized_pixel.y(),
                     camera_parameters_[RADIAL_DISTORTION_1],
                     camera_parameters_[RADIAL_DISTORTION_2],
                     distorted_pixel.data(),
                     distorted_pixel.data() + 1);

  Eigen::Vector2d pixel;
  pixel.x() = FocalLength() * distorted_pixel[0] + Skew() * distorted_pixel[1] +
              PrincipalPointX();
  pixel.y() =
      FocalLength() * AspectRatio() * distorted_pixel[1] + PrincipalPointY();
  return pixel;
}

Vector3d PinholeCameraModel::ImageToCameraCoordinates(
    const Vector2d& pixel) const {
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
void PinholeCameraModel::SetFocalLength(const double focal_length) {
  camera_parameters_[FOCAL_LENGTH] = focal_length;
}

double PinholeCameraModel::FocalLength() const {
  return camera_parameters_[FOCAL_LENGTH];
}

void PinholeCameraModel::SetAspectRatio(const double aspect_ratio) {
  camera_parameters_[ASPECT_RATIO] = aspect_ratio;
}
double PinholeCameraModel::AspectRatio() const {
  return camera_parameters_[ASPECT_RATIO];
}

void PinholeCameraModel::SetSkew(const double skew) {
  camera_parameters_[SKEW] = skew;
}

double PinholeCameraModel::Skew() const {
  return camera_parameters_[SKEW];
}

void PinholeCameraModel::SetPrincipalPoint(const double principal_point_x,
                               const double principal_point_y) {
  camera_parameters_[PRINCIPAL_POINT_X] = principal_point_x;
  camera_parameters_[PRINCIPAL_POINT_Y] = principal_point_y;
}

double PinholeCameraModel::PrincipalPointX() const {
  return camera_parameters_[PRINCIPAL_POINT_X];
}

double PinholeCameraModel::PrincipalPointY() const {
  return camera_parameters_[PRINCIPAL_POINT_Y];
}

void PinholeCameraModel::SetRadialDistortion(const double radial_distortion_1,
                                 const double radial_distortion_2) {
  camera_parameters_[RADIAL_DISTORTION_1] = radial_distortion_1;
  camera_parameters_[RADIAL_DISTORTION_2] = radial_distortion_2;
}

double PinholeCameraModel::RadialDistortion1() const {
  return camera_parameters_[RADIAL_DISTORTION_1];
}

double PinholeCameraModel::RadialDistortion2() const {
  return camera_parameters_[RADIAL_DISTORTION_2];
}

// Directly get and set the parameters directly. Each derived class will
// define a set of indices for the intrinsic parameters as a public enum.
void PinholeCameraModel::SetParameter(const int parameter_index,
                                      const double parameter_value) {
  CHECK_GE(parameter_index, 0);
  CHECK_LT(parameter_index, kIntrinsicsSize);
  camera_parameters_[parameter_index] = parameter_value;
}

const double PinholeCameraModel::GetParameter(
    const int parameter_index) const {
  CHECK_GE(parameter_index, 0);
  CHECK_LT(parameter_index, kIntrinsicsSize);
  return camera_parameters_[parameter_index];
}

const double* PinholeCameraModel::parameters() const {
  return camera_parameters_;
}
double* PinholeCameraModel::mutable_parameters() { return camera_parameters_; }

}  // namespace theia
