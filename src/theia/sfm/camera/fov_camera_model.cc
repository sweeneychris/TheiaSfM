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

#include "theia/sfm/camera/fov_camera_model.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/camera/projection_matrix_utils.h"

namespace theia {

using Eigen::AngleAxisd;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

FOVCameraModel::FOVCameraModel() {
  parameters_.resize(kIntrinsicsSize);
  SetFocalLength(1.0);
  SetPrincipalPoint(0.0, 0.0);
  SetParameter(ASPECT_RATIO, 1.0);
  SetParameter(RADIAL_DISTORTION_1, kDefaultOmega);
}

int FOVCameraModel::NumParameters() const {return kIntrinsicsSize;}

// Returns the camera model type of the object.
CameraIntrinsicsModelType FOVCameraModel::Type() const {
  return CameraIntrinsicsModelType::FOV;
}

// Set the intrinsic camera parameters from the priors.
void FOVCameraModel::SetFromCameraIntrinsicsPriors(
    const CameraIntrinsicsPrior& prior) {
  // Set the focal length.
  if (prior.focal_length.is_set) {
    SetFocalLength(prior.focal_length.value[0]);
  } else if (prior.image_width != 0.0 && prior.image_height != 0.0) {
    SetFocalLength(1.2 * static_cast<double>(std::max(
        prior.image_width, prior.image_height)));
  }

  // Set the principal point.
  if (prior.principal_point.is_set) {
    SetPrincipalPoint(prior.principal_point.value[0],
                      prior.principal_point.value[1]);
  } else if (prior.image_width != 0.0 && prior.image_height != 0.0) {
    SetPrincipalPoint(prior.image_width / 2.0, prior.image_height / 2.0);
  }

  // Set aspect ratio if available.
  if (prior.aspect_ratio.is_set) {
    SetParameter(ASPECT_RATIO, prior.aspect_ratio.value[0]);
  }

  // Set radial distortion if available.
  if (prior.radial_distortion.is_set) {
    SetParameter(RADIAL_DISTORTION_1, prior.radial_distortion.value[0]);
  } else {
    SetParameter(RADIAL_DISTORTION_1, kDefaultOmega);
  }
}

CameraIntrinsicsPrior FOVCameraModel::CameraIntrinsicsPriorFromIntrinsics()
    const {
  CameraIntrinsicsPrior prior;
  prior.camera_intrinsics_model_type =
      CameraIntrinsicsModelTypeToString(Type());
  prior.focal_length.is_set = true;
  prior.focal_length.value[0] = FocalLength();
  prior.principal_point.is_set = true;
  prior.principal_point.value[0] = PrincipalPointX();
  prior.principal_point.value[1] = PrincipalPointY();
  prior.aspect_ratio.is_set = true;
  prior.aspect_ratio.value[0] = AspectRatio();
  prior.radial_distortion.is_set = true;
  prior.radial_distortion.value[0] = RadialDistortion1();

  return prior;
}


// Returns the indices of the parameters that will be optimized during bundle
// adjustment.
std::vector<int> FOVCameraModel::GetSubsetFromOptimizeIntrinsicsType(
    const OptimizeIntrinsicsType& intrinsics_to_optimize) const {
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
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::PRINCIPAL_POINTS) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(PRINCIPAL_POINT_X);
    constant_intrinsics.emplace_back(PRINCIPAL_POINT_Y);
  }
  if ((intrinsics_to_optimize & OptimizeIntrinsicsType::RADIAL_DISTORTION) ==
      OptimizeIntrinsicsType::NONE) {
    constant_intrinsics.emplace_back(RADIAL_DISTORTION_1);
  }
  return constant_intrinsics;
}

void FOVCameraModel::GetCalibrationMatrix(Matrix3d* kmatrix) const {
  IntrinsicsToCalibrationMatrix(parameters_[FOCAL_LENGTH],
                                0.0,
                                parameters_[ASPECT_RATIO],
                                parameters_[PRINCIPAL_POINT_X],
                                parameters_[PRINCIPAL_POINT_Y],
                                kmatrix);
}

void FOVCameraModel::PrintIntrinsics() const {
  LOG(INFO) << "Camera model type: "
            << CameraIntrinsicsModelTypeToString(Type())
            << "\nFocal length (pixels): " << FocalLength()
            << "\nPrincipal Point (px, py) = (" << PrincipalPointX() << ", "
            << PrincipalPointY() << ")"
            << "\nAspect Ratio: " << AspectRatio()
            << "\nRadialDistortion: " << RadialDistortion1();
}

// ----------------------- Getter and Setter methods ---------------------- //

void FOVCameraModel::SetAspectRatio(const double aspect_ratio) {
  parameters_[ASPECT_RATIO] = aspect_ratio;
}
double FOVCameraModel::AspectRatio() const {
  return parameters_[ASPECT_RATIO];
}

void FOVCameraModel::SetRadialDistortion(const double radial_distortion_1) {
  parameters_[RADIAL_DISTORTION_1] = radial_distortion_1;
}

double FOVCameraModel::RadialDistortion1() const {
  return parameters_[RADIAL_DISTORTION_1];
}

}  // namespace theia
