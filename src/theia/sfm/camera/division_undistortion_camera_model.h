// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#ifndef THEIA_SFM_CAMERA_DIVISION_UNDISTORTION_CAMERA_MODEL_H_
#define THEIA_SFM_CAMERA_DIVISION_UNDISTORTION_CAMERA_MODEL_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <ceres/jet.h>
#include <stdint.h>
#include <vector>

#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/types.h"

namespace theia {

// This class contains the camera intrinsic information for a 1-parameter radial
// distortion "division model" as proposed in:
//
//   "Simultaneous linear estimation of multiple view geometry and lens
//   distortion" by Andrew Fitzgibbon, CVPR 2001.
//
// This model is capable of handling a wide range of radial distortions, and is
// particularly simple to implement. The division formulation also lends itself
// well to minimal pose solvers due to the simpler formulation.
//
// NOTE: The value of the radial distortion parameter is typically quite small,
// but also depends on the image resolution of the camera. csweeney found some
// papers to briefly mention that scaling the image coordinates to be between
// [-1, 1] before applying distortion can improve numeric stability. However, we
// find that this can be avoided if care is taken to initialize the radial
// distortion parameters properly.
class DivisionUndistortionCameraModel : public CameraIntrinsicsModel {
 public:
  DivisionUndistortionCameraModel();
  ~DivisionUndistortionCameraModel() {}

  static const int kIntrinsicsSize = 5;

  enum InternalParametersIndex {
    FOCAL_LENGTH = 0,
    ASPECT_RATIO = 1,
    PRINCIPAL_POINT_X = 2,
    PRINCIPAL_POINT_Y = 3,
    RADIAL_DISTORTION_1 = 4,
  };

  int NumParameters() const override;

  // Returns the camera model type of the object.
  CameraIntrinsicsModelType Type() const override;

  // Set the intrinsic camera parameters from the priors.
  void SetFromCameraIntrinsicsPriors(
      const CameraIntrinsicsPrior& prior) override;

  // Return a CameraIntrinsicsPrior that can be used to initialize a camera with
  // the same parameters with the SetFromCameraIntrinsicsPriors method.
  CameraIntrinsicsPrior CameraIntrinsicsPriorFromIntrinsics() const override;

  // Returns the indices of the parameters that will be optimized during bundle
  // adjustment.
  std::vector<int> GetSubsetFromOptimizeIntrinsicsType(
      const OptimizeIntrinsicsType& intrinsics_to_optimize) const override;

  // Returns the calibration matrix in the form specified above.
  void GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const override;

  // Prints the camera intrinsics in a human-readable format.
  void PrintIntrinsics() const override;

  // Given a point in the camera coordinate system, apply the camera intrinsics
  // (e.g., focal length, principal point, distortion) to transform the point
  // into pixel coordinates.
  template <typename T>
  static void CameraToPixelCoordinates(const T* intrinsic_parameters,
                                       const T* point,
                                       T* pixel);

  // Given a pixel in the image coordinates, remove the effects of camera
  // intrinsics parameters and lens distortion to produce a point in the camera
  // coordinate system. The point output by this method is effectively a ray in
  // the direction of the pixel in the camera coordinate system.
  template <typename T>
  static void PixelToCameraCoordinates(const T* intrinsic_parameters,
                                       const T* pixel,
                                       T* point);

  // Given an undistorted point, apply lens distortion to the point to get a
  // distorted point. The type of distortion (i.e. radial, tangential, fisheye,
  // etc.) will depend on the camera intrinsics model.
  template <typename T>
  static void DistortPoint(const T* intrinsic_parameters,
                           const T* undistorted_point,
                           T* distorted_point);

  // Given a distorted point, apply lens undistortion to the point to get an
  // undistorted point. The type of distortion (i.e. radial, tangential,
  // fisheye, etc.) will depend on the camera intrinsics model.
  template <typename T>
  static void UndistortPoint(const T* intrinsic_parameters,
                             const T* distorted_point,
                             T* undistorted_point);

  // ----------------------- Getter and Setter methods ---------------------- //
  void SetAspectRatio(const double aspect_ratio);
  double AspectRatio() const;

  void SetRadialDistortion(const double radial_distortion_1);
  double RadialDistortion1() const;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    ar(cereal::base_class<CameraIntrinsicsModel>(this));
  }
};

// In the division model, the distortion is applied after the focal length is
// first applied to the undistorted pixel.
template <typename T>
void DivisionUndistortionCameraModel::CameraToPixelCoordinates(
    const T* intrinsic_parameters, const T* point, T* pixel) {
  // Get normalized pixel projection at image plane depth = 1.
  const T& depth = point[2];
  const T normalized_pixel[2] = {point[0] / depth, point[1] / depth};

  // Apply calibration parameters to transform normalized units into pixels.
  const T& focal_length =
      intrinsic_parameters[DivisionUndistortionCameraModel::FOCAL_LENGTH];
  const T& aspect_ratio =
      intrinsic_parameters[DivisionUndistortionCameraModel::ASPECT_RATIO];
  const T focal_length_y = focal_length * aspect_ratio;
  const T& principal_point_x =
      intrinsic_parameters[DivisionUndistortionCameraModel::PRINCIPAL_POINT_X];
  const T& principal_point_y =
      intrinsic_parameters[DivisionUndistortionCameraModel::PRINCIPAL_POINT_Y];

  // Apply the focal length and aspect ratio.
  T undistorted_pixel[2];
  undistorted_pixel[0] = focal_length * normalized_pixel[0];
  undistorted_pixel[1] = focal_length_y * normalized_pixel[1];

  // Apply radial distortion.
  DivisionUndistortionCameraModel::DistortPoint(
      intrinsic_parameters, undistorted_pixel, pixel);

  // Add the principal point.
  pixel[0] += principal_point_x;
  pixel[1] += principal_point_y;
}

template <typename T>
void DivisionUndistortionCameraModel::PixelToCameraCoordinates(
    const T* intrinsic_parameters, const T* pixel, T* point) {
  const T& focal_length =
      intrinsic_parameters[DivisionUndistortionCameraModel::FOCAL_LENGTH];
  const T& aspect_ratio =
      intrinsic_parameters[DivisionUndistortionCameraModel::ASPECT_RATIO];
  const T focal_length_y = focal_length * aspect_ratio;
  const T& principal_point_x =
      intrinsic_parameters[DivisionUndistortionCameraModel::PRINCIPAL_POINT_X];
  const T& principal_point_y =
      intrinsic_parameters[DivisionUndistortionCameraModel::PRINCIPAL_POINT_Y];

  // Normalize the y coordinate first.
  T distorted_point[2];
  distorted_point[0] = pixel[0] - principal_point_x;
  distorted_point[1] = pixel[1] - principal_point_y;

  // Undo the radial distortion.
  T undistorted_point[2];
  DivisionUndistortionCameraModel::UndistortPoint(
      intrinsic_parameters, distorted_point, point);

  point[0] /= focal_length;
  point[1] /= focal_length_y;
  point[2] = T(1.0);
}

template <typename T>
void DivisionUndistortionCameraModel::DistortPoint(
    const T* intrinsic_parameters,
    const T* undistorted_point,
    T* distorted_point) {
  static const T kVerySmallNumber = std::numeric_limits<T>::epsilon();

  // The undistorted point may be simply computed from the distorted point as:
  //
  //   x_u = x_d / (1 + k * r_d^2)               (1)
  //
  // Thus, the distorted point is given as:
  //
  //   x_d = x_u * (1 + k * r_d^2)               (2)
  //
  // If we can determine the distorted radius r_d, then we can easily determine
  // the distorted point. From Eq (1), it can be easily shown that the
  // relationship between the undistorted radius and the distorted radius is
  // similarly given by:
  //
  //   r_u = r_d / (1 + k * r_d^2)
  //
  // This may be transformed into a quadratic equation in terms of the unknown
  // distorted radius r_d:
  //
  //   r_d^2 * (k * r_u) - r_d + r_u = 0
  //
  // Solving for r_d we get a single root in the valid range of r_d > 0:
  //
  //  r_d = (1 - sqrt(1 - 4 * k * r_u^2)) / (2 * k * r_u)
  //
  // We can then plug this into Eq (2) to obtain the distorted point.
  const T r_u_sq = undistorted_point[0] * undistorted_point[0] +
                   undistorted_point[1] * undistorted_point[1];
  const T r_u = ceres::sqrt(r_u_sq);
  const T& k = intrinsic_parameters
      [DivisionUndistortionCameraModel::RADIAL_DISTORTION_1];
  const T denom = 2.0 * k * r_u;

  // If the denominator is nearly zero, use L'Hopital's rule to compute r_d by
  // taking the derivatives of the numerator and denominator.
  T r_d_sq;
  if (denom < kVerySmallNumber && denom > -kVerySmallNumber) {
    // We can directly compute r_d_sq to avoid a sqrt call.
    r_d_sq = 1.0 / (1.0 - 4.0 * k * r_u_sq);
  } else {
    const T r_d = (1.0 - ceres::sqrt(1.0 - 4.0 * k * r_u_sq)) / denom;
    r_d_sq = r_d * r_d;
  }

  // Plug in r_d into Eq (2) to obtain the distorted point.
  distorted_point[0] = undistorted_point[0] * (1.0 + k * r_d_sq);
  distorted_point[1] = undistorted_point[1] * (1.0 + k * r_d_sq);
}

template <typename T>
void DivisionUndistortionCameraModel::UndistortPoint(
    const T* intrinsic_parameters,
    const T* distorted_point,
    T* undistorted_point) {
  // The squared radius of the distorted image point.
  const T r_d_sq = distorted_point[0] * distorted_point[0] +
                   distorted_point[1] * distorted_point[1];

  // The undistorted point may be simply computed from the distorted point as:
  //
  //   x_u = x_d / (1 + k * r_d^2)
  //
  // While care must be taken to avoid dividing by a zero term, this is almost
  // never the case in practice so we ignore it here.
  const T& k = intrinsic_parameters
      [DivisionUndistortionCameraModel::RADIAL_DISTORTION_1];
  const T undistortion = 1.0 / (1.0 + k * r_d_sq);
  undistorted_point[0] = distorted_point[0] * undistortion;
  undistorted_point[1] = distorted_point[1] * undistortion;
}

}  // namespace theia

#include <cereal/archives/portable_binary.hpp>

CEREAL_CLASS_VERSION(theia::DivisionUndistortionCameraModel, 0)
// Register the polymorphic relationship for serialization.
CEREAL_REGISTER_TYPE(theia::DivisionUndistortionCameraModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(theia::CameraIntrinsicsModel,
                                     theia::DivisionUndistortionCameraModel)

#endif  // THEIA_SFM_CAMERA_DIVISION_UNDISTORTION_CAMERA_MODEL_H_
