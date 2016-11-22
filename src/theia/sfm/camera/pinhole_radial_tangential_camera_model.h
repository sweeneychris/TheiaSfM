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

#ifndef THEIA_SFM_CAMERA_PINHOLE_RADIAL_TANGENTIAL_CAMERA_MODEL_H_
#define THEIA_SFM_CAMERA_PINHOLE_RADIAL_TANGENTIAL_CAMERA_MODEL_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <ceres/ceres.h>
#include <stdint.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/types.h"

namespace theia {

// This class contains the full camera pose information including extrinsic
// parameters as well as intrinsic parameters. Extrinsic parameters include the
// camera orientation (as angle-axis) and position, and intrinsic parameters
// include focal length, aspect ratio, skew, principal points, radial
// distortion, and tangential distortion. Methods are provided for common
// transformations and projections.
//
// Intrinsics of the camera are modeled such that:
//
//  K = [f     s     px]
//      [0   f * a   py]
//      [0     0      1]
//
// where f = focal length, px and py is the principal point, s = skew, and
// a = aspect ratio.
//
// Extrinsic parametser transform the homogeneous 3D point X to the image point
// p such that:
//   p = R * (X[0..2] / X[3] - C);
//   p = p[0,1] / p[2];
//   r = p[0] * p[0] + p[1] * p[1];
//   d = 1 + k1 * r + k2 * r * r + k3 * r * r * r;
//   td[0] = 2 * t1 * p[0] * p[1] + t2 * (r + 2 * p[0] * p[0])
//   td[0] = 2 * t2 * p[0] * p[1] + t1 * (r + 2 * p[1] * p[1])
//   p[0] = p[0] * d + td[0];
//   p[1] = p[1] * d + td[1];
//   p = K * p;
//
//  where R = orientation, C = camera position, k1 k2 k3 are the radial
//  distortion parameters, and t1 t2 are the tangential distortion parameters.
class PinholeRadialTangentialCameraModel : public CameraIntrinsicsModel {
 public:
  PinholeRadialTangentialCameraModel();
  ~PinholeRadialTangentialCameraModel() {}

  static const int kIntrinsicsSize = 10;

  enum InternalParametersIndex{
    FOCAL_LENGTH = 0,
    ASPECT_RATIO = 1,
    SKEW = 2,
    PRINCIPAL_POINT_X = 3,
    PRINCIPAL_POINT_Y = 4,
    RADIAL_DISTORTION_1 = 5,
    RADIAL_DISTORTION_2 = 6,
    RADIAL_DISTORTION_3 = 7,
    TANGENTIAL_DISTORTION_1 = 8,
    TANGENTIAL_DISTORTION_2 = 9,
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

  void SetSkew(const double skew);
  double Skew() const;

  void SetRadialDistortion(const double radial_distortion_1,
                           const double radial_distortion_2,
                           const double radial_distortion_3);
  double RadialDistortion1() const;
  double RadialDistortion2() const;
  double RadialDistortion3() const;

  void SetTangentialDistortion(const double tangential_distortion_1,
                               const double tangential_distortion_2);
  double TangentialDistortion1() const;
  double TangentialDistortion2() const;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    ar(cereal::base_class<CameraIntrinsicsModel>(this));
  }
};

template <typename T>
void PinholeRadialTangentialCameraModel::CameraToPixelCoordinates(
    const T* intrinsic_parameters, const T* point, T* pixel) {
  // Get normalized pixel projection at image plane depth = 1.
  const T& depth = point[2];
  const T normalized_pixel[2] = { point[0] / depth,
                                  point[1] / depth };

  // Apply lens distortion.
  T distorted_pixel[2];
  PinholeRadialTangentialCameraModel::DistortPoint(
      intrinsic_parameters, normalized_pixel, distorted_pixel);

  // Apply calibration parameters to transform normalized units into pixels.
  const T& focal_length =
      intrinsic_parameters[PinholeRadialTangentialCameraModel::FOCAL_LENGTH];
  const T& skew =
      intrinsic_parameters[PinholeRadialTangentialCameraModel::SKEW];
  const T& aspect_ratio =
      intrinsic_parameters[PinholeRadialTangentialCameraModel::ASPECT_RATIO];
  const T& principal_point_x = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::PRINCIPAL_POINT_X];
  const T& principal_point_y = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::PRINCIPAL_POINT_Y];

  pixel[0] = focal_length * distorted_pixel[0] + skew * distorted_pixel[1] +
             principal_point_x;
  pixel[1] = focal_length * aspect_ratio * distorted_pixel[1] +
             principal_point_y;
}

template <typename T>
void PinholeRadialTangentialCameraModel::PixelToCameraCoordinates(
    const T* intrinsic_parameters, const T* pixel, T* point) {
  const T& focal_length =
      intrinsic_parameters[PinholeRadialTangentialCameraModel::FOCAL_LENGTH];
  const T& aspect_ratio =
      intrinsic_parameters[PinholeRadialTangentialCameraModel::ASPECT_RATIO];
  const T& focal_length_y = focal_length * aspect_ratio;
  const T& skew =
      intrinsic_parameters[PinholeRadialTangentialCameraModel::SKEW];
  const T& principal_point_x = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::PRINCIPAL_POINT_X];
  const T& principal_point_y = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::PRINCIPAL_POINT_Y];

  // Normalize the y coordinate first so that the skew may be applied
  // appropriately.
  T distorted_point[2];
  distorted_point[1] = (pixel[1] - principal_point_y) / focal_length_y;
  distorted_point[0] =
      (pixel[0] - principal_point_x - distorted_point[1] * skew) / focal_length;

  // Undo the lens distortion.
  T undistorted_point[2];
  PinholeRadialTangentialCameraModel::UndistortPoint(intrinsic_parameters,
                                                     distorted_point, point);
  point[2] = T(1.0);
}

template <typename T>
void PinholeRadialTangentialCameraModel::DistortPoint(
    const T* intrinsic_parameters,
    const T* undistorted_point,
    T* distorted_point) {
  const T& radial_distortion1 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_1];
  const T& radial_distortion2 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_2];
  const T& radial_distortion3 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_3];
  const T& tangential_distortion1 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_1];
  const T& tangential_distortion2 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_2];

  // Compute the radius of the undistorted point.
  const T r_sq = undistorted_point[0] * undistorted_point[0] +
                 undistorted_point[1] * undistorted_point[1];
  // Find
  const T rd = T(1.0) + radial_distortion1 * r_sq +
               radial_distortion2 * r_sq * r_sq +
               radial_distortion3 * r_sq * r_sq * r_sq;

  // Tangential distortion x = T2 * (r_sq + 2 * x * x) + 2 * T1 * x * y.
  const T tangential_x =
      tangential_distortion2 *
          (r_sq + T(2.0) * undistorted_point[0] * undistorted_point[0]) +
      T(2.0) * tangential_distortion1 * undistorted_point[0] *
          undistorted_point[1];

  // Tangential distortion y = T1 * (r_sq + 2 * y * y) + 2 * T2 * x * y.
  const T tangential_y =
      tangential_distortion1 *
          (r_sq + T(2.0) * undistorted_point[1] * undistorted_point[1]) +
      T(2.0) * tangential_distortion2 * undistorted_point[0] *
          undistorted_point[1];

  // Apply the distortion to the undistorted point.
  distorted_point[0] = undistorted_point[0] * rd + tangential_x;
  distorted_point[1] = undistorted_point[1] * rd + tangential_y;
}

template <typename T>
void PinholeRadialTangentialCameraModel::UndistortPoint(
    const T* intrinsic_parameters,
    const T* distorted_point,
    T* undistorted_point) {
  const int kNumUndistortionIterations = 100;
  const T kUndistortionEpsilon = T(1e-10);

  const T& radial_distortion1 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_1];
  const T& radial_distortion2 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_2];
  const T& radial_distortion3 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_3];
  const T& tangential_distortion1 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_1];
  const T& tangential_distortion2 = intrinsic_parameters
      [PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_2];

  T prev_undistorted_point[2];
  undistorted_point[0] = distorted_point[0];
  undistorted_point[1] = distorted_point[1];
  for (size_t i = 0; i < kNumUndistortionIterations; ++i) {
    prev_undistorted_point[0] = undistorted_point[0];
    prev_undistorted_point[1] = undistorted_point[1];

    // Compute the radius of the undistorted point.
    const T r_sq = undistorted_point[0] * undistorted_point[0] +
                   undistorted_point[1] * undistorted_point[1];
    // Find
    const T rd = T(1.0) + radial_distortion1 * r_sq +
                 radial_distortion2 * r_sq * r_sq +
                 radial_distortion3 * r_sq * r_sq * r_sq;

    // Tangential distortion x = T2 * (r_sq + 2 * x * x) + 2 * T1 * x * y.
    const T tangential_x =
        tangential_distortion2 *
            (r_sq + T(2.0) * undistorted_point[0] * undistorted_point[0]) +
        T(2.0) * tangential_distortion1 * undistorted_point[0] *
            undistorted_point[1];

    // Tangential distortion y = T1 * (r_sq + 2 * y * y) + 2 * T2 * x * y.
    const T tangential_y =
        tangential_distortion1 *
            (r_sq + T(2.0) * undistorted_point[1] * undistorted_point[1]) +
        T(2.0) * tangential_distortion2 * undistorted_point[0] *
            undistorted_point[1];

    // We know that the distorted point = d * undistorted point, so we can solve
    // for a better estimate of the undistorted point by taking the inverse of
    // this equation: undistorted_point = distorted_point / d.
    undistorted_point[0] = (distorted_point[0] - tangential_x) / rd;
    undistorted_point[1] = (distorted_point[1] - tangential_y) / rd;

    // Repeat until convergence.
    if (ceres::abs(undistorted_point[0] - prev_undistorted_point[0]) <
            kUndistortionEpsilon &&
        ceres::abs(undistorted_point[1] - prev_undistorted_point[1]) <
            kUndistortionEpsilon) {
      break;
    }
  }
}

}  // namespace theia

#include <cereal/archives/portable_binary.hpp>

CEREAL_CLASS_VERSION(theia::PinholeRadialTangentialCameraModel, 0)
// Register the polymorphic relationship for serialization.
CEREAL_REGISTER_TYPE(theia::PinholeRadialTangentialCameraModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(theia::CameraIntrinsicsModel,
                                     theia::PinholeRadialTangentialCameraModel)

#endif  // THEIA_SFM_CAMERA_PINHOLE_RADIAL_TANGENTIAL_CAMERA_MODEL_H_
