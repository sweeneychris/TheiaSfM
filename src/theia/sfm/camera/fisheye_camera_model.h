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

#ifndef THEIA_SFM_CAMERA_FISHEYE_CAMERA_MODEL_H_
#define THEIA_SFM_CAMERA_FISHEYE_CAMERA_MODEL_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <stdint.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/types.h"

namespace theia {

// This class contains the camera intrinsic information for fisheye
// cameras. This camera model is more adept at modelling camera lenses with
// large amounts of radial distortion (such as with GoPro cameras) than the
// standard pinhole + radial distortion model. The 4-parameters lens distortion
// model is based on OpenCV's fisheye lens distortion:
//
// http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#fisheye
class FisheyeCameraModel : public CameraIntrinsicsModel {
 public:
  FisheyeCameraModel();
  ~FisheyeCameraModel() {}

  static const int kIntrinsicsSize = 9;

  enum InternalParametersIndex{
    FOCAL_LENGTH = 0,
    ASPECT_RATIO = 1,
    SKEW = 2,
    PRINCIPAL_POINT_X = 3,
    PRINCIPAL_POINT_Y = 4,
    RADIAL_DISTORTION_1 = 5,
    RADIAL_DISTORTION_2 = 6,
    RADIAL_DISTORTION_3 = 7,
    RADIAL_DISTORTION_4 = 8,
  };

  int NumParameters() const override;

  // Returns the camera model type of the object.
  CameraIntrinsicsModelType Type() const override;

  // Set the intrinsic camera parameters from the priors.
  void SetFromCameraIntrinsicsPriors(
      const CameraIntrinsicsPrior& prior) override;

  // Returns the indices of the parameters that will be optimized during bundle
  // adjustment.
  std::vector<int> GetSubsetFromOptimizeIntrinsicsType(
      const OptimizeIntrinsicsType& intrinsics_to_optimize) const override;

  // Returns the calibration matrix in the form specified above.
  void GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const override;

  // Given a point in the camera coordinate system, apply the camera intrinsics
  // (e.g., focal length, principal point, distortion) to transform the point
  // into pixel coordinates.
  //
  // NOTE: This method should transform to pixel coordinates and so lens
  // distortion should be applied.
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
                           const double radial_distortion_3,
                           const double radial_distortion_4);
  double RadialDistortion1() const;
  double RadialDistortion2() const;
  double RadialDistortion3() const;
  double RadialDistortion4() const;

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
void FisheyeCameraModel::CameraToPixelCoordinates(
    const T* intrinsic_parameters, const T* point, T* pixel) {
  // Get normalized pixel projection at image plane depth = 1.
  const T& depth = point[2];
  const T normalized_pixel[2] = { point[0] / depth,
                                  point[1] / depth };

  // Apply radial distortion.
  T distorted_pixel[2];
  FisheyeCameraModel::DistortPoint(intrinsic_parameters,
                                   normalized_pixel,
                                   distorted_pixel);

  // Apply calibration parameters to transform normalized units into pixels.
  const T& focal_length =
      intrinsic_parameters[FisheyeCameraModel::FOCAL_LENGTH];
  const T& skew = intrinsic_parameters[FisheyeCameraModel::SKEW];
  const T& aspect_ratio =
      intrinsic_parameters[FisheyeCameraModel::ASPECT_RATIO];
  const T& principal_point_x =
      intrinsic_parameters[FisheyeCameraModel::PRINCIPAL_POINT_X];
  const T& principal_point_y =
      intrinsic_parameters[FisheyeCameraModel::PRINCIPAL_POINT_Y];

  pixel[0] = focal_length * distorted_pixel[0] + skew * distorted_pixel[1] +
             principal_point_x;
  pixel[1] = focal_length * aspect_ratio * distorted_pixel[1] +
             principal_point_y;
}

template <typename T>
void FisheyeCameraModel::PixelToCameraCoordinates(const T* intrinsic_parameters,
                                                  const T* pixel,
                                                  T* point) {
  const T& focal_length =
      intrinsic_parameters[FisheyeCameraModel::FOCAL_LENGTH];
  const T& aspect_ratio =
      intrinsic_parameters[FisheyeCameraModel::ASPECT_RATIO];
  const T& focal_length_y = focal_length * aspect_ratio;
  const T& skew = intrinsic_parameters[FisheyeCameraModel::SKEW];
  const T& principal_point_x =
      intrinsic_parameters[FisheyeCameraModel::PRINCIPAL_POINT_X];
  const T& principal_point_y =
      intrinsic_parameters[FisheyeCameraModel::PRINCIPAL_POINT_Y];

  // Normalize the y coordinate first.
  T distorted_point[2];
  distorted_point[1] = (pixel[1] - principal_point_y) / focal_length_y;
  distorted_point[0] =
      (pixel[0] - principal_point_x - distorted_point[1] * skew) / focal_length;

  // Undo the radial distortion.
  T undistorted_point[2];
  FisheyeCameraModel::UndistortPoint(intrinsic_parameters,
                                     distorted_point,
                                     point);
  point[2] = T(1.0);
}

template <typename T>
void FisheyeCameraModel::DistortPoint(const T* intrinsic_parameters,
                                      const T* undistorted_point,
                                      T* distorted_point) {
  static const T kVerySmallNumber = T(1e-8);
  const T& radial_distortion1 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_1];
  const T& radial_distortion2 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_2];
  const T& radial_distortion3 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_3];
  const T& radial_distortion4 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_4];

  const T r_sq = undistorted_point[0] * undistorted_point[0] +
                 undistorted_point[1] * undistorted_point[1];
  const T r = ceres::sqrt(r_sq);

  // If the radius of the undistorted is too small then the divide by r below is
  // unstable. In this case, the point is very close to the center of distortion
  // and so we can assume there is no distortion.
  if (r < kVerySmallNumber) {
    distorted_point[0] = undistorted_point[0];
    distorted_point[1] = undistorted_point[1];
    return;
  }

  const T theta = ceres::atan2(r, T(1.0));
  const T theta_sq = theta * theta;
  const T theta_d =
      theta * (T(1.0) + radial_distortion1 * theta_sq +
               radial_distortion2 * theta_sq * theta_sq +
               radial_distortion3 * theta_sq * theta_sq * theta_sq +
               radial_distortion4 * theta_sq * theta_sq * theta_sq * theta_sq);

  distorted_point[0] = undistorted_point[0] * theta_d / r;
  distorted_point[1] = undistorted_point[1] * theta_d / r;
}

template <typename T>
void FisheyeCameraModel::UndistortPoint(const T* intrinsic_parameters,
                                        const T* distorted_point,
                                        T* undistorted_point) {
  static const T kVerySmallNumber = T(1e-8);
  const int kNumUndistortionIterations = 100;
  const T kUndistortionEpsilon = 1e-10;

  // Get convenient references to the lens distortion camera parameters.
  const T& radial_distortion1 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_1];
  const T& radial_distortion2 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_2];
  const T& radial_distortion3 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_3];
  const T& radial_distortion4 =
      intrinsic_parameters[FisheyeCameraModel::RADIAL_DISTORTION_4];

  T prev_undistorted_point[2];
  undistorted_point[0] = distorted_point[0];
  undistorted_point[1] = distorted_point[1];
  for (size_t i = 0; i < kNumUndistortionIterations; ++i) {
    prev_undistorted_point[0] = undistorted_point[0];
    prev_undistorted_point[1] = undistorted_point[1];

    // Compute an estimate of the radius of the undistorted point from the
    // center of distortion (which is assumed to be the principal point).
    const T r_sq = undistorted_point[0] * undistorted_point[0] +
                   undistorted_point[1] * undistorted_point[1];
    const T r = ceres::sqrt(r_sq);

    // If the radius of the undistorted is too small then the divide by r below
    // is unstable. In this case, the point is very close to the center of
    // distortion and so we can assume there is no distortion.
    if (r < kVerySmallNumber) {
      undistorted_point[0] = distorted_point[0];
      undistorted_point[1] = distorted_point[1];
      return;
    }

    const T theta = ceres::atan2(r, T(1.0));
    const T theta_sq = theta * theta;

    // Compute the distortion factor.
    const T theta_d =
        theta *
        (T(1.0) + radial_distortion1 * theta_sq +
         radial_distortion2 * theta_sq * theta_sq +
         radial_distortion3 * theta_sq * theta_sq * theta_sq +
         radial_distortion4 * theta_sq * theta_sq * theta_sq * theta_sq);

    // We know that the distorted point = theta_d / r * undistorted point, so we
    // can solve for a better estimate of the undistorted point by taking the
    // inverse of this equation:
    //   undistorted_point = r * distorted_point / theta_d.
    undistorted_point[0] = r * distorted_point[0] / theta_d;
    undistorted_point[1] = r * distorted_point[1] / theta_d;

    // Repeat until convergence.
    if (std::abs(undistorted_point[0] - prev_undistorted_point[0]) <
            kUndistortionEpsilon &&
        std::abs(undistorted_point[1] - prev_undistorted_point[1]) <
            kUndistortionEpsilon) {
      break;
    }
  }
}

}  // namespace theia

#include <cereal/archives/portable_binary.hpp>

CEREAL_CLASS_VERSION(theia::FisheyeCameraModel, 0)
// Register the polymorphic relationship for serialization.
CEREAL_REGISTER_TYPE(theia::FisheyeCameraModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(theia::CameraIntrinsicsModel,
                                     theia::FisheyeCameraModel)

#endif  // THEIA_SFM_CAMERA_FISHEYE_CAMERA_MODEL_H_
