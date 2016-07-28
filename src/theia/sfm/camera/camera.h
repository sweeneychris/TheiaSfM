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

#ifndef THEIA_SFM_CAMERA_CAMERA_H_
#define THEIA_SFM_CAMERA_CAMERA_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <glog/logging.h>
#include <stdint.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <vector>

#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/types.h"

namespace theia {

// This class contains the full camera pose information including extrinsic
// parameters as well as intrinsic parameters. Extrinsic parameters include the
// camera orientation (as angle-axis) and position, and intrinsic parameters
// include focal length, aspect ratio, skew, principal points, and (up to
// 2-parameter) radial distortion. Methods are provided for common
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
//   d = 1 + k1 * r + k2 * r * r;
//   p *= d;
//   p = K * p;
//
//  where R = orientation, C = camera position, and k1 k2 are the radial
//  distortion parameters.
class Camera {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Camera();
  explicit Camera(const CameraIntrinsicsModelType& camera_type);
  Camera(const Camera& camera);
  ~Camera() {}

  Camera& operator=(const Camera& rhs);

  // Initializes the camera intrinsic and extrinsic parameters from the
  // projection matrix by decomposing the matrix.
  //
  // NOTE: The projection matrix does not contain information about radial
  // distortion, so those parameters will need to be set separately.
  bool InitializeFromProjectionMatrix(
      const int image_width,
      const int image_height,
      const Matrix3x4d projection_matrix);

  // Returns the type corresponding to the camera intrinsics model used the
  // describe the lens of this camera.
  CameraIntrinsicsModelType GetCameraIntrinsicsModelType() const;

  // Sets the camera to use the designated camera intrinsics model type. If the
  // camera model type is the same as what is currently being used then this is
  // a no-op. Otherwise, the camera intrinsics parameters are cleared and the
  // camera intrinsics model that is used by this camera is appropriately
  // changed.
  void SetCameraIntrinsicsModelType(
      const CameraIntrinsicsModelType& camera_model_type);

  // ---------------------------- Helper methods ---------------------------- //
  // Returns the projection matrix. Does not include radial distortion.
  void GetProjectionMatrix(Matrix3x4d* pmatrix) const;

  // Returns the calibration matrix in the form specified above.
  void GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const;

  // Projects the homogeneous 3D point into the image plane and undistorts the
  // point according to the radial distortion parameters. The function returns
  // the depth of the point so that points that project behind the camera (i.e.,
  // negative depth) can be determined. Points at infinity return a depth of
  // infinity.
  double ProjectPoint(const Eigen::Vector4d& point,
                      Eigen::Vector2d* pixel) const;

  // Converts the pixel point to a ray in 3D space such that the origin of the
  // ray is at the camera center and the direction is the pixel direction
  // rotated according to the camera orientation in 3D space.
  //
  // NOTE: The depth of the ray is set to 1. This is so that we can remain
  // consistent with ProjectPoint. That is:
  //      if we have:
  //         d = ProjectPoint(X, &x);
  //         r = PixelToUnitDepthRay(x);
  //      then it will be the case that
  //         X = c + r * d;
  //    X is the 3D point
  //    x is the image projection of X
  //    c is the camera position
  //    r is the ray obtained from PixelToRay
  //    d is the depth of the 3D point with respect to the image
  Eigen::Vector3d PixelToUnitDepthRay(const Eigen::Vector2d& pixel) const;

  // Converts image pixel coordinates to normalized coordinates in the camera
  // coordinates by removing the effect of camera intrinsics/calibration. This
  // method is similar to PixelToUnitDepthRay except that it only removes the
  // effect of camera calibration and does not account for the camera pose.
  Eigen::Vector3d PixelToNormalizedCoordinates(
      const Eigen::Vector2d& pixel) const;

  // ----------------------- Getter and Setter methods ---------------------- //
  void SetPosition(const Eigen::Vector3d& position);
  Eigen::Vector3d GetPosition() const;

  void SetOrientationFromRotationMatrix(const Eigen::Matrix3d& rotation);
  void SetOrientationFromAngleAxis(const Eigen::Vector3d& angle_axis);
  Eigen::Matrix3d GetOrientationAsRotationMatrix() const;
  Eigen::Vector3d GetOrientationAsAngleAxis() const;

  void SetFocalLength(const double focal_length);
  double FocalLength() const;

  void SetPrincipalPoint(const double principal_point_x,
                         const double principal_point_y);
  double PrincipalPointX() const;
  double PrincipalPointY() const;

  void SetImageSize(const int image_width, const int image_height);
  int ImageWidth() const { return image_size_[0]; }
  int ImageHeight() const { return image_size_[1]; }

  const CameraIntrinsicsModel& CameraIntrinsics() const {
    return *camera_intrinsics_;
  }

  CameraIntrinsicsModel* MutableCameraIntrinsics() {
    return camera_intrinsics_.get();
  }

  const double* parameters() const { return camera_parameters_; }
  double* mutable_parameters() { return camera_parameters_; }

  const double* extrinsics() const { return camera_parameters_; }
  double* mutable_extrinsics() { return camera_parameters_; }

  const double* intrinsics() const { return camera_intrinsics_->parameters(); }
  double* mutable_intrinsics() {
    return camera_intrinsics_->mutable_parameters();
  }

  // Indexing for the location of parameters. Collecting the extrinsics and
  // intrinsics into a single array makes the interface to bundle adjustment
  // with Ceres much simpler.
  enum ExternalParametersIndex {
    POSITION = 0,
    ORIENTATION = 3
  };

  static const int kExtrinsicsSize = 6;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    if (version > 0) {
      ar(cereal::binary_data(camera_parameters_,
                             sizeof(double) * kExtrinsicsSize),
         camera_intrinsics_,
         cereal::binary_data(image_size_, sizeof(int) * 2));
    } else {
      CHECK(GetCameraIntrinsicsModelType() ==
            CameraIntrinsicsModelType::PINHOLE)
          << "the theia::Camera class version " << version
          << " can only serialize Pinhole cameras. Please make sure all "
             "cameras are set as pinhole cameras";
      const int num_parameters =
          kExtrinsicsSize + camera_intrinsics_->NumParameters();
      std::vector<double> parameters(num_parameters);

      // Copy the extrinsics and intrinsics into the vector.
      std::copy(camera_parameters_,
                camera_parameters_ + kExtrinsicsSize,
                parameters.data());
      std::copy(camera_intrinsics_->parameters(),
                camera_intrinsics_->parameters() +
                    camera_intrinsics_->NumParameters(),
                parameters.data() + kExtrinsicsSize);

      // I/O with the serialization.
      ar(cereal::binary_data(parameters.data(),
                             sizeof(double) * num_parameters),
         cereal::binary_data(image_size_, sizeof(int) * 2));

      // Copy the extrinsics and intrinsics back into the local variables.
      std::copy(parameters.data(),
                parameters.data() + kExtrinsicsSize,
                camera_parameters_);
      std::copy(parameters.data() + kExtrinsicsSize,
                parameters.data() + num_parameters,
                camera_intrinsics_->mutable_parameters());
    }
  }

  double camera_parameters_[kExtrinsicsSize];

  std::unique_ptr<CameraIntrinsicsModel> camera_intrinsics_;

  // The image size as width then height.
  int image_size_[2];
};

}  // namespace theia

CEREAL_CLASS_VERSION(theia::Camera, 1);

#endif  // THEIA_SFM_CAMERA_CAMERA_H_
