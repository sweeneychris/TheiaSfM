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

#ifndef THEIA_SFM_CAMERA_PINHOLE_CAMERA_MODEL_H_
#define THEIA_SFM_CAMERA_PINHOLE_CAMERA_MODEL_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.hpp>
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
class PinholeCameraModel : public CameraIntrinsicsModel {
 public:
  PinholeCameraModel();
  ~PinholeCameraModel() {}

  int NumParameters() const;

  // ---------------------------- Helper methods ---------------------------- //
  // Returns the calibration matrix in the form specified above.
  void GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const;

  // Returns the camera model type of the object.
  CameraIntrinsicsModelType Type() const;

  // Set the intrinsic camera parameters from the priors.
  void SetFromCameraIntrinsicsPriors(
      const CameraIntrinsicsPrior& prior);

  // Returns the indices of the parameters that will be optimized during bundle
  // adjustment.
  std::vector<int> GetSubsetFromOptimizeIntrinsicsType(
      const OptimizeIntrinsicsType& intrinsics_to_optimize);

  // Projects the homogeneous 3D point into the image plane and undistorts the
  // point according to the radial distortion parameters. The function returns
  // the depth of the point so that points that project behind the camera (i.e.,
  // negative depth) can be determined. Points at infinity return a depth of
  // infinity.
  Eigen::Vector2d CameraToImageCoordinates(
      const Eigen::Vector3d& point) const;

  // Converts image pixel coordinates to normalized coordinates in the camera
  // coordinates by removing the effect of camera intrinsics/calibration.
  Eigen::Vector3d ImageToCameraCoordinates(
      const Eigen::Vector2d& pixel) const;

  // ----------------------- Getter and Setter methods ---------------------- //
  void SetFocalLength(const double focal_length);
  double FocalLength() const;

  void SetPrincipalPoint(const double principal_point_x,
                                 const double principal_point_y);
  double PrincipalPointX() const;
  double PrincipalPointY() const;

  void SetAspectRatio(const double aspect_ratio);
  double AspectRatio() const;

  void SetSkew(const double skew);
  double Skew() const;

  void SetRadialDistortion(const double radial_distortion_1,
                           const double radial_distortion_2);
  double RadialDistortion1() const;
  double RadialDistortion2() const;

  // Directly get and set the parameters directly. Each derived class will
  // define a set of indices for the intrinsic parameters as a public enum.
  void SetParameter(const int parameter_index, double parameter_value);
  const double GetParameter(const int parameter_index) const;

  const double* parameters() const;
  double* mutable_parameters();

  enum InternalParametersIndex{
    FOCAL_LENGTH = 0,
    ASPECT_RATIO = 1,
    SKEW = 2,
    PRINCIPAL_POINT_X = 3,
    PRINCIPAL_POINT_Y = 4,
    RADIAL_DISTORTION_1 = 5,
    RADIAL_DISTORTION_2 = 6
  };

  static const int kIntrinsicsSize = 7;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    ar(cereal::binary_data(camera_parameters_,
                           sizeof(double) * kIntrinsicsSize));
  }

  double camera_parameters_[kIntrinsicsSize];
};

}  // namespace theia

#include <cereal/archives/portable_binary.hpp>

CEREAL_CLASS_VERSION(theia::PinholeCameraModel, 0)
// Register the polymorphic relationship for serialization.
CEREAL_REGISTER_TYPE(theia::PinholeCameraModel)
CEREAL_REGISTER_POLYMORPHIC_RELATION(theia::CameraIntrinsicsModel,
                                     theia::PinholeCameraModel)

#endif  // THEIA_SFM_CAMERA_PINHOLE_CAMERA_MODEL_H_
