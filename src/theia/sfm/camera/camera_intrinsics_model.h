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

#ifndef THEIA_SFM_CAMERA_CAMERA_INTRINSICS_MODEL_H_
#define THEIA_SFM_CAMERA_CAMERA_INTRINSICS_MODEL_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <stdint.h>
#include <Eigen/Core>
#include <vector>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/camera/camera_intrinsics_model_type.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/types.h"

namespace theia {

// This class encapsulates the camera lens model used for projecting points in
// space onto the pixels in images. We utilize two coordinate systems:
//
//   1) Camera coordinate system: This is the 3D coordinate system centered at
//      the camera with the z-axis pointing directly forward (i.e. identity
//      orientation).
//
//   2) Image/pixel coordinate system: This 2D coordinate system has the origin
//      at the top-right of the image with the positive x-axis going towards the
//      right and the positive y-axis pointing down.
//
// The CameraIntrinsicsModel describes the mapping between camera and image
// coordinate systems. This may include parameters such as focal length,
// principal point, lens distortion, and others.
//
// To implement a new camera model please take the following steps:
//
//   1) Create a derived class from this one, and implement all of the pure
//   virtual methods and the static methods that are used for camera projection.
//
//   2) Add an enum to CameraIntrinsicsModelType and add an "else if" to the
//   Create method in this class to allow your camera model to be created.
//
//   3) Add the new class and its CameraIntrinsicsType to the
//   CAMERA_MODEL_SWITCH_STATEMENT macro in camera_intrinsics_model.cc
//
//   4) Add a switch/case in create_reprojection_error_cost_function.h to handle
//   the new camera model.
//
//   5) Create unit tests to ensure that your new camera model is functioning
//   properly!
class CameraIntrinsicsModel {
 public:
  CameraIntrinsicsModel() {}
  virtual CameraIntrinsicsModel& operator=(const CameraIntrinsicsModel& camera);
  virtual ~CameraIntrinsicsModel() {}
  // Creates a camera model object based on the model type.
  static std::unique_ptr<CameraIntrinsicsModel> Create(
      const CameraIntrinsicsModelType& camera_type);

  // All derived clases must implement this enum class to describe how the
  // camera parameters are organized. For all derived classes, the FOCAL_LENGTH,
  // PRINCIPAL_POINT_X, and PRINCIPAL_POINT_Y values must be set.
  enum InternalParametersIndex {
    FOCAL_LENGTH = -1,
    PRINCIPAL_POINT_X = -2,
    PRINCIPAL_POINT_Y = -3
  };

  // Returns the camera model type of the object.
  virtual CameraIntrinsicsModelType Type() const = 0;

  // Number of parameters that the camera model uses (i.e. the size of the
  // parameters_ container).
  virtual int NumParameters() const = 0;

  // Set the intrinsic camera parameters from the priors.
  virtual void SetFromCameraIntrinsicsPriors(
      const CameraIntrinsicsPrior& prior) = 0;

  // Returns the indices of the parameters that will be optimized during bundle
  // adjustment.
  virtual std::vector<int> GetSubsetFromOptimizeIntrinsicsType(
      const OptimizeIntrinsicsType& intrinsics_to_optimize) const = 0;

  // Returns the calibration matrix.
  virtual void GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const = 0;

  // ------------------------------------------------------------------------ //
  //  Although the methods below are static and templated they must be
  //  implemented by all derived classes. This is because these methods are used
  //  to compute the reprojection error used by Ceres for bundle adjustment.
  // ------------------------------------------------------------------------ //

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

  // Given a distorted point, apply lens undistortion to the pixel to get an
  // undistorted point. The type of distortion (i.e. radial, tangential,
  // fisheye, etc.) will depend on the camera intrinsics model.
  template <typename T>
  static void UndistortPoint(const T* intrinsic_parameters,
                             const T* distorted_point,
                             T* undistorted_point);

  // ------------------------------------------------------------------------ //
  //  All methods below this point are implemented in this base class and do
  //  not need to be implemented in the derived classes.
  // ------------------------------------------------------------------------ //

  // Projects the homogeneous 3D point in the camera coordinate system into the
  // image plane and undistorts the point according to the radial distortion
  // parameters.
  virtual Eigen::Vector2d CameraToImageCoordinates(
      const Eigen::Vector3d& point) const;

  // Converts image pixel coordinates to normalized coordinates in the camera
  // coordinates by removing the effect of camera intrinsics/calibration.
  virtual Eigen::Vector3d ImageToCameraCoordinates(
      const Eigen::Vector2d& pixel) const;

  // Apply or remove radial distortion to the given point. Points should be
  // given in *normalized* coordinates such that the effects of camera
  // intrinsics are not present.
  virtual Eigen::Vector2d DistortPoint(
      const Eigen::Vector2d& undistorted_point) const;
  virtual Eigen::Vector2d UndistortPoint(
      const Eigen::Vector2d& distorted_point) const;

  // ----------------------- Getter and Setter methods ---------------------- //
  virtual void SetFocalLength(const double focal_length);
  virtual double FocalLength() const;

  virtual void SetPrincipalPoint(const double principal_point_x,
                                 const double principal_point_y);
  virtual double PrincipalPointX() const;
  virtual double PrincipalPointY() const;

  // Directly get and set the parameters. Each derived class will define a set
  // of indices for the intrinsic parameters as a public enum.
  virtual void SetParameter(const int parameter_index,
                            const double parameter_value);
  virtual const double GetParameter(const int parameter_index) const;

  virtual const double* parameters() const;
  virtual double* mutable_parameters();

 protected:
  std::vector<double> parameters_;

  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    ar(parameters_);
  }
};

}  // namespace theia

CEREAL_CLASS_VERSION(theia::CameraIntrinsicsModel, 0)

#endif  // THEIA_SFM_CAMERA_CAMERA_INTRINSICS_MODEL_H_
