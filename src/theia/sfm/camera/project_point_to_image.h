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

#ifndef THEIA_SFM_CAMERA_PROJECT_POINT_TO_IMAGE_H_
#define THEIA_SFM_CAMERA_PROJECT_POINT_TO_IMAGE_H_

#include <Eigen/Core>
#include <glog/logging.h>
#include <ceres/rotation.h>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/radial_distortion.h"

namespace theia {

// Projects a homogeneous 3D point to an image by assuming a camera model
// defined by the Camera class. This function is templated so that it can be
// used by Ceres for bundle adjument in addition to the standard reprojection
// with doubles.
//
// Returns the depth of the 3D point subject to the projective scale (i.e. the
// depth of the point assuming the image plane is at a depth of 1). The depth is
// useful, for instance, to determine if the point reprojects behind the image.
//
// NOTE: The unit test for this method is included in
// theia/sfm/camera/camera_test.cc
template <typename T>
T ProjectPointToImage(const T* extrinsic_parameters,
                      const T* intrinsic_parameters,
                      const T* point,
                      T* pixel) {
  typedef Eigen::Matrix<T, 3, 1> Matrix3T;
  typedef Eigen::Map<const Matrix3T> ConstMap3T;

  // Remove the translation.
  Eigen::Matrix<T, 3, 1> adjusted_point =
      ConstMap3T(point) -
      point[3] * ConstMap3T(extrinsic_parameters + Camera::POSITION);

  // Rotate the point.
  T rotated_point[3];
  ceres::AngleAxisRotatePoint(extrinsic_parameters + Camera::ORIENTATION,
                              adjusted_point.data(),
                              rotated_point);

  // Get normalized pixel projection at image plane depth = 1.
  const T& depth = rotated_point[2];
  const T normalized_pixel[2] = { rotated_point[0] / depth,
                                  rotated_point[1] / depth };

  // Apply radial distortion.
  T distorted_pixel[2];
  RadialDistortPoint(normalized_pixel[0],
                     normalized_pixel[1],
                     intrinsic_parameters[Camera::RADIAL_DISTORTION_1],
                     intrinsic_parameters[Camera::RADIAL_DISTORTION_2],
                     distorted_pixel,
                     distorted_pixel + 1);

  // Apply calibration parameters to transform normalized units into pixels.
  const T& focal_length = intrinsic_parameters[Camera::FOCAL_LENGTH];
  const T& skew = intrinsic_parameters[Camera::SKEW];
  const T& aspect_ratio = intrinsic_parameters[Camera::ASPECT_RATIO];
  const T& principal_point_x = intrinsic_parameters[Camera::PRINCIPAL_POINT_X];
  const T& principal_point_y = intrinsic_parameters[Camera::PRINCIPAL_POINT_Y];

  pixel[0] = focal_length * distorted_pixel[0] + skew * distorted_pixel[1] +
             principal_point_x;
  pixel[1] = focal_length * aspect_ratio * distorted_pixel[1] +
             principal_point_y;

  return depth / point[3];
}

}  // namespace theia

#endif  // THEIA_SFM_CAMERA_PROJECT_POINT_TO_IMAGE_H_
