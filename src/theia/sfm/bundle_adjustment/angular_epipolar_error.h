// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_ANGULAR_EPIPOLAR_ERROR_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_ANGULAR_EPIPOLAR_ERROR_H_

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <Eigen/Core>

namespace theia {

struct AngularEpipolarError {
 public:
  AngularEpipolarError(const Eigen::Vector2d& feature1,
                       const Eigen::Vector2d& feature2)
      : feature1_(feature1), feature2_(feature2) {}

  template<typename T> bool operator()(const T* rotation,
                                       const T* translation,
                                       T* angular_epipolar_error) const {
    // Convert features to T.
    const Eigen::Matrix<T, 3, 1> feature1(T(feature1_(0)),
                                          T(feature1_(1)),
                                          T(1.0));
    const Eigen::Matrix<T, 3, 1> feature2(T(feature2_(0)),
                                          T(feature2_(1)),
                                          T(1.0));

    // Obtain rotation matrix.
    Eigen::Matrix<T, 3, 3> rotation_matrix;
    ceres::AngleAxisToRotationMatrix(
        rotation, ceres::ColumnMajorAdapter3x3(rotation_matrix.data()));

    // Compute values A and B (Eq. 11 in the paper).
    Eigen::Map<const Eigen::Matrix<T, 3, 1> > translation_map(translation);
    const Eigen::Matrix<T, 3, 3> translation_term =
        Eigen::Matrix<T, 3, 3>::Identity() -
        translation_map * translation_map.transpose();
    const T a =
        feature1.dot(translation_term * feature1) +
        (rotation_matrix * feature2)
            .dot(translation_term * (rotation_matrix.transpose() * feature2));
    const T b_sqrt = translation_map.dot(
        feature1.cross(rotation_matrix.transpose() * feature2));

    // Ensure the square root is real.
    const T sqrt_term = a * a / T(4.0) - b_sqrt * b_sqrt;
    if (sqrt_term < T(0.0)) {
      return false;
    }

    angular_epipolar_error[0] = a / T(2.0) - sqrt(sqrt_term);
    return true;
  }

  static ceres::CostFunction* Create(const Eigen::Vector2d& feature1,
                                     const Eigen::Vector2d& feature2) {
    static const int kParameterSize = 3;
    return new ceres::AutoDiffCostFunction<AngularEpipolarError,
                                           1,
                                           kParameterSize,
                                           kParameterSize>(
        new AngularEpipolarError(feature1, feature2));
  }

 private:
  const Eigen::Vector2d feature1_;
  const Eigen::Vector2d feature2_;
};

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_ANGULAR_EPIPOLAR_ERROR_H_
