// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#include "theia/sfm/transformation/align_point_clouds.h"

#include <glog/logging.h>
#include <Eigen/Dense>

namespace theia {

void AlignPointCloudsUmeyama(const std::vector<Eigen::Vector3d>& left,
                             const std::vector<Eigen::Vector3d>& right,
                             Eigen::Matrix3d* rotation,
                             Eigen::Vector3d* translation, double* scale) {
  std::vector<double> weights(left.size(), 1.0);
  AlignPointCloudsUmeyamaWithWeights(left, right, weights, rotation,
                                     translation, scale);
}

void AlignPointCloudsUmeyamaWithWeights(
    const std::vector<Eigen::Vector3d>& left,
    const std::vector<Eigen::Vector3d>& right,
    const std::vector<double>& weights, Eigen::Matrix3d* rotation,
    Eigen::Vector3d* translation, double* scale) {
  CHECK_EQ(left.size(), right.size());
  CHECK_EQ(left.size(), weights.size());
  CHECK_NOTNULL(rotation);
  CHECK_NOTNULL(translation);
  CHECK_NOTNULL(scale);

  // Fill outputs (useful when it fails)
  *scale = 1.0;
  *translation = Eigen::Vector3d::Zero();
  *rotation = Eigen::Matrix3d::Identity();

  const size_t num_points = left.size();
  Eigen::Map<const Eigen::Matrix<double, 3, Eigen::Dynamic> > left_points(
      left[0].data(), 3, num_points);
  Eigen::Map<const Eigen::Matrix<double, 3, Eigen::Dynamic> > right_points(
      right[0].data(), 3, num_points);

  Eigen::Vector3d left_centroid, right_centroid;
  left_centroid.setZero();
  right_centroid.setZero();
  double weights_sum = 0.0;
  for (size_t i = 0; i < num_points; i++) {
    CHECK_GE(weights[i], 0)
        << "The point weight must be greater or equal to zero.";
    weights_sum += weights[i];
    left_centroid += left[i] * weights[i];
    right_centroid += right[i] * weights[i];
  }
  // Check if the sum is valid
  CHECK_GT(weights_sum, 0) << "The sum of weights must be greater than zero.";

  left_centroid /= weights_sum;
  right_centroid /= weights_sum;

  double sigma = 0.0;
  for (size_t i = 0; i < num_points; i++) {
    sigma += (left[i] - left_centroid).squaredNorm() * weights[i];
  }
  sigma /= weights_sum;

  // Calculate cross correlation matrix based on the points shifted about the
  // centroid.
  Eigen::Matrix3d cross_correlation = Eigen::Matrix3d::Zero();
  for (int i = 0; i < num_points; i++) {
    cross_correlation += weights[i] * (left_points.col(i) - left_centroid) *
                         (right_points.col(i) - right_centroid).transpose();
  }
  cross_correlation /= weights_sum;

  // Compute SVD decomposition of the cross correlation.
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(
      cross_correlation.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);

  const Eigen::Matrix3d& umatrix = svd.matrixU();
  const Eigen::Matrix3d& vtmatrix = svd.matrixV().transpose();
  const Eigen::Vector3d& singular_values = svd.singularValues();

  const double det = umatrix.determinant() * vtmatrix.determinant();
  Eigen::Matrix3d s = Eigen::Matrix3d::Identity();
  s(2, 2) = det > 0 ? 1 : -1;

  *scale =
      (singular_values(0) + singular_values(1) + s(2, 2) * singular_values(2)) /
      sigma;
  *rotation = umatrix * s * vtmatrix;
  *translation = right_centroid - (*scale) * (*rotation) * left_centroid;
}

}  // namespace theia
