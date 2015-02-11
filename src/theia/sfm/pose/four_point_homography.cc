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

#include "theia/sfm/pose/four_point_homography.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <glog/logging.h>
#include <vector>

#include "theia/sfm/pose/util.h"

namespace theia {
using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;

namespace {

inline Matrix<double, 2, 9> CreateActionConstraint(const Vector2d& img1_point,
                                                   const Vector2d& img2_point) {
  Matrix<double, 2, 9> constraint;
  constraint << Eigen::RowVector3d::Zero(), -img1_point.transpose(), -1.0,
      img1_point.transpose() * img2_point.y(), img2_point.y(),
      img1_point.transpose(), 1.0, Eigen::RowVector3d::Zero(),
      -img1_point.transpose() * img2_point.x(), -img2_point.x();

  return constraint;
}

}  // namespace

// Normalized DLT method to compute the homography H that maps image points in
// image_1 to image_2 via x' = Hx (where x is in image 1 and x' is in image
// 2). The DLT algorithm implemented is from Algorithm 4.2 in Hartley and
// Zisserman (page 109).
bool FourPointHomography(const std::vector<Vector2d>& image_1_points,
                         const std::vector<Vector2d>& image_2_points,
                         Matrix3d* homography) {
  CHECK_GE(image_1_points.size(), 4);
  CHECK_EQ(image_1_points.size(), image_2_points.size());

  // Normalize the image points.
  std::vector<Vector2d> norm_image_1_points, norm_image_2_points;
  Matrix3d norm_image_1_mat, norm_image_2_mat;
  NormalizeImagePoints(image_1_points, &norm_image_1_points, &norm_image_1_mat);
  NormalizeImagePoints(image_2_points, &norm_image_2_points, &norm_image_2_mat);

  // Create the constraint matrix based on x' = Hx (Eq. 4.1 in Hartley and
  // Zisserman).
  Matrix<double, Eigen::Dynamic, 9> action_matrix(
      2 * image_1_points.size(), 9);
  for (int i = 0; i < image_1_points.size(); i++) {
    action_matrix.block<2, 9>(2 * i, 0) =
        CreateActionConstraint(norm_image_1_points[i], norm_image_2_points[i]);
  }

  const Matrix<double, 9, 1> null_vector =
      (action_matrix.transpose() * action_matrix).jacobiSvd(Eigen::ComputeFullV)
          .matrixV().rightCols<1>();

  *homography = norm_image_2_mat.inverse() *
                Eigen::Map<const Matrix3d>(null_vector.data()).transpose() *
                norm_image_1_mat;
  return true;
}

}  // namespace theia
