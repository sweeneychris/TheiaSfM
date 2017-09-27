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

#include "theia/sfm/pose/relative_pose_from_two_points_with_known_rotation.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <glog/logging.h>
#include <vector>

namespace theia {

// We use the epipolar constraint to solve for the unknown relative position:
//
//   p1^t * [t]_x * p2 = 0
//
// By stacking this constraint for both correspondences, we can constraint the 2
// degrees of freedom in our unknown relative position and solve for the
// solution easily.
bool RelativePoseFromTwoPointsWithKnownRotation(
    const Eigen::Vector2d rotated_features1[2],
    const Eigen::Vector2d rotated_features2[2],
    Eigen::Vector3d* relative_position2) {
  CHECK_NOTNULL(relative_position2);

  // The epipolar constraint can be rewritten in terms of the rotated features p
  // an d q as:
  //
  // p2 * q1 * t3 - p2 * t1 - q1 * t2 + q2 * t1 - p1 * q2 * t3 + p1 * t2 = 0
  //
  // Note that p3 and q3 are both 1.0, so they do not appear in this
  // equation. These constraints are stacked together to form a linear system
  // whose null space reveals the relative position.
  Eigen::Matrix<double, 2, 3> epipolar_constraint;
  for (int i = 0; i < 2; i++) {
    epipolar_constraint(i, 0) =
        -rotated_features1[i].y() + rotated_features2[i].y();
    epipolar_constraint(i, 1) =
        -rotated_features2[i].x() + rotated_features1[i].x();
    epipolar_constraint(i, 2) =
        rotated_features1[i].y() * rotated_features2[i].x() -
        rotated_features1[i].x() * rotated_features2[i].y();
  }

  // Compute the null space of this matrix.
  Eigen::FullPivLU<Eigen::Matrix<double, 2, 3> > linear_solver(
      epipolar_constraint);
  // If the null space is greater than 1 dimension then the linear system is
  // ill-conditioned or otherwise degenerate and the relative position cannot be
  // solved for.
  if (linear_solver.dimensionOfKernel() != 1) {
    return false;
  }
  // Extract the relative position as the null space kernel.
  *relative_position2 = linear_solver.kernel().normalized();
  return true;
}

}  // namespace theia
