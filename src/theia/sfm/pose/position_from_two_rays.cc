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

#include "theia/sfm/pose/position_from_two_rays.h"

#include <Eigen/Core>
#include <Eigen/QR>
#include <glog/logging.h>
#include <vector>

namespace theia {

// We use the reprojection error constraint:
//
//   rotated_feature = [x - cx; y - cy] / (z - cz)
//
// where the 3D point is [x y z]^t and the unknown position is [cx cy cz].
// This constraint can be rearranged into a linear system:
//
//   [1  0  -u] * c = [x - u * z]
//   [0  1  -v]     = [y - v * z]
//
// where rotated_feature = [u v]. By stacking this constraint for both features
// we obtain a 4x3 linear system whose solution is the camera position.
bool PositionFromTwoRays(const Eigen::Vector2d& rotated_feature1,
                         const Eigen::Vector3d& point1,
                         const Eigen::Vector2d& rotated_feature2,
                         const Eigen::Vector3d& point2,
                         Eigen::Vector3d* position) {
  CHECK_NOTNULL(position);
  // Create the left hand side of the linear system above.
  Eigen::Matrix<double, 4, 3> lhs;
  lhs.block<2, 2>(0, 0).setIdentity();
  lhs.block<2, 2>(2, 0).setIdentity();
  lhs.block<2, 1>(0, 2) = -rotated_feature1;
  lhs.block<2, 1>(2, 2) = -rotated_feature2;

  // Create the right hand side of the linear system above.
  Eigen::Vector4d rhs;
  rhs.head<2>() = point1.head<2>() - rotated_feature1 * point1.z();
  rhs.tail<2>() = point2.head<2>() - rotated_feature2 * point2.z();

  Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 4, 3> > linear_solver(lhs);
  // If the linear solver is not well-conditioned then the position cannot be
  // reliably computed.
  if (linear_solver.rank() != 3) {
    return false;
  }

  *position = linear_solver.solve(rhs);
  return true;
}

}  // namespace theia
