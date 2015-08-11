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

#include "theia/sfm/camera/radial_distortion.h"

#include <Eigen/Core>
#include <glog/logging.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "theia/math/polynomial.h"

namespace theia {

// Solves for the undistorted image points. We assume that the radial distortion
// model is:
//   r = x * x + y * y;
//   d = 1 + r * (k1 + k2 * r);
//   xp = x * d;
//   yp = y * d;
//
// Given this model, we know that:
//   xp / x  =  yp / y
// So we can rewrite r such that
//   r = (1 + yp * yp / (yp * xp)) * x * x
// and now we have r in terms of a single unknown x. Plugging this into
// xp = x * d we have a 5th order polynomial in x that we can solve with an
// iterative solver.
void RadialUndistortPoint(const Eigen::Vector2d& distorted_point,
                          const double radial_distortion1,
                          const double radial_distortion2,
                          Eigen::Vector2d* undistorted_point) {
  const double kMinRadius = 1e-5;
  const double kEpsilon = 1e-8;
  const int kMaxIter = 10;

  if (std::max(std::abs(distorted_point.x()), std::abs(distorted_point.y())) <
      kMinRadius) {
    *undistorted_point = distorted_point;
    return;
  }

  // Choose which variable we solve around to improve stability.
  if (std::abs(distorted_point.x()) > std::abs(distorted_point.y())) {
    const double point_ratio = distorted_point.y() / distorted_point.x();
    const double r_sq = 1 + point_ratio * point_ratio;
    Eigen::VectorXd quintic_polynomial(6);
    quintic_polynomial(0) = radial_distortion2 * r_sq * r_sq;
    quintic_polynomial(1) = 0;
    quintic_polynomial(2) = radial_distortion1 * r_sq;
    quintic_polynomial(3) = 0;
    quintic_polynomial(4) = 1;
    quintic_polynomial(5) = -distorted_point.x();

    undistorted_point->x() = FindRootIterativeLaguerre(
        quintic_polynomial, distorted_point.x(), kEpsilon, kMaxIter);
    undistorted_point->y() = point_ratio * undistorted_point->x();
  } else {
    const double point_ratio = distorted_point.x() / distorted_point.y();
    const double r_sq = 1 + point_ratio * point_ratio;
    Eigen::VectorXd quintic_polynomial(6);
    quintic_polynomial(0) = radial_distortion2 * r_sq * r_sq;
    quintic_polynomial(1) = 0;
    quintic_polynomial(2) = radial_distortion1 * r_sq;
    quintic_polynomial(3) = 0;
    quintic_polynomial(4) = 1;
    quintic_polynomial(5) = -distorted_point.y();

    undistorted_point->y() = FindRootIterativeLaguerre(
        quintic_polynomial, distorted_point.y(), kEpsilon, kMaxIter);
    undistorted_point->x() = point_ratio * undistorted_point->y();
  }
}

}  // namespace theia
