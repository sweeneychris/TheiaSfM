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

#ifndef THEIA_SFM_CAMERA_RADIAL_DISTORTION_H_
#define THEIA_SFM_CAMERA_RADIAL_DISTORTION_H_

#include <Eigen/Core>
#include <glog/logging.h>

namespace theia {

// Applies radial distortion to the point given the 2 radial distortion
// parameters. This method is templated so that it can be used with Ceres bundle
// adjustment.
template <typename T>
void RadialDistortPoint(const T undistorted_point_x,
                        const T undistorted_point_y,
                        const T radial_distortion1,
                        const T radial_distortion2,
                        T* distorted_point_x,
                        T* distorted_point_y) {
  const T r_sq = undistorted_point_x * undistorted_point_x +
                 undistorted_point_y * undistorted_point_y;
  const T d =
      T(1.0) + r_sq * (radial_distortion1 + radial_distortion2 * r_sq);

  (*distorted_point_x) = undistorted_point_x * d;
  (*distorted_point_y) = undistorted_point_y * d;
}

// Undistorts the point given the two radial distortion parameters.
void RadialUndistortPoint(const Eigen::Vector2d& distorted_point,
                          const double radial_distortion1,
                          const double radial_distortion2,
                          Eigen::Vector2d* undistorted_point);

}  // namespace theia

#endif  // THEIA_SFM_CAMERA_RADIAL_DISTORTION_H_
