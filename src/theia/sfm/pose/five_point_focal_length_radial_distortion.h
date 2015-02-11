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

#ifndef THEIA_SFM_POSE_FIVE_POINT_FOCAL_LENGTH_RADIAL_DISTORTION_H_
#define THEIA_SFM_POSE_FIVE_POINT_FOCAL_LENGTH_RADIAL_DISTORTION_H_

#include <Eigen/Core>
#include <vector>

#include "theia/alignment/alignment.h"

namespace theia {

// Description: Compute the absolute pose, focal length, and radial distortion
// of a camera using five 3D-to-2D correspondences from the ICCV paper: "Real
// time solution to the absolute pose problem with unknown radial distortion and
// focal length" by Kukelova et. al.
//
// The method solves for the projection matrix (up to scale) by using a cross
// product constraint on the standard projection equation. This allows for
// simple solution to the first two rows of the projection matrix, and the third
// row (which contains the focal length and distortion parameters) can then be
// solved with SVD on the remaining constraint equations from the first row of
// the projection matrix. See the paper for more details.
//
// Input:
//   feature_positions: feature positions (must be 5 features).
//   world_points: 5 3-vectors with corresponding 3D world points
//   num_radial_distortion_params: The number of radial distortion paramters to
//     solve for. Must be 1, 2, or 3. Experiments by cmsweeney showed that the
//     solution method is unstable when no radial distortion is solved for.
//   solutions: Camera projection matrices (that encapsulate focal
//     length) and radial distortion parameters. Thes are only up to scale.
// Output: true if success, false if not.

bool FivePointFocalLengthRadialDistortion(
    const std::vector<Eigen::Vector2d>& feature_positions,
    const std::vector<Eigen::Vector3d>& world_points,
    const int num_radial_distortion_params,
    std::vector<Eigen::Matrix<double, 3, 4> >* projection_matrices,
    std::vector<std::vector<double> >* radial_distortions);


}  // namespace theia

#endif  // THEIA_SFM_POSE_FIVE_POINT_FOCAL_LENGTH_RADIAL_DISTORTION_H_
