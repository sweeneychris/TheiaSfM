// Copyright (C) 2019 The Regents of the University of California (Regents).
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

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// December 2018

#ifndef THEIA_SFM_POSE_FOUR_POINT_FOCAL_LENGTH_RADIAL_DISTORTION_H_
#define THEIA_SFM_POSE_FOUR_POINT_FOCAL_LENGTH_RADIAL_DISTORTION_H_

#include <Eigen/Core>
#include <vector>

namespace theia {

// Computes the absolute pose, focal length and one radial distortion parameter
// of a camera with the P4Pfr
// algorithm from the paper "Making Minimal Solvers for Absolute Pose
// Estimation Compact and Robust" by Larsson et al. (ICCV 2017)
// this is a Matlab port from
// http://www.maths.lth.se/matematiklth/personal/viktorl/

//
// Input:
//   featureVectors: 4 vectors with image positions.
//   worldPoints: 4 3-vectors with corresponding 3D world points
//     (each entry is a point)
// Output: bool: returns true if a solution was found
// Return: rotations - world to camera rotations (R*X+t)
//         translations - world to camera translations (R*X+t)
//         radial_distortions - radial distortion corresponding to the "Division
//         Undistortion Camera Model"
//         focal_lengths - camera's focal length
bool FourPointsPoseFocalLengthRadialDistortion(
    const std::vector<Eigen::Vector2d>& feature_vectors,
    const std::vector<Eigen::Vector3d>& world_points,
    std::vector<Eigen::Matrix3d>* rotations,
    std::vector<Eigen::Vector3d>* translations,
    std::vector<double>* radial_distortions,
    std::vector<double>* focal_lengths);

}  // namespace theia

#endif  // THEIA_SFM_POSE_FOUR_POINT_FOCAL_LENGTH_RADIAL_DISTORTION_H_
