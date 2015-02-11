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

#ifndef THEIA_SFM_POSE_PERSPECTIVE_THREE_POINT_H_
#define THEIA_SFM_POSE_PERSPECTIVE_THREE_POINT_H_

#include <Eigen/Core>
#include <vector>

namespace theia {
// Computes camera pose using the three point algorithm and returns all possible
// solutions (up to 4). Follows steps from the paper "A Novel Parameterization
// of the Perspective-Three-Point Problem for a direct computation of Absolute
// Camera position and Orientation" by Kneip et. al.
//
// Params:
//   feature_point: Feature points corresponding to model points. NOTE: these
//     points should be calibrated with the camera intrinsics as opposed to raw
//     pixel coordinates.
//   points_3d: 3D location of features. Must correspond to the image_ray
//     of the same index.
//   solution_rotations: the rotation matrix of the candidate solutions
//   solution_translation: the translation of the candidate solutions
// Return: the number of poses computed.
bool PoseFromThreePoints(const Eigen::Vector2d feature_point[3],
                         const Eigen::Vector3d world_point[3],
                         std::vector<Eigen::Matrix3d>* solution_rotations,
                         std::vector<Eigen::Vector3d>* solution_translations);

}  // namespace theia

#endif  // THEIA_SFM_POSE_PERSPECTIVE_THREE_POINT_H_
