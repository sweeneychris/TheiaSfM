// Copyright (C) 2014 The Regents of the University of California (Regents)
// and Google, Inc. All rights reserved.
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
//     * Neither the name of The Regents or University of California, Google,
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu), John Flynn (jflynn@google.com)

#ifndef THEIA_SFM_POSE_FOUR_POINT_RELATIVE_POSE_PARTIAL_ROTATION_H_
#define THEIA_SFM_POSE_FOUR_POINT_RELATIVE_POSE_PARTIAL_ROTATION_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/alignment/alignment.h"

namespace theia {
// Solves for the limited transformation between correspondences from two sets
// of image rays. These image rays may come from generalized cameras or from
// different cameras that are in the same known coordinate system.
// The transformation is limited in that it only solves for a single
// rotation around a known axis. The translation is solved exactly if at least
// one of the rays in one of the cameras has a different origin than the other
// rays in that camera
//
// Details can be found in the paper "Solving Relative Pose with a
// Partially Known Rotation is a Quadratic Eigenvalue Problem" by Sweeney et al
// (3DV 2014).
//
// This implementation is intended to form the core of a RANSAC routine, and as
// such has an optimized interface for this use case.
//
// Computes the limited pose between the two sets of image rays, that may come
// from different cameras or generalized cameras that do not have a single
// center of projection. Places the rotation and translation solutions in
// soln_rotations and soln_translations. The translations are computed up to
// unit length. There are at most 6 solutions. The rotations and translations
// are defined such that the ray in image one are transformed according to:
//
//     ray_in_image_2 = Q * ray_in_image_1 + t
//
// The computed rotations are guaranteed to be rotations around the passed
// axis only.
void FourPointRelativePosePartialRotation(
    const Eigen::Vector3d& rotation_axis,
    const Eigen::Vector3d image_one_ray_directions[4],
    const Eigen::Vector3d image_one_ray_origins[4],
    const Eigen::Vector3d image_two_ray_directions[4],
    const Eigen::Vector3d image_two_ray_origins[4],
    std::vector<Eigen::Quaterniond>* soln_rotations,
    std::vector<Eigen::Vector3d>* soln_translations);

}  // namespace theia

#endif  // THEIA_SFM_POSE_FOUR_POINT_RELATIVE_POSE_PARTIAL_ROTATION_H_
