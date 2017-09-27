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

#ifndef THEIA_SFM_POSE_RELATIVE_POSE_FROM_TWO_POINTS_WITH_KNOWN_ROTATION_H_
#define THEIA_SFM_POSE_RELATIVE_POSE_FROM_TWO_POINTS_WITH_KNOWN_ROTATION_H_

#include <Eigen/Core>
#include <vector>

namespace theia {

// Computes the relative position of two cameras when the relative position is
// known. This assumes that the features have been rotated into 3d space and
// homogenized. For a pixel observation [u v], the rotated feature is:
//
//   rf = hnormalize(R^t * K^{-1} * [u v 1]^t)
//
// where hnormalize is the homogeneous normalization that divides by the 3rd
// coordinate to achieve a 2d coordinate. The features in rotated_features1
// correspond to features in the first camera, and rotated_features2 corresponds
// to the second camera.
//
// The relative position of the second camera is returned. In other words, an
// estimate of the unit-norm direction: (c2 - c1) / ||c2 - c1|| is
// provided. This means that the relative position will be in the coordinate
// system defined by the rotated features.
//
// Returns true on success and false if the relative position could not be
// estimated.
bool RelativePoseFromTwoPointsWithKnownRotation(
    const Eigen::Vector2d rotated_features1[2],
    const Eigen::Vector2d rotated_features2[2],
    Eigen::Vector3d* relative_position2);

}  // namespace theia

#endif  // THEIA_SFM_POSE_RELATIVE_POSE_FROM_TWO_POINTS_WITH_KNOWN_ROTATION_H_
