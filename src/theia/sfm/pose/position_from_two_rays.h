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

#ifndef THEIA_SFM_POSE_POSITION_FROM_TWO_RAYS_H_
#define THEIA_SFM_POSE_POSITION_FROM_TWO_RAYS_H_

#include <Eigen/Core>
#include <vector>

namespace theia {

// Computes the camera position given two oriented features seen by the camera
// and their corresponding 3d points. The features should have the effect of
// intrinsics removed and be oriented in the same coordinate system as the 3d
// points. The rays are constrained by:
//
//   depth * rotated_feature = X - c
//
// Such that rotated_feature = R^t * [u v 1]^t, with R being the known
// world-to-camera rotation.
bool PositionFromTwoRays(const Eigen::Vector2d& rotated_feature1,
                         const Eigen::Vector3d& point1,
                         const Eigen::Vector2d& rotated_feature2,
                         const Eigen::Vector3d& point2,
                         Eigen::Vector3d* position);

}  // namespace theia

#endif  // THEIA_SFM_POSE_POSITION_FROM_TWO_RAYS_H_
