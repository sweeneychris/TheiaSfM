// Copyright (C) 2018 The Regents of the University of California (Regents).
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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#ifndef THEIA_SFM_ESTIMATORS_NON_CENTRAL_CAMERA_FEATURE_CORRESPONDENCE_H_
#define THEIA_SFM_ESTIMATORS_NON_CENTRAL_CAMERA_FEATURE_CORRESPONDENCE_H_

#include <Eigen/Core>

namespace theia {

struct NonCentralCameraFeatureCorrespondence {
  NonCentralCameraFeatureCorrespondence(const Eigen::Vector3d _ray_direction,
                                        const Eigen::Vector3d _ray_origin,
                                        const Eigen::Vector3d _world_point) {
    ray_direction = _ray_direction;
    ray_origin = _ray_origin;
    world_point = _world_point;
  }
  NonCentralCameraFeatureCorrespondence() = default;
  ~NonCentralCameraFeatureCorrespondence() = default;
  // Ray direction is a vector indicating the direction from the origin of a
  // camera to a 3D point. This is wrt to the coordinate system of the
  // non-central camera or camera rig.
  Eigen::Vector3d ray_direction;
  // Ray origin is the center of a camera. This is wrt to the coordinate system
  // of the non-central camera or camera rig.
  Eigen::Vector3d ray_origin;
  // The 3D point being observed by ray_direction.
  Eigen::Vector3d world_point;
};
  
}  // namespace theia

#endif  // THEIA_SFM_ESTIMATORS_NON_CENTRAL_CAMERA_FEATURE_CORRESPONDENCE_H_
