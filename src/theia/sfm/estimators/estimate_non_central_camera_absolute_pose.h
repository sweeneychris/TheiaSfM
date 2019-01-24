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

#ifndef THEIA_SFM_ESTIMATORS_ESTIMATE_NON_CENTRAL_CAMERA_ABSOLUTE_POSE_H_
#define THEIA_SFM_ESTIMATORS_ESTIMATE_NON_CENTRAL_CAMERA_ABSOLUTE_POSE_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "theia/sfm/create_and_initialize_ransac_variant.h"

namespace theia {

struct FeatureCorrespondence2D3D;
struct NonCentralCameraFeatureCorrespondence;
struct RansacParameters;
struct RansacSummary;

struct NonCentralCameraAbsolutePose {
  Eigen::Quaterniond rotation;
  Eigen::Vector3d translation;
};

// Estimates the pose of a non-central camera via Upnp and a RANSAC variant
// (e.g., Ransac, Prosac, etc.). Correspondences comprise ray origins or camera
// centers, ray directions or bearing vectors, and world points or 3D points.
// Returns true if a pose could be succesfully estimated, and false otherwise.
// The quality of the result depends on the quality of the input data.
//
// Params:
//   ransac_params:  Ransac parameters.
//   ransac_type:  Ransac variant, e.g., Ransac, Prosac, etc..
//   correspondences:  Ray origins, ray directions, and world points.
//   estimated_poses:  The estimated poses for the camera rig.
//   ransac_summary:  A summary of the Ransac iterations.
bool EstimateNonCentralCameraAbsolutePose(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<NonCentralCameraFeatureCorrespondence>& correspondences,
    NonCentralCameraAbsolutePose* estimated_poses,
    RansacSummary* ransac_summary);

// Estimates the calibrated absolute pose via Upnp, a non-central-camera pose
// estimator, and a RANSAC variant (e.g., Ransac, Prosac, etc.). Correspondences
// must be normalized by the camera intrinsics. Returns true if a pose could be
// succesfully estimated, and false otherwise. The quality of the result depends
// on the quality of the input data.
//
// Params:
//   ransac_params:  Ransac parameters.
//   ransac_type:  Ransac variant, e.g., Ransac, Prosac, etc..
//   normalized_correspondences:  2D-3D or 2D feature to 3D point.
//     correspondences. Note that the 2D features must be normalized by the
//     camera intrinsics.
//   estimated_poses:  The estimated poses for the camera rig.
//   ransac_summary:  A summary of the Ransac iterations.
bool EstimateNonCentralCameraAbsolutePose(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<FeatureCorrespondence2D3D>& normalized_correspondences,
    NonCentralCameraAbsolutePose* estimated_poses,
    RansacSummary* ransac_summary);

}  // namespace theia

#endif  // THEIA_SFM_ESTIMATORS_ESTIMATE_NON_CENTRAL_CAMERA_ABSOLUTE_POSE_H_
