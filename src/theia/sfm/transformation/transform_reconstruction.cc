// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#include "theia/sfm/transformation/transform_reconstruction.h"

#include <Eigen/Core>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {
namespace {

void TransformPoint(const Eigen::Matrix3d rotation,
                    const Eigen::Vector3d translation,
                    const double scale,
                    Eigen::Vector3d* point) {
  *point = scale * rotation * (*point) + translation;
}

void TransformCamera(const Eigen::Matrix3d rotation,
                     const Eigen::Vector3d translation,
                     const double scale,
                     Camera* camera) {
  const Eigen::Matrix3d camera_rotation =
      camera->GetOrientationAsRotationMatrix();
  camera->SetOrientationFromRotationMatrix(rotation.transpose() *
                                           camera_rotation);

  Eigen::Vector3d camera_position = camera->GetPosition();
  TransformPoint(rotation, translation, scale, &camera_position);
  camera->SetPosition(camera_position);
}

}  // namespace

// Applies the similarity transformation to the reconstruction, transforming the
// 3d points and the cameras poses appropriately.
void TransformReconstruction(const Eigen::Matrix3d rotation,
                             const Eigen::Vector3d translation,
                             const double scale,
                             Reconstruction* reconstruction) {
  const auto& view_ids = reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    View* view = reconstruction->MutableView(view_id);
    if (view->IsEstimated()) {
      TransformCamera(rotation, translation, scale, view->MutableCamera());
    }
  }

  const auto& track_ids = reconstruction->TrackIds();
  for (const TrackId track_id : track_ids) {
    Track* track = reconstruction->MutableTrack(track_id);
    if (track->IsEstimated()) {
      Eigen::Vector3d point = track->Point().hnormalized();
      TransformPoint(rotation, translation, scale, &point);
      *track->MutablePoint() = point.homogeneous();
    }
  }
}

}  // namespace theia
