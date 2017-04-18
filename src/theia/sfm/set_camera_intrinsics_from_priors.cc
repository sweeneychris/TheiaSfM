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

#include "theia/sfm/set_camera_intrinsics_from_priors.h"

#include <glog/logging.h>

#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {
namespace {

// Choose the representative camera for the intrinsics group, and initialize its
// intrinsics.
ViewId InitializeRepresentativeCameraInGroup(
    const std::unordered_set<ViewId>& views_in_intrinsics_group,
    Reconstruction* reconstruction) {
  ViewId representative_view_id = kInvalidViewId;
  for (const ViewId view_id : views_in_intrinsics_group) {
    View* view = reconstruction->MutableView(view_id);
    // Skip this view if it does not exist in the reconstruction.
    if (view == nullptr) {
      continue;
    }

    // Set the representative view to the current view. It should not matter
    // much which view we choose if none of them are already estimated. We
    // simply assign the representative views here to avoid any conditionals.
    representative_view_id = view_id;

    // If the views is already estimated, use its intrinsics as the
    // initialization point for the group's intrinsics.
    if (view->IsEstimated()) {
      break;
    }
  }

  return representative_view_id;
}

}  // namespace

// Sets the camera intrinsics from the CameraIntrinsicsPrior of each view. Views
// that do not have a focal length prior will set a value corresponding to a
// median viewing angle. Principal points that are not provided by the priors
// are simply initialized as half of the corresponding image size dimension.
void SetCameraIntrinsicsFromPriors(Reconstruction* reconstruction) {
  // Set the camera intrinsics from the priors one group at a time.
  const std::unordered_set<CameraIntrinsicsGroupId>
      camera_intrinsics_group_ids = reconstruction->CameraIntrinsicsGroupIds();
  for (const CameraIntrinsicsGroupId intrinsics_group_id :
       camera_intrinsics_group_ids) {
    // Get all views in this camera intrinsics group.
    const std::unordered_set<ViewId> views_in_intrinsics_group =
        reconstruction->GetViewsInCameraIntrinsicGroup(intrinsics_group_id);

    // We choose a "representative view" for the intrinsics group. We set the
    // intrinsics for this view from the priors, then set the intrinsics of all
    // other views in the same intrinsics group to point to the representative
    // intrinsics. Since shared_ptrs are used, the shared intrinsics remain
    // alive until all cameras in the group go out of context.
    const ViewId representative_view_id = InitializeRepresentativeCameraInGroup(
        views_in_intrinsics_group, reconstruction);
    // Set the representative camera intrinsics.
    View* representative_view =
      reconstruction->MutableView(representative_view_id);
    Camera* representative_camera = representative_view->MutableCamera();
    representative_camera->SetFromCameraIntrinsicsPriors(
        representative_view->CameraIntrinsicsPrior());

    CHECK_NOTNULL(representative_camera);
    CHECK_NOTNULL(representative_camera->CameraIntrinsics().get());

    // Set all intrinsics for this group.
    for (const ViewId view_in_intrinsics_group : views_in_intrinsics_group) {
      View* view = reconstruction->MutableView(view_in_intrinsics_group);
      // Skip this view if it does not exist in the reconstruction or it is the
      // representative view.
      if (view == nullptr ||
          representative_view_id == view_in_intrinsics_group) {
        continue;
      }

      // Set the view's intrinsics to point to the shared intrinsics. This
      // includes estimated views who may have estimated intrinsics parameters
      // that are different than the shared intrinsics.
      CHECK_NOTNULL(representative_camera->CameraIntrinsics().get());
      view->MutableCamera()->MutableCameraIntrinsics() =
        representative_camera->CameraIntrinsics();
    }
  }
}

}  // namespace theia
