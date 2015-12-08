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

void SetViewCameraIntrinsicsFromPriors(View* view) {
  Camera* camera = view->MutableCamera();
  const CameraIntrinsicsPrior prior = view->CameraIntrinsicsPrior();

  // Set the image dimensions.
  camera->SetImageSize(prior.image_width, prior.image_height);

  // Set the focal length.
  if (prior.focal_length.is_set) {
    camera->SetFocalLength(prior.focal_length.value);
  } else {
    camera->SetFocalLength(1.2 * static_cast<double>(std::max(
        prior.image_width, prior.image_height)));
  }

  // Set the principal point.
  if (prior.principal_point[0].is_set && prior.principal_point[1].is_set) {
    camera->SetPrincipalPoint(prior.principal_point[0].value,
                              prior.principal_point[1].value);
  } else {
    camera->SetPrincipalPoint(prior.image_width / 2.0,
                              prior.image_height / 2.0);
  }

  // Set aspect ratio if available.
  if (prior.aspect_ratio.is_set) {
    camera->SetAspectRatio(prior.aspect_ratio.value);
  }

  // Set skew if available.
  if (prior.skew.is_set) {
    camera->SetSkew(prior.skew.value);
  }

  // Set radial distortion if available.
  if (prior.radial_distortion[0].is_set &&
      prior.radial_distortion[1].is_set) {
    camera->SetRadialDistortion(prior.radial_distortion[0].value,
                                prior.radial_distortion[1].value);
  }
}

// Sets the camera intrinsics from the CameraIntrinsicsPrior of each view. Views
// that do not have a focal length prior will set a value corresponding to a
// median viewing angle. Principal points that are not provided by the priors
// are simply initialized as half of the corresponding image size dimension.
void SetCameraIntrinsicsFromPriors(Reconstruction* reconstruction) {
  const auto& view_ids = reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
    SetViewCameraIntrinsicsFromPriors(view);
  }
}

}  // namespace theia
