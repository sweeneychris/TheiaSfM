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

#ifndef THEIA_SFM_SET_CAMERA_INTRINSICS_FROM_PRIORS_H_
#define THEIA_SFM_SET_CAMERA_INTRINSICS_FROM_PRIORS_H_

#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

class Camera;
class Reconstruction;
class View;

// Sets the camera intrinsics from the CameraIntrinsicsPrior.. If the prior does
// not have a focal length prior and focal_length_from_median_fov is true, then
// the focal length will be set to a value corresponding to a median viewing
// angle. Principal points that are not provided by the priors are simply
// initialized as half of the corresponding image size dimension.
void SetCameraIntrinsicsFromPriors(const CameraIntrinsicsPrior& prior,
                                   const bool focal_length_from_median_fov,
                                   const Camera* camera);

// Sets the camera intrinsics for every view in the reconstruction.
void SetCameraIntrinsicsFromPriors(Reconstruction* reconstruction);

}  // namespace theia

#endif  // THEIA_SFM_SET_CAMERA_INTRINSICS_FROM_PRIORS_H_
