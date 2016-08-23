// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_UNDISTORT_IMAGE_H_
#define THEIA_SFM_UNDISTORT_IMAGE_H_

#include "theia/sfm/camera/camera.h"
#include "theia/image/image.h"

#include "theia/sfm/camera/camera_intrinsics_model.h"

namespace theia {

// Given an image with lens distortion distortion described by the camera
// parameters, undistort the image to produce an image free of lens
// distortion. This is accomplished by mapping distorted pixels to undistorted
// pixels, then cropping the image so that no unmapped pixels will be present
// (i.e. no black pixels are present). Cropping and undistorting the image in
// this way will cause the principal point and focal length to change, and so
// the undistorted camera is returned.
//
// The implementation of this method was inspired by the library
// COLMAP: https://colmap.github.io/
bool UndistortImage(const Camera& distorted_camera,
                    const FloatImage& distorted_image,
                    Camera* undistorted_camera,
                    FloatImage* undistorted_image);

}  // namespace theia

#endif  // THEIA_SFM_UNDISTORT_IMAGE_H_
