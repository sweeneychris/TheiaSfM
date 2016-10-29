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

namespace theia {
class Camera;
class FloatImage;
class Reconstruction;

// Given an image with lens distortion distortion described by the camera
// parameters, undistort the image according to the parameters of the
// undistorted camera to produce an image free of lens distortion. This is
// accomplished by mapping distorted pixels to undistorted pixels, then cropping
// the image so that no unmapped pixels will be present (i.e. no black pixels
// are present).
//
// The implementation of this method was inspired by the library
// COLMAP: https://colmap.github.io/
bool UndistortImage(const Camera& distorted_camera,
                    const FloatImage& distorted_image,
                    const Camera& undistorted_camera,
                    FloatImage* undistorted_image);

// Create the undistorted camera by removing radial distortion parameters.
bool UndistortCamera(const Camera& distorted_camera,
                     Camera* undistorted_camera);

// Undistorts the entire reconstruction. All features in all views are
// undistorted, but only the features which would survive the crop of the
// undistorted image (see description above) are kept. This will modify the
// reconstruction in place and potentially remove observations or tracks.
bool UndistortReconstruction(Reconstruction* reconstruction);

}  // namespace theia

#endif  // THEIA_SFM_UNDISTORT_IMAGE_H_
