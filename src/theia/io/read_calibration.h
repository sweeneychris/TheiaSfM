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

#ifndef THEIA_IO_READ_CALIBRATION_H_
#define THEIA_IO_READ_CALIBRATION_H_

#include <string>
#include <unordered_map>

namespace theia {
struct CameraIntrinsicsPrior;

// Reads calibration data for images that are to be reconstructed. The
// calibration file should have the following JSON format:
//
//   {
//     "priors" : [
//       {"CameraIntrinsicsPrior" : {
//         "image_name" : "view_1.jpg",
//         "focal_length" : 300,
//         "width" : 480,
//         "height" : 480,
//         "principal_point" : [240, 240],
//         "aspect_ratio" : 1.0,
//         "skew" : 0.0,
//         "radial_distortion_coeffs" : [0.1, 0.01],
//         "camera_intrinsics_type" : "PINHOLE"
//        }},
//       {"CameraIntrinsicsPrior" : {
//         "image_name" : "view_2.jpg",
//         "focal_length" : 300,
//         "principal_point" : [240, 240],
//         "aspect_ratio" : 1.0,
//         "skew" : 0.0,
//         "radial_distortion_coeffs" : [0.1, 0.01],
//         "tangential_distortion_coeffs" : [0.05, 0.05],
//         "latitude" : 120.0,
//         "longitude" : 58.0,
//         "altitude" : 64,
//         "camera_intrinsics_type" : "PINHOLE_RADIAL_TANGENTIAL"
//        }}
//     ]
//   }
//
//
// When the image width and/or height are not set, Theia assumeds that the
// principal point lies at the center of the image, so the width and height are
// set to be twice those values.
//
// Notes:
//  1. See theia/sfm/camera/camera_intrinsics_model_type.h for the camera
//     intrinsic types.
//  2. A calibration file is optional and it is not required that all images
//     have calibration.
bool ReadCalibration(const std::string& calibration_file,
                     std::unordered_map<std::string, CameraIntrinsicsPrior>*
                         camera_intrinsics_prior);

// Extracts the camera intrinsics priors from a loaded JSON string, and creates
// a map between view_name to camera intrinsics priors. The function returns
// true upon successful parsing, and false otherwise.
//
// Params:
//   json_str:  The JSON file into a C string.
//   camera_intrinsics_prior:  A map from view name to camera intrinsic prior.
bool ExtractCameraIntrinsicPriorsFromJson(
    const char* json_str,
    std::unordered_map<std::string, CameraIntrinsicsPrior>*
      camera_intrinsics_prior);

}  // namespace theia

#endif  // THEIA_IO_READ_CALIBRATION_H_
