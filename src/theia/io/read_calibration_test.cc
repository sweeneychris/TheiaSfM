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
// Author: Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include <string>
#include <unordered_map>
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/io/read_calibration.h"

namespace theia {
namespace {

static const char* kCameraIntrinsicsPriorsJson =
    "{\"priors\" : ["
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_1.jpg\","
    " \"focal_length\" : 300,"
    " \"width\" : 480,"
    " \"height\" : 480,"
    " \"principal_point_x\" : 240,"
    " \"principal_point_y\" : 240,"
    " \"aspect_ratio\" : 1.0,"
    " \"skew\" : 0.0,"
    " \"radial_distortion_coeff_1\" : 0.1,"
    " \"radial_distortion_coeff_2\" : 0.01, "
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_2.jpg\","
    " \"focal_length\" : 300,"
    " \"principal_point_x\" : 240,"
    " \"principal_point_y\" : 240,"
    " \"aspect_ratio\" : 1.0,"
    " \"skew\" : 0.0,"
    " \"radial_distortion_coeff_1\" : 0.1,"
    " \"radial_distortion_coeff_2\" : 0.01, "
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_3.jpg\","
    " \"focal_length\" : 300,"
    " \"principal_point_x\" : 240,"
    " \"principal_point_y\" : 240,"
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}"
    "]}";

TEST(ReadCalibrationTest, ParseIntrinsicPriorsFromJsonStr) {
  VLOG(1) << "Input JSON: \n" << kCameraIntrinsicsPriorsJson;
  std::unordered_map<std::string, theia::CameraIntrinsicsPrior> view_to_prior;
  EXPECT_TRUE(ExtractCameraIntrinsicPriorsFromJson(
      kCameraIntrinsicsPriorsJson, &view_to_prior));
}

}  // namespace
}  // namespace theia
