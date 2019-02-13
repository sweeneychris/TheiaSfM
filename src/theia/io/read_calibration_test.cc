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
    " \"radial_distortion_coeffs\" : [0.1, 0.1], "
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_2.jpg\","
    " \"focal_length\" : 350,"
    " \"principal_point_x\" : 240,"
    " \"principal_point_y\" : 240,"
    " \"aspect_ratio\" : 1.5,"
    " \"skew\" : 0.25,"
    " \"radial_distortion_coeffs\" : [0.1], "
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_3.jpg\","
    " \"principal_point_x\" : 240,"
    " \"principal_point_y\" : 240,"
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_4.jpg\","
    " \"focal_length\" : 300,"
    " \"width\" : 480,"
    " \"height\" : 480,"
    " \"principal_point_x\" : 240,"
    " \"principal_point_y\" : 240,"
    " \"aspect_ratio\" : 1.0,"
    " \"skew\" : 0.0,"
    " \"radial_distortion_coeffs\" : [0.1, 0.1, 0.01], "
    " \"tangential_distortion_coeffs\" : [0.05, 0.05], "
    " \"camera_intrinsics_type\" : \"PINHOLE_RADIAL_TANGENTIAL\""
    "}} "
    "]}";

TEST(ReadCalibrationTest, ParseIntrinsicPriorsFromJsonStr) {
  VLOG(3) << "Input JSON: \n" << kCameraIntrinsicsPriorsJson;
  std::unordered_map<std::string, CameraIntrinsicsPrior> view_to_prior;
  EXPECT_TRUE(ExtractCameraIntrinsicPriorsFromJson(
      kCameraIntrinsicsPriorsJson, &view_to_prior));
  // Check first camera.
  EXPECT_TRUE(view_to_prior.find("view_1.jpg") != view_to_prior.end());
  const CameraIntrinsicsPrior prior1 = view_to_prior["view_1.jpg"];
  EXPECT_TRUE(prior1.principal_point.is_set);
  EXPECT_NEAR(prior1.image_width / 2.0, prior1.principal_point.value[0], 1e-6);
  EXPECT_NEAR(prior1.image_height / 2.0, prior1.principal_point.value[1], 1e-6);
  EXPECT_TRUE(prior1.focal_length.is_set);
  EXPECT_NEAR(prior1.focal_length.value[0], 300, 1e-6);
  EXPECT_TRUE(prior1.aspect_ratio.is_set);
  EXPECT_NEAR(prior1.aspect_ratio.value[0], 1.0, 1e-6);
  EXPECT_TRUE(prior1.skew.is_set);
  EXPECT_NEAR(prior1.skew.value[0], 0.0, 1e-6);
  EXPECT_TRUE(prior1.radial_distortion.is_set);
  EXPECT_NEAR(prior1.radial_distortion.value[0], 0.1, 1e-6);
  EXPECT_NEAR(prior1.radial_distortion.value[1], 0.1, 1e-6);
  EXPECT_TRUE(view_to_prior.find("view_1.jpg") != view_to_prior.end());

  // Check second camera.
  const CameraIntrinsicsPrior prior2 = view_to_prior["view_2.jpg"];
  EXPECT_TRUE(prior2.principal_point.is_set);
  EXPECT_NEAR(prior2.image_width / 2.0, prior2.principal_point.value[0], 1e-6);
  EXPECT_NEAR(prior2.image_height / 2.0, prior2.principal_point.value[1], 1e-6);
  EXPECT_TRUE(prior2.focal_length.is_set);
  EXPECT_NEAR(prior2.focal_length.value[0], 350, 1e-6);
  EXPECT_TRUE(prior2.aspect_ratio.is_set);
  EXPECT_NEAR(prior2.aspect_ratio.value[0], 1.5, 1e-6);
  EXPECT_TRUE(prior2.skew.is_set);
  EXPECT_NEAR(prior2.skew.value[0], 0.25, 1e-6);
  EXPECT_TRUE(prior2.radial_distortion.is_set);
  EXPECT_NEAR(prior2.radial_distortion.value[0], 0.1, 1e-6);
  EXPECT_NEAR(prior2.radial_distortion.value[1], 0.0, 1e-6);
  EXPECT_TRUE(view_to_prior.find("view_2.jpg") != view_to_prior.end());

  // Check third camera.
  const CameraIntrinsicsPrior prior3 = view_to_prior["view_3.jpg"];
  EXPECT_TRUE(prior3.principal_point.is_set);
  EXPECT_NEAR(prior3.image_width / 2.0, prior3.principal_point.value[0], 1e-6);
  EXPECT_NEAR(prior3.image_height / 2.0, prior3.principal_point.value[1], 1e-6);
  EXPECT_FALSE(prior3.aspect_ratio.is_set);
  EXPECT_FALSE(prior3.skew.is_set);
  EXPECT_FALSE(prior3.radial_distortion.is_set);
  EXPECT_FALSE(prior3.focal_length.is_set);
  EXPECT_TRUE(view_to_prior.find("view_3.jpg") != view_to_prior.end());

  // Fourth camera is pinhole-radial-tangential.
  const CameraIntrinsicsPrior prior4 = view_to_prior["view_4.jpg"];
  EXPECT_TRUE(prior4.principal_point.is_set);
  EXPECT_NEAR(prior4.image_width / 2.0, prior4.principal_point.value[0], 1e-6);
  EXPECT_NEAR(prior4.image_height / 2.0, prior4.principal_point.value[1], 1e-6);
  EXPECT_TRUE(prior4.focal_length.is_set);
  EXPECT_NEAR(prior4.focal_length.value[0], 300, 1e-6);
  EXPECT_TRUE(prior4.aspect_ratio.is_set);
  EXPECT_NEAR(prior4.aspect_ratio.value[0], 1.0, 1e-6);
  EXPECT_TRUE(prior4.skew.is_set);
  EXPECT_NEAR(prior4.skew.value[0], 0.0, 1e-6);
  EXPECT_TRUE(prior4.radial_distortion.is_set);
  EXPECT_NEAR(prior4.radial_distortion.value[0], 0.1, 1e-6);
  EXPECT_NEAR(prior4.radial_distortion.value[1], 0.1, 1e-6);
  EXPECT_NEAR(prior4.radial_distortion.value[2], 0.01, 1e-6);
  EXPECT_TRUE(prior4.tangential_distortion.is_set);
  EXPECT_NEAR(prior4.tangential_distortion.value[0], 0.05, 1e-6);
  EXPECT_NEAR(prior4.tangential_distortion.value[1], 0.05, 1e-6);
  EXPECT_TRUE(view_to_prior.find("view_4.jpg") != view_to_prior.end());
}

}  // namespace
}  // namespace theia
