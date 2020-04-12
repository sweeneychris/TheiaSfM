// Copyright (C) 2019 The Regents of the University of California (Regents).
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
#include "theia/io/write_calibration.h"
#include "theia/util/map_util.h"

DEFINE_string(json_file, "io/write_calibration_test.json",
              "JSON writing testing file.");

namespace theia {
namespace {

std::string json_filepath = THEIA_DATA_DIR + std::string("/") + FLAGS_json_file;

static const char* kCameraIntrinsicsPriorsJson =
    "{\"priors\" : ["
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_1.jpg\","
    " \"focal_length\" : 300,"
    " \"width\" : 480,"
    " \"height\" : 480,"
    " \"principal_point\" : [240, 240],"
    " \"aspect_ratio\" : 1.0,"
    " \"skew\" : 0.0,"
    " \"radial_distortion_coeffs\" : [0.1, 0.1], "
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_2.jpg\","
    " \"focal_length\" : 350,"
    " \"principal_point\" : [240, 240],"
    " \"aspect_ratio\" : 1.5,"
    " \"skew\" : 0.25,"
    " \"radial_distortion_coeffs\" : [0.1], "
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_3.jpg\","
    " \"principal_point\" : [240, 240],"
    " \"camera_intrinsics_type\" : \"PINHOLE\""
    "}}, "
    "{\"CameraIntrinsicsPrior\": {"
    " \"image_name\" : \"view_4.jpg\","
    " \"focal_length\" : 300,"
    " \"width\" : 480,"
    " \"height\" : 480,"
    " \"principal_point\" : [240, 240],"
    " \"aspect_ratio\" : 1.0,"
    " \"skew\" : 0.0,"
    " \"radial_distortion_coeffs\" : [0.1, 0.1, 0.01], "
    " \"tangential_distortion_coeffs\" : [0.05, 0.05], "
    " \"orientation\" : [0.1, 0.1, 0.1], "
    " \"position\" : [1, 2.0, -3.0], "
    " \"latitude\" : 128.0, "
    " \"longitude\" : 256.0, "
    " \"altitude\" : 512.0, "
    " \"camera_intrinsics_type\" : \"PINHOLE_RADIAL_TANGENTIAL\""
    "}} "
    "]}";

TEST(WriteCalibrationTest, WriteAndParseIntrinsicPriors) {
  // Read the expected priors.
  std::unordered_map<std::string, CameraIntrinsicsPrior> expected_priors;
  EXPECT_TRUE(ExtractCameraIntrinsicPriorsFromJson(
      kCameraIntrinsicsPriorsJson, &expected_priors));

  // Write JSON file.
  VLOG(1) << "Writing calibration priors to: " << json_filepath;
  EXPECT_TRUE(WriteCalibration(json_filepath, expected_priors));

  // Read the JSON file and compare with expected_priors.
  std::unordered_map<std::string, CameraIntrinsicsPrior> read_priors;
  EXPECT_TRUE(ReadCalibration(json_filepath, &read_priors));

  // Compare the two priors.
  EXPECT_EQ(expected_priors.size(), read_priors.size());
  for (const auto& prior : expected_priors) {
    CameraIntrinsicsPrior* camera_prior = FindOrNull(read_priors, prior.first);
    EXPECT_TRUE(camera_prior != nullptr);
    // Check basic entries.
    EXPECT_EQ(camera_prior->camera_intrinsics_model_type,
              prior.second.camera_intrinsics_model_type);
    EXPECT_EQ(camera_prior->image_width, prior.second.image_width);
    EXPECT_EQ(camera_prior->image_height, prior.second.image_height);
    EXPECT_EQ(camera_prior->focal_length.is_set,
              prior.second.focal_length.is_set);
    EXPECT_EQ(camera_prior->focal_length.value[0],
              prior.second.focal_length.value[0]);
    EXPECT_EQ(camera_prior->aspect_ratio.is_set,
              prior.second.aspect_ratio.is_set);
    EXPECT_EQ(camera_prior->aspect_ratio.value[0],
              prior.second.aspect_ratio.value[0]);
    EXPECT_EQ(camera_prior->skew.is_set, prior.second.skew.is_set);
    EXPECT_EQ(camera_prior->skew.value[0], prior.second.skew.value[0]);
    EXPECT_EQ(camera_prior->principal_point.is_set,
              prior.second.principal_point.is_set);
    EXPECT_EQ(camera_prior->principal_point.value[0],
              prior.second.principal_point.value[0]);
    EXPECT_EQ(camera_prior->principal_point.value[1],
              prior.second.principal_point.value[1]);
    EXPECT_EQ(camera_prior->radial_distortion.is_set,
              prior.second.radial_distortion.is_set);
    EXPECT_EQ(camera_prior->radial_distortion.value[0],
              prior.second.radial_distortion.value[0]);
    EXPECT_EQ(camera_prior->radial_distortion.value[1],
              prior.second.radial_distortion.value[1]);
    EXPECT_EQ(camera_prior->radial_distortion.value[2],
              prior.second.radial_distortion.value[2]);
    EXPECT_EQ(camera_prior->radial_distortion.value[3],
              prior.second.radial_distortion.value[3]);
    EXPECT_EQ(camera_prior->tangential_distortion.is_set,
              prior.second.tangential_distortion.is_set);
    EXPECT_EQ(camera_prior->tangential_distortion.value[0],
              prior.second.tangential_distortion.value[0]);
    EXPECT_EQ(camera_prior->tangential_distortion.value[1],
              prior.second.tangential_distortion.value[1]);
    EXPECT_EQ(camera_prior->position.is_set,
              prior.second.position.is_set);
    EXPECT_EQ(camera_prior->position.value[0],
              prior.second.position.value[0]);
    EXPECT_EQ(camera_prior->position.value[1],
              prior.second.position.value[1]);
    EXPECT_EQ(camera_prior->position.value[2],
              prior.second.position.value[2]);
    EXPECT_EQ(camera_prior->orientation.is_set,
              prior.second.orientation.is_set);
    EXPECT_EQ(camera_prior->orientation.value[0],
              prior.second.orientation.value[0]);
    EXPECT_EQ(camera_prior->orientation.value[1],
              prior.second.orientation.value[1]);
    EXPECT_EQ(camera_prior->orientation.value[2],
              prior.second.orientation.value[2]);
    EXPECT_EQ(camera_prior->latitude.is_set,
              prior.second.latitude.is_set);
    EXPECT_EQ(camera_prior->latitude.value[0],
              prior.second.latitude.value[0]);
    EXPECT_EQ(camera_prior->longitude.is_set,
              prior.second.longitude.is_set);
    EXPECT_EQ(camera_prior->longitude.value[0],
              prior.second.longitude.value[0]);
    EXPECT_EQ(camera_prior->altitude.is_set,
              prior.second.altitude.is_set);
    EXPECT_EQ(camera_prior->altitude.value[0],
              prior.second.altitude.value[0]);
  }
}

}  // namespace
}  // namespace theia
