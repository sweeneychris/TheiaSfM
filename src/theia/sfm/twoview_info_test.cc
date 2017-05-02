// Copyright (C) 2014 The Regents of the University of California (Regents).
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

#include "gtest/gtest.h"

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/twoview_info.h"
#include "theia/util/random.h"

namespace theia {

namespace {

static const int kNumTrials = 100;

RandomNumberGenerator rng(45);

void CheckSwappedCamera(const TwoViewInfo& info,
                        const TwoViewInfo& swapped_info) {
  // Check that the intrinsics (focal lengths) were swapped.
  EXPECT_EQ(info.focal_length_1, swapped_info.focal_length_2);
  EXPECT_EQ(info.focal_length_2, swapped_info.focal_length_1);

  // Swapping the camera pair twice should result in the original setup.
  TwoViewInfo inverse_of_inverse_info = swapped_info;
  SwapCameras(&inverse_of_inverse_info);

  // Check rotations.
  static const double kPoseTolerance = 1e-10;
  EXPECT_LT((inverse_of_inverse_info.rotation_2 - info.rotation_2).norm(),
            kPoseTolerance);
  EXPECT_LT((inverse_of_inverse_info.position_2 - info.position_2).norm(),
            kPoseTolerance);
}

}  // namespace

TEST(TwoViewInfo, SwapCameras) {
  for (int i = 0; i < kNumTrials; i++) {
    TwoViewInfo info;
    const Eigen::Vector2d focal_lengths = rng.RandVector2d(500.0, 1000.0);
    info.focal_length_1 = focal_lengths[0];
    info.focal_length_2 = focal_lengths[1];
    info.rotation_2 = rng.RandVector3d();
    info.position_2 = rng.RandVector3d();

    TwoViewInfo info2 = info;
    SwapCameras(&info2);
    CheckSwappedCamera(info, info2);
  }
}

// Test constructing a TwoViewInfo when camera1 has an identity pose. The
// TwoViewInfo should exactly contain camera2's "global" pose.
TEST(TwoViewInfo, TwoViewInfoFromTwoCamerasIdentity) {
  static const double kFocalLength1 = 800.0;
  static const double kFocalLength2 = 1200.0;
  static const Eigen::Vector3d kOrientation1(0.0, 0.0, 0.0);
  static const Eigen::Vector3d kOrientation2(0.1, 0.7, -0.2);
  static const Eigen::Vector3d kPosition1(0.0, 0.0, 0.0);
  static const Eigen::Vector3d kPosition2(10.4, 11.3, 4.2);
  Camera camera1, camera2;
  camera1.SetFocalLength(kFocalLength1);
  camera1.SetOrientationFromAngleAxis(kOrientation1);
  camera1.SetPosition(kPosition1);
  camera2.SetFocalLength(kFocalLength2);
  camera2.SetOrientationFromAngleAxis(kOrientation2);
  camera2.SetPosition(kPosition2);

  TwoViewInfo info_camera1_reference;
  TwoViewInfoFromTwoCameras(camera1, camera2, &info_camera1_reference);
  EXPECT_EQ(info_camera1_reference.focal_length_1, kFocalLength1);
  EXPECT_EQ(info_camera1_reference.focal_length_2, kFocalLength2);
  const Eigen::Vector3d normalized_position2 = kPosition2.normalized();
  for (int i = 0; i < 3; i++) {
    EXPECT_DOUBLE_EQ(info_camera1_reference.rotation_2[i], kOrientation2[i]);
    EXPECT_DOUBLE_EQ(info_camera1_reference.position_2[i],
                     normalized_position2[i]);
  }
}

// Test constructing a TwoViewInfo when camera1 has an identity pose. The
// TwoViewInfo should exactly contain camera2's "global" pose. Unlike the
// previous test, we first construct the TwoViewInfo with camera 2 as the
// reference, then swap the cameras. This is a nice test for both the swap and
// TwoViewInfoFromTwoCameras function to ensure they work in harmony as
// expected.
TEST(TwoViewInfo, TwoViewInfoFromTwoCamerasIdentitySwap) {
  static const double kFocalLength1 = 800.0;
  static const double kFocalLength2 = 1200.0;
  static const Eigen::Vector3d kOrientation1(0.0, 0.0, 0.0);
  static const Eigen::Vector3d kOrientation2(0.1, 0.7, -0.2);
  static const Eigen::Vector3d kPosition1(0.0, 0.0, 0.0);
  static const Eigen::Vector3d kPosition2(10.4, 11.3, 4.2);
  Camera camera1, camera2;
  camera1.SetFocalLength(kFocalLength1);
  camera1.SetOrientationFromAngleAxis(kOrientation1);
  camera1.SetPosition(kPosition1);
  camera2.SetFocalLength(kFocalLength2);
  camera2.SetOrientationFromAngleAxis(kOrientation2);
  camera2.SetPosition(kPosition2);

  TwoViewInfo info_camera2_reference;
  TwoViewInfoFromTwoCameras(camera2, camera1, &info_camera2_reference);
  TwoViewInfo info_camera1_reference = info_camera2_reference;
  SwapCameras(&info_camera1_reference);

  EXPECT_EQ(info_camera1_reference.focal_length_1, kFocalLength1);
  EXPECT_EQ(info_camera1_reference.focal_length_2, kFocalLength2);
  const Eigen::Vector3d normalized_position2 = kPosition2.normalized();
  for (int i = 0; i < 3; i++) {
    EXPECT_DOUBLE_EQ(info_camera1_reference.rotation_2[i], kOrientation2[i]);
    EXPECT_DOUBLE_EQ(info_camera1_reference.position_2[i],
                     normalized_position2[i]);
  }
}

// Construct TwoViewInfos with each camera as the reference. Then check that the
// swap function works in both directions.
TEST(TwoViewInfo, TwoViewInfoFromTwoCameras) {
  static const double kFocalLength1 = 800.0;
  static const double kFocalLength2 = 1200.0;
  static const Eigen::Vector3d kOrientation1(0.3, 1.2, -0.5);
  static const Eigen::Vector3d kOrientation2(0.1, 0.7, -0.2);
  static const Eigen::Vector3d kPosition1(10.1, 12.2, 3.3);
  static const Eigen::Vector3d kPosition2(10.4, 11.3, 4.2);
  Camera camera1, camera2;
  camera1.SetFocalLength(kFocalLength1);
  camera1.SetOrientationFromAngleAxis(kOrientation1);
  camera1.SetPosition(kPosition1);
  camera2.SetFocalLength(kFocalLength2);
  camera2.SetOrientationFromAngleAxis(kOrientation2);
  camera2.SetPosition(kPosition2);

  TwoViewInfo info_camera1_reference, info_camera2_reference;
  TwoViewInfoFromTwoCameras(camera1, camera2, &info_camera1_reference);
  TwoViewInfoFromTwoCameras(camera2, camera1, &info_camera2_reference);
  CheckSwappedCamera(info_camera1_reference, info_camera2_reference);
  CheckSwappedCamera(info_camera2_reference, info_camera1_reference);
}

}  // namespace theia
