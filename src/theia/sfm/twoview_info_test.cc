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

#include "theia/sfm/twoview_info.h"

namespace theia {

namespace {

static const int kNumTrials = 100;

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
    const Eigen::Vector2d focal_lengths =
        500.0 * (Eigen::Vector2d::Random() + Eigen::Vector2d::Ones());
    info.focal_length_1 = focal_lengths[0];
    info.focal_length_2 = focal_lengths[1];
    info.rotation_2 = Eigen::Vector3d::Random();
    info.position_2 = Eigen::Vector3d::Random();

    TwoViewInfo info2 = info;
    SwapCameras(&info2);
    CheckSwappedCamera(info, info2);
  }
}

}  // namespace theia
