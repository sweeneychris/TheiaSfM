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

#include <Eigen/Core>
#include <glog/logging.h>

#include "gtest/gtest.h"
#include "theia/sfm/camera/radial_distortion.h"

namespace theia {

using Eigen::Vector2d;

static const int kNumTrials = 1;
static const double kTolerance = 1e-8;

void TestRadialDistortion(const double k1, const double k2) {
  for (int i = 0; i < kNumTrials; i++) {
    const Vector2d undistorted_point = Vector2d::Random();
    Vector2d distorted_point;
    RadialDistortPoint(undistorted_point.x(),
                       undistorted_point.y(),
                       k1,
                       k2,
                       distorted_point.data(),
                       distorted_point.data() + 1);

    Vector2d new_undistorted_point;
    RadialUndistortPoint(distorted_point,
                         k1,
                         k2,
                         &new_undistorted_point);

    EXPECT_NEAR(undistorted_point.x(), new_undistorted_point.x(), kTolerance);
    EXPECT_NEAR(undistorted_point.y(), new_undistorted_point.y(), kTolerance);
  }
}

TEST(RadialUndistort, ZeroDistortion) {
  TestRadialDistortion(0.0, 0.0);
}

TEST(RadialUndistort, OneParameterDistortion) {
  TestRadialDistortion(0.1, 0.0);
}

TEST(RadialUndistort, TwoParameterDistortion) {
  TestRadialDistortion(0.1, 0.05);
}

}  // namespace theia
