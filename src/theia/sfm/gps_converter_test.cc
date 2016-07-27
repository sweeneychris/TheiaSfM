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

#include <Eigen/Core>
#include <glog/logging.h>
#include "gtest/gtest.h"

#include "theia/sfm/gps_converter.h"
#include "theia/util/random.h"

namespace theia {

TEST(GPSConverter, ECEFToLLA) {
  static const double kTolerance = 1e-8;
  const Eigen::Vector3d taj_mahal_lla(27.173891, 78.042068, 168.0);
  const Eigen::Vector3d ecef(1176498.769459714,
                             5555043.905503586,
                             2895446.8901510699);
  const Eigen::Vector3d lla = GPSConverter::ECEFToLLA(ecef);
  EXPECT_NEAR(taj_mahal_lla[0], lla[0], kTolerance);
  EXPECT_NEAR(taj_mahal_lla[1], lla[1], kTolerance);
  EXPECT_NEAR(taj_mahal_lla[2], lla[2], kTolerance);
}

TEST(GPSConverter, LLAToECEF) {
  static const double kTolerance = 1e-8;
  const Eigen::Vector3d taj_mahal_lla(27.173891, 78.042068, 168.0);
  const Eigen::Vector3d gt_ecef(1176498.769459714,
                                5555043.905503586,
                                2895446.8901510699);
  const Eigen::Vector3d ecef = GPSConverter::LLAToECEF(taj_mahal_lla);
  EXPECT_NEAR(gt_ecef[0], ecef[0], kTolerance);
  EXPECT_NEAR(gt_ecef[1], ecef[1], kTolerance);
  EXPECT_NEAR(gt_ecef[2], ecef[2], kTolerance);
}

TEST(GPSConverter, RoundTrip) {
  InitRandomGenerator();
  static const double kTolerance = 1e-8;
  static const int kNumTrials = 1000.0;
  for (int i = 0; i < kNumTrials; i++) {
    // Use the same random configuration as in the original paper: Olson,
    // D.K. "Converting earth-Centered, Earth-Fixed Coordinates to Geodetic
    // Coordinates," IEEE Transactions on Aerospace and Electronic Systems,
    // Vol. 32, No. 1, January 1996, pp. 473-476.
    const Eigen::Vector3d gt_lla(RandDouble(-90.0, 90.0),
                                 RandDouble(-180.0, 180.0),
                                 RandDouble(-10000, 100000));

    const Eigen::Vector3d ecef = GPSConverter::LLAToECEF(gt_lla);
    const Eigen::Vector3d lla = GPSConverter::ECEFToLLA(ecef);

    EXPECT_NEAR(gt_lla[0], lla[0], kTolerance);
    EXPECT_NEAR(gt_lla[1], lla[1], kTolerance);
    EXPECT_NEAR(gt_lla[2], lla[2], kTolerance);
  }
}

}  // namespace theia
