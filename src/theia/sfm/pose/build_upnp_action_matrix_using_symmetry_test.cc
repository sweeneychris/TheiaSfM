// Copyright (C) 2018 The Regents of the University of California (Regents).
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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include "gtest/gtest.h"

#include <vector>

#include "theia/math/util.h"
#include "theia/sfm/pose/build_upnp_action_matrix_using_symmetry.h"

namespace theia {
namespace {

TEST(BuildUpnpActionMatrixTests,
     BuildActionMatrixUsingSymmetryForNonCentralMinimalCase) {
  // The upnp cost parmaeters are from
  // UpnpTests.MinimalNonCentralCameraPoseEstimation.
  Eigen::Matrix<double, 10, 10> a_matrix;
  a_matrix <<
      006.71742, 001.85533, 0-2.17812, 0-6.39463, 00-2.5622, 0-2.51401,
      001.00919, 0-5.13897, 004.51181, -0.658492,
      001.85533, 008.78853, 0-10.4981, -0.145733, 003.58971, 0-3.30778,
      000015.55, 0-5.03747, 001.45259, 00.299855,
      0-2.17812, 0-10.4981, 00013.386, -0.709788, 0-6.60231, 003.38488,
      00-16.878, 004.98896, 0-2.39375, 0-1.37423,
      0-6.39463, -0.145733, -0.709788, 007.25015, 0005.5748, 002.43691,
      00.318745, 005.18747, 0-3.57065, 001.73287,
      00-2.5622, 003.58971, 0-6.60231, 0005.5748, 0014.3744, 002.36133,
      007.77514, 00-6.9803, 0-7.06722, 005.97198,
      0-2.51401, 0-3.30778, 003.38488, 002.43691, 002.36133, 004.78982,
      001.28329, 007.11912, 0-5.71837, 005.08662,
      001.00919, 000015.55, 00-16.878, 00.318745, 007.77514, 001.28329,
      0053.9373, 0-3.14292, 0-14.5032, 0011.3284,
      0-5.13897, 0-5.03747, 004.98896, 005.18747, 00-6.9803, 007.11912,
      0-3.14292, 0040.8286, 005.32323, 005.22061,
      004.51181, 001.45259, 0-2.39375, 0-3.57065, 0-7.06722, 0-5.71837,
      0-14.5032, 005.32323, 0017.2046, 0-7.62227,
      -0.658492, 00.299855, 0-1.37423, 001.73287, 005.97198, 005.08662,
      0011.3284, 005.22061, 0-7.62227, 008.45574; 
  Eigen::Matrix<double, 10, 1> b_vector;
  b_vector <<
      00-6.6629,
      0-3.57868,
      004.05766,
      006.18392,
      001.58341,
      002.30622,
      0-7.06697,
      005.36014,
      0-2.77699,
      -0.646325;
  const double gamma = 7.29313;

  // Expected action matrix;
  Eigen::Matrix<double, 8, 8> action_matrix;
  action_matrix <<
      -0.45047, 0.390688, -1.15408, -0.0828841, -0.019765, 0.0389208, 0.0736624,
      0.226646,
      2.50121, 0.00267812, 2.12052, 1.27115, -0.0042492, -0.0786345, 0.0675796,
      -0.260323,
      2.27437, -0.648776, 1.05155, 1.24298, 0.0243789, -0.0625738, -0.0237413,
      -0.485711,
      -1.80143, -0.212471, -0.548587, -0.993304, 0.0345941, 0.08877, -0.105789,
      -0.0881157,
      1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0;
  const Eigen::Matrix<double, 8, 8> computed_action_matrix =
      BuildActionMatrixUsingSymmetry(a_matrix, b_vector, gamma);
  EXPECT_NEAR((computed_action_matrix - action_matrix).squaredNorm(),
              0.0, 1e-6);
}

}  // namespace
}  // namespace theia
