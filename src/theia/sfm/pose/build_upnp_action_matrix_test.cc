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
#include "theia/sfm/pose/build_upnp_action_matrix.h"

namespace theia {
namespace {

TEST(BuildUpnpActionMatrixTests, GaussJordanEliminationOnSquaredMatrix) {
  const int kNumRows = 32;
  RowMajorMatrixXd mat = RowMajorMatrixXd::Random(kNumRows, kNumRows);
  //  GaussJordanElimination(kNumRows - 1, &mat);
  GaussJordanElimination(0, &mat);
  // Trace of matrix must be equals to the number of rows.
  EXPECT_NEAR(mat.trace(), static_cast<double>(kNumRows), 1e-6);
  // Verify that the lower triangular part sums to the trace.
  EXPECT_NEAR(mat.sum(), mat.trace(), 1e-6);
}

TEST(BuildUpnpActionMatrixTests, GaussJordanEliminationOnFatMatrix) {
  const int kNumRows = 32;
  const int kNumCols = kNumRows + 4;
  RowMajorMatrixXd mat = RowMajorMatrixXd::Random(kNumRows, kNumCols);
  GaussJordanElimination(0, &mat);
  // Verify that the left-block (rows, rows) is diagonalized.
  EXPECT_NEAR(mat.block(0, 0, kNumRows, kNumRows).sum(),
              mat.block(0, 0, kNumRows, kNumRows).trace(), 1e-6);
}

TEST(BuildUpnpActionMatrixTests, PartialDiagonalizationOnFatMatrix) {
  const int kNumRows = 32;
  const int kNumCols = kNumRows + 4;
  const int kLastRowToProcess = 2;
  RowMajorMatrixXd mat = RowMajorMatrixXd::Random(kNumRows, kNumCols);
  GaussJordanElimination(kLastRowToProcess, &mat);
  // Verify that the left-block (rows, rows) is partially diagonalized.
  EXPECT_NEAR(mat.block(kLastRowToProcess, kLastRowToProcess,
                        kNumRows - kLastRowToProcess,
                        kNumRows - kLastRowToProcess).sum(),
              kNumRows - kLastRowToProcess, 1e-6);
  EXPECT_NEAR(mat.block(kLastRowToProcess, kLastRowToProcess,
                        kNumRows - kLastRowToProcess,
                        kNumRows - kLastRowToProcess).sum(),
              mat.block(kLastRowToProcess, kLastRowToProcess,
                        kNumRows - kLastRowToProcess,
                        kNumRows - kLastRowToProcess).trace(),
              1e-6);
  EXPECT_NE(mat.block(0, 0, kNumRows, kNumRows).sum(),
            mat.block(kLastRowToProcess, kLastRowToProcess,
                      kNumRows - kLastRowToProcess,
                      kNumRows - kLastRowToProcess).sum());
}

TEST(BuildUpnpActionMatrixTests, BuildActionMatrixForCentralCameraMinimalCase) {
  // The upnp cost parmaeters are from
  // UpnpTests.MinimalCentralCameraPoseEstimation.
  Eigen::Matrix<double, 10, 10> a_matrix;
  a_matrix <<
      00.646311, 0-1.80655, 001.32834, -0.168107, 001.41862, 00.730516,
      0-5.65344, -0.481733, 00.621246, 001.14998,
      0-1.80655, 0012.6466, 0-11.5022, 00.662098, 0-3.57447, 0-4.33037,
      0015.7804, 00.247601, 004.23048, 0-9.49427,
      001.32834, 0-11.5022, 0011.8782, 0-1.70435, -0.352207, 002.68112,
      0-11.4645, 00-2.4114, 0-5.02743, 006.33111,
      -0.168107, 00.662098, 0-1.70435, 001.21036, 002.50805, 00.918736,
      001.33756, 002.64553, 00.175711, 002.01318,
      001.41862, 0-3.57447, -0.352207, 002.50805, 000013.43, 007.37281,
      0-12.7369, 005.69076, 0-1.12633, 0010.7572,
      00.730516, 0-4.33037, 002.68112, 00.918736, 007.37281, 000010.18,
      0-6.51633, 0014.5573, 0-4.31214, 0012.1559,
      0-5.65344, 0015.7804, 0-11.4645, 001.33756, 0-12.7369, 0-6.51633,
      0049.4673, 003.92669, 0-5.47262, 0-10.3226,
      -0.481733, 00.247601, 00-2.4114, 002.64553, 005.69076, 0014.5573,
      003.92669, 0035.5648, 0-2.34604, 0016.3519,
      00.621246, 004.23048, 0-5.02743, 00.175711, 0-1.12633, 0-4.31214,
      0-5.47262, 0-2.34604, 008.63544, 0-6.23962,
      001.14998, 0-9.49427, 006.33111, 002.01318, 0010.7572, 0012.1559,
      0-10.3226, 0016.3519, 0-6.23962, 0018.1126;
  Eigen::Matrix<double, 10, 1> b_vector;
  b_vector <<
      -1.82424e-30, 
      -1.97215e-31,
      01.18329e-30,
      01.04771e-31,
      01.52842e-30,
      -2.71171e-30,
      01.47911e-30,
      -3.15544e-30,
      -1.77494e-30,
      -5.91646e-31;
  const double gamma = 3.41921e-30;
  // Expected action matrix;
  Eigen::Matrix<double, 16, 16> action_matrix;
  action_matrix <<
      0, -0.0590812, -0.755901, -0.615635, 1.63155, 0, 0, 0, 0, 0, 0, 0.0297073,
      1.17114, 0.813402,-0.274845, 0,
      -2.85509, 0, 0, 0, 0, 0.980003, -1.85333, 4.59084, -0.188048, 1.93117,
      4.71699, 0, 0, 0, 0, -0.574897,
      -1.59155, 0, 0, 0, 0, 0.893845, -0.709019, 2.50022, 1.06107, 0.945678,
      2.56959, 0, 0, 0, 0, -0.313883,
      -12.2561, 0, 0, 0, 0, 4.3892, -7.79232, 19.8807, -1.91047, 8.64316,
      20.0887, 0, 0, 0, 0, -2.49153,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0.790771, 0.614409, 0.297029, -0.294195, 0, 0, 0, 0, 0, 0, -0.0101929,
      -0.942889, -0.660286, 0.0042962, 0,
      0, 4.29232, 0.316754, 0.469403,-0.268756, 0, 0, 0, 0, 0, 0, -0.0715431,
      -1.29034, -0.926158, 0.148596, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
  const Eigen::Matrix<double, 16, 16> computed_action_matrix =
      BuildActionMatrix(a_matrix, b_vector, gamma);
  EXPECT_NEAR((computed_action_matrix - action_matrix).squaredNorm(),
              0.0, 1e-6);
}

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
