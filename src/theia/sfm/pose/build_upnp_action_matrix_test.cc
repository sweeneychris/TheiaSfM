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

TEST(BuildUpnpActionMatrixTests, BuildActionMatrixForCentralMinimalCase) {
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

}  // namespace
}  // namespace theia
