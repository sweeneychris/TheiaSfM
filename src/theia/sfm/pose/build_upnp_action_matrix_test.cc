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

}  // namespace
}  // namespace theia
