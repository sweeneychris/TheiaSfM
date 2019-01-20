// Copyright (C) 2013 The Regents of the University of California (Regents).
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
// Modified by Victor Fragoso (victor.fragoso@mail.wvu.edu)

#include <Eigen/Dense>
#include "gtest/gtest.h"

#include "theia/math/matrix/gauss_jordan.h"

namespace theia {
using RowMajorMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

TEST(GaussJordan, FullDiagonalizationOnSquaredRowMajorMatrix) {
  const int kNumRows = 32;
  RowMajorMatrixXd mat = RowMajorMatrixXd::Random(kNumRows, kNumRows);
  GaussJordan(&mat);
  // Trace of matrix must be equals to the number of rows.
  EXPECT_NEAR(mat.trace(), static_cast<double>(kNumRows), 1e-6);
  // Verify that the lower triangular part sums to the trace.
  EXPECT_NEAR(mat.sum(), mat.trace(), 1e-6);
}

TEST(GaussJordan, EliminationOnFatMatrix) {
  const int kNumRows = 32;
  const int kNumCols = kNumRows + 4;
  RowMajorMatrixXd mat = RowMajorMatrixXd::Random(kNumRows, kNumCols);
  GaussJordan(&mat);
  // Verify that the left-block (rows, rows) is diagonalized.
  EXPECT_NEAR(mat.block(0, 0, kNumRows, kNumRows).sum(),
              mat.block(0, 0, kNumRows, kNumRows).trace(), 1e-6);
}

}  // namespace theia
