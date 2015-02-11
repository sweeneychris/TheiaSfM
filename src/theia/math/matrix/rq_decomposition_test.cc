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
#include <Eigen/LU>

#include <algorithm>
#include <cmath>

#include "gtest/gtest.h"
#include "theia/math/matrix/rq_decomposition.h"

namespace theia {

using Eigen::Matrix;
using Eigen::Matrix3d;
using Eigen::Matrix3f;
using Eigen::Matrix4d;
using Eigen::Matrix4f;
using Eigen::MatrixXd;

const double kSmallTolerance = 1e-12;
const double kLargeTolerance = 1e-6;

template <typename MatrixType>
void TestRQDecomposition(const int rows, const int cols, const double tol) {
  for (int i = 0; i < 100; ++i) {
    const MatrixType& m = MatrixType::Random(rows, cols);
    RQDecomposition<MatrixType> rq(m);
    EXPECT_TRUE(rq.matrixR().bottomRows(std::min(rows, cols))
                    .isUpperTriangular());
    EXPECT_TRUE(rq.matrixQ().isUnitary());
    EXPECT_NEAR(rq.matrixQ().determinant(), 1.0, kLargeTolerance);

    const MatrixType& diff = rq.matrixR() * rq.matrixQ() - m;
    EXPECT_LE(diff.cwiseAbs().maxCoeff(), tol)
        << rq.matrixR() << "\n\n" << rq.matrixQ() << "\n\n"
        << rq.matrixR() * rq.matrixQ() << "\n\n" << m;
  }
}

TEST(RQTest, FixedSizeMatrices) {
  // Square matrices.
  TestRQDecomposition<Matrix3d>(3, 3, kSmallTolerance);
  TestRQDecomposition<Matrix4d>(4, 4, kSmallTolerance);
  TestRQDecomposition<Matrix<double, 5, 5> >(5, 5, kSmallTolerance);
  TestRQDecomposition<Matrix<double, 6, 6> >(6, 6, kSmallTolerance);

  // Non-square matrices.
  TestRQDecomposition<Matrix<double, 3, 4> >(3, 4, kSmallTolerance);
  TestRQDecomposition<Matrix<double, 4, 3> >(4, 3, kSmallTolerance);
}

TEST(RQTest, DynamicSizedMatrices) {
  // Square matrices.
  TestRQDecomposition<MatrixXd>(25, 25, kSmallTolerance);
  TestRQDecomposition<MatrixXd>(50, 50, kSmallTolerance);
  TestRQDecomposition<MatrixXd>(75, 75, kSmallTolerance);
  TestRQDecomposition<MatrixXd>(100, 100, kSmallTolerance);

  // Non-square matrices.
  TestRQDecomposition<MatrixXd>(30, 40, kSmallTolerance);
  TestRQDecomposition<MatrixXd>(40, 30, kSmallTolerance);
}

TEST(RQTest, FloatMatrices) {
  TestRQDecomposition<Matrix3f>(3, 3, kLargeTolerance);
  TestRQDecomposition<Matrix4f>(4, 4, kLargeTolerance);
}

TEST(RQTest, RowMajorMatrices) {
  TestRQDecomposition<Matrix<double, 3, 3, Eigen::RowMajor> >(3, 3,
                                                              kSmallTolerance);
  TestRQDecomposition<Matrix<float, 3, 3, Eigen::RowMajor> >(3, 3,
                                                             kLargeTolerance);
}

}  // namespace theia
