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

#ifndef THEIA_TEST_TEST_UTILS_H_
#define THEIA_TEST_TEST_UTILS_H_

#include <Eigen/Dense>
#include <glog/logging.h>
#include "gtest/gtest.h"

namespace theia {
namespace test {

// Assert that values of the two matrices are nearly the same.
template <typename Derived>
void ExpectMatricesNear(const Eigen::MatrixBase<Derived>& a,
                        const Eigen::MatrixBase<Derived>& b,
                        double tolerance) {
  ASSERT_EQ(a.rows(), b.rows());
  ASSERT_EQ(a.cols(), b.cols());
  for (int i = 0; i < a.rows(); i++)
    for (int j = 0; j < a.cols(); j++)
      ASSERT_NEAR(a(i, j), b(i, j), tolerance)
          << "Entry (" << i << ", " << j << ") did not meet the tolerance!";
}

void ExpectArraysNear(int n,
                      const double* a,
                      const double* b,
                      double tolerance) {
  ASSERT_GT(n, 0);
  CHECK(a);
  CHECK(b);
  for (int i = 0; i < n; i++) {
    EXPECT_NEAR(a[i], b[i], tolerance) << "i = " << i;
  }
}

// Expects that for all i = 1,.., n - 1
//
//   |p[i] / max_norm_p - q[i] / max_norm_q| < tolerance
//
// where max_norm_p and max_norm_q are the max norms of the arrays p
// and q respectively.
bool ArraysEqualUpToScale(int n, const double* p, const double* q,
                          double tolerance) {
  Eigen::Map<const Eigen::VectorXd> p_vec(p, n);
  Eigen::Map<const Eigen::VectorXd> q_vec(q, n);

  // Use the cos term in the dot product to determine equality normalized for
  // scale.
  const double cos_diff = p_vec.dot(q_vec) / (p_vec.norm() * q_vec.norm());
  return std::abs(cos_diff) >= 1.0 - tolerance;
}

}  // namespace test
}  // namespace theia
#endif  // THEIA_TEST_TEST_UTILS_H_
