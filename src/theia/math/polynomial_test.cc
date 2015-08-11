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

// These functions are largely based off the the Ceres solver polynomial
// functions which are not available through the public interface. The license
// is below:
//
// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2012 Google Inc. All rights reserved.
// http://code.google.com/p/ceres-solver/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: moll.markus@arcor.de (Markus Moll)
//         sameeragarwal@google.com (Sameer Agarwal)

#include <glog/logging.h>
#include <algorithm>
#include <vector>
#include "gtest/gtest.h"

#include "theia/math/polynomial.h"
#include "theia/test/test_utils.h"

namespace theia {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// For IEEE-754 doubles, machine precision is about 2e-16.
const double kEpsilon = 1e-13;

TEST(Polynomial, FindRootIterativeLaguerreTest) {
  VectorXd polynomial(6);
  // (x - 3) * (x + 4) * (x + 5) * (x - 6) * (x + 7)
  polynomial(0) = 1;
  polynomial(1) = 7;
  polynomial(2) = -43;
  polynomial(3) = -319;
  polynomial(4) = 234;
  polynomial(5) = 2520;

  const double kTolerance = 1e-10;
  const double kEpsilon = 1e-8;
  const int kMaxIter = 10;
  EXPECT_NEAR(FindRootIterativeLaguerre(polynomial, 3.1, kEpsilon, kMaxIter),
              3.0, kTolerance);
  EXPECT_NEAR(FindRootIterativeLaguerre(polynomial, -4.1, kEpsilon, kMaxIter),
              -4.0, kTolerance);
  EXPECT_NEAR(FindRootIterativeLaguerre(polynomial, -5.1, kEpsilon, kMaxIter),
              -5.0, kTolerance);
  EXPECT_NEAR(FindRootIterativeLaguerre(polynomial, 6.1, kEpsilon, kMaxIter),
              6.0, kTolerance);
  EXPECT_NEAR(FindRootIterativeLaguerre(polynomial, -7.1, kEpsilon, kMaxIter),
              -7.0, kTolerance);
}

TEST(Polynomial, FindRootIterativeNewtonTest) {
  VectorXd polynomial(6);
  // (x - 3) * (x + 4) * (x + 5) * (x - 6) * (x + 7)
  polynomial(0) = 1;
  polynomial(1) = 7;
  polynomial(2) = -43;
  polynomial(3) = -319;
  polynomial(4) = 234;
  polynomial(5) = 2520;

  const double kTolerance = 1e-10;
  const double kEpsilon = 1e-8;
  const int kMaxIter = 20;
  EXPECT_NEAR(FindRootIterativeNewton(polynomial, 3.1, kEpsilon, kMaxIter), 3.0,
              kTolerance);
  EXPECT_NEAR(FindRootIterativeNewton(polynomial, -4.1, kEpsilon, kMaxIter),
              -4.0, kTolerance);
  EXPECT_NEAR(FindRootIterativeNewton(polynomial, -5.1, kEpsilon, kMaxIter),
              -5.0, kTolerance);
  EXPECT_NEAR(FindRootIterativeNewton(polynomial, 6.1, kEpsilon, kMaxIter), 6.0,
              kTolerance);
  EXPECT_NEAR(FindRootIterativeNewton(polynomial, -7.1, kEpsilon, kMaxIter),
              -7.0, kTolerance);
}

TEST(Polynomial, DifferentiateConstantPolynomial) {
  // p(x) = 1;
  VectorXd polynomial(1);
  polynomial(0) = 1.0;
  const VectorXd derivative = DifferentiatePolynomial(polynomial);
  EXPECT_EQ(derivative.rows(), 1);
  EXPECT_EQ(derivative(0), 0);
}

TEST(Polynomial, DifferentiateQuadraticPolynomial) {
  // p(x) = x^2 + 2x + 3;
  VectorXd polynomial(3);
  polynomial(0) = 1.0;
  polynomial(1) = 2.0;
  polynomial(2) = 3.0;

  const VectorXd derivative = DifferentiatePolynomial(polynomial);
  EXPECT_EQ(derivative.rows(), 2);
  EXPECT_EQ(derivative(0), 2.0);
  EXPECT_EQ(derivative(1), 2.0);
}

TEST(Polynomial, MultiplyPolynomials) {
  VectorXd poly1(3);
  poly1[0] = 2;
  poly1[1] = 1;
  poly1[2] = 1;

  VectorXd poly2 = VectorXd::Zero(3);
  poly2[0] = 1;

  const VectorXd multiplied_poly = MultiplyPolynomials(poly1, poly2);
  VectorXd expected_poly(5);
  expected_poly.setZero();
  expected_poly[0] = 2;
  expected_poly[1] = 1;
  expected_poly[2] = 1;

  const double kTolerance = 1e-12;
  test::ExpectMatricesNear(expected_poly, multiplied_poly, kTolerance);
}

TEST(Polynomial, DividePolynomialSameDegree) {
  VectorXd poly1(3);
  poly1[0] = 2;
  poly1[1] = 1;
  poly1[2] = 1;

  VectorXd poly2 = VectorXd::Zero(3);
  poly2[0] = 1;

  VectorXd quotient, remainder;
  DividePolynomial(poly1, poly2, &quotient, &remainder);
  VectorXd reconstructed_poly =
      MultiplyPolynomials(quotient, poly2);
  reconstructed_poly.tail(remainder.size()) += remainder;

  const double kTolerance = 1e-12;
  test::ExpectMatricesNear(poly1, reconstructed_poly, kTolerance);
}

TEST(Polynomial, DividePolynomialLowerDegree) {
  VectorXd poly1(3);
  poly1[0] = 2;
  poly1[1] = 1;
  poly1[2] = 1;

  VectorXd poly2 = VectorXd::Zero(2);
  poly2[0] = 1;

  VectorXd quotient, remainder;
  DividePolynomial(poly1, poly2, &quotient, &remainder);
  VectorXd reconstructed_poly =
      MultiplyPolynomials(quotient, poly2);
  reconstructed_poly.tail(remainder.size()) += remainder;

  const double kTolerance = 1e-12;
  test::ExpectMatricesNear(poly1, reconstructed_poly, kTolerance);
}

TEST(Polynomial, DividePolynomialHigherDegree) {
  VectorXd poly1(3);
  poly1[0] = 2;
  poly1[1] = 1;
  poly1[2] = 1;

  VectorXd poly2 = VectorXd::Zero(3);
  poly2[0] = 1;

  VectorXd quotient, remainder;
  DividePolynomial(poly1, poly2, &quotient, &remainder);
  VectorXd reconstructed_poly =
      MultiplyPolynomials(quotient, poly2);
  reconstructed_poly.tail(remainder.size()) += remainder;

  const double kTolerance = 1e-12;
  test::ExpectMatricesNear(poly1, reconstructed_poly, kTolerance);
}

TEST(Polynomial, DividePolynomial) {
  VectorXd poly1(3);
  poly1[0] = 2;
  poly1[1] = 1;
  poly1[2] = 1;

  VectorXd poly2 = VectorXd::Zero(3);
  poly2[0] = 1;

  VectorXd quotient, remainder;
  DividePolynomial(poly1, poly2, &quotient, &remainder);
  VectorXd reconstructed_poly =
      MultiplyPolynomials(quotient, poly2);
  reconstructed_poly.tail(remainder.size()) += remainder;

  const double kTolerance = 1e-12;
  test::ExpectMatricesNear(poly1, reconstructed_poly, kTolerance);
}

TEST(Polynomial, MinimizeConstantPolynomial) {
  // p(x) = 1;
  VectorXd polynomial(1);
  polynomial(0) = 1.0;

  double optimal_x = 0.0;
  double optimal_value = 0.0;
  double min_x = 0.0;
  double max_x = 1.0;
  MinimizePolynomial(polynomial, min_x, max_x, &optimal_x, &optimal_value);

  EXPECT_EQ(optimal_value, 1.0);
  EXPECT_LE(optimal_x, max_x);
  EXPECT_GE(optimal_x, min_x);
}

TEST(Polynomial, MinimizeLinearPolynomial) {
  // p(x) = x - 2
  VectorXd polynomial(2);

  polynomial(0) = 1.0;
  polynomial(1) = 2.0;

  double optimal_x = 0.0;
  double optimal_value = 0.0;
  double min_x = 0.0;
  double max_x = 1.0;
  MinimizePolynomial(polynomial, min_x, max_x, &optimal_x, &optimal_value);

  EXPECT_EQ(optimal_x, 0.0);
  EXPECT_EQ(optimal_value, 2.0);
}

TEST(Polynomial, MinimizeQuadraticPolynomial) {
  // p(x) = x^2 - 3 x + 2
  // min_x = 3/2
  // min_value = -1/4;
  VectorXd polynomial(3);
  polynomial(0) = 1.0;
  polynomial(1) = -3.0;
  polynomial(2) = 2.0;

  double optimal_x = 0.0;
  double optimal_value = 0.0;
  double min_x = -2.0;
  double max_x = 2.0;
  MinimizePolynomial(polynomial, min_x, max_x, &optimal_x, &optimal_value);
  EXPECT_EQ(optimal_x, 3.0 / 2.0);
  EXPECT_EQ(optimal_value, -1.0 / 4.0);

  min_x = -2.0;
  max_x = 1.0;
  MinimizePolynomial(polynomial, min_x, max_x, &optimal_x, &optimal_value);
  EXPECT_EQ(optimal_x, 1.0);
  EXPECT_EQ(optimal_value, 0.0);

  min_x = 2.0;
  max_x = 3.0;
  MinimizePolynomial(polynomial, min_x, max_x, &optimal_x, &optimal_value);
  EXPECT_EQ(optimal_x, 2.0);
  EXPECT_EQ(optimal_value, 0.0);
}

}  // namespace theia
