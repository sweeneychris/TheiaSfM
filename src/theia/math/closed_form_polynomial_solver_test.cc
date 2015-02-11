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

#include <algorithm>
#include <complex>
#include "gtest/gtest.h"

#include "theia/math/closed_form_polynomial_solver.h"

namespace theia {

const double kTolerance = 1e-12;

TEST(SolveQuadraticPolynomial, DegenerateSolution) {
  // - 2x + 1 = 0
  const double a = 0.0;
  const double b = -2.0;
  const double c = 1.0;
  double roots[2];
  int num_roots = SolveQuadraticReals(a, b, c, kTolerance, roots);
  EXPECT_EQ(num_roots, 1);
  EXPECT_DOUBLE_EQ(roots[0], 0.5);
}

TEST(SolveQuadraticPolynomial, SolveReals) {
  // x^2 - 2x + 1 = 0
  double a = 1.0;
  double b = -2.0;
  double c = 1.0;
  double roots[2];
  int num_roots = SolveQuadraticReals(a, b, c, kTolerance, roots);
  EXPECT_EQ(num_roots, 2);
  EXPECT_DOUBLE_EQ(roots[0], 1.0);
  EXPECT_DOUBLE_EQ(roots[1], 1.0);

  // x^2 - 11x + 30 = 0
  a = 1.0;
  b = -11.0;
  c = 30.0;
  num_roots = SolveQuadraticReals(a, b, c, kTolerance, roots);
  EXPECT_EQ(num_roots, 2);

  const double soln_roots[2] = { 5.0, 6.0 };
  std::sort(std::begin(roots), std::end(roots));
  for (int i = 0; i < num_roots; i++) {
    EXPECT_DOUBLE_EQ(roots[i], soln_roots[i]);
  }
}

TEST(SolveQuadraticPolynomial, SolveComplex) {
  // x^2 - 2x + 5 = 0 should yield 1 + 2i, 1 - 2i
  const double a = 1.0;
  const double b = -2.0;
  const double c = 5.0;
  std::complex<double> roots[2];
  int num_roots = SolveQuadratic(a, b, c, roots);
  EXPECT_EQ(num_roots, 2);
  EXPECT_DOUBLE_EQ(roots[0].real(), 1.0);
  EXPECT_DOUBLE_EQ(roots[1].real(), 1.0);
  EXPECT_DOUBLE_EQ(fabs(roots[0].imag()), 2.0);
  EXPECT_DOUBLE_EQ(fabs(roots[1].imag()), 2.0);
}

TEST(SolveCubicPolynomial, SolveReals) {
  // x^3 - 6x^2 + 11x - 6 = 0
  const double a = 1.0;
  const double b = -6.0;
  const double c = 11.0;
  const double d = -6.0;
  double roots[3];
  int num_roots = SolveCubicReals(a, b, c, d, kTolerance, roots);
  EXPECT_EQ(num_roots, 3);

  // Check that each root is valid.
  std::sort(std::begin(roots), std::end(roots));
  const double soln_roots[3] = { 1.0, 2.0, 3.0 };
  for (int i = 0; i < num_roots; i++) {
    EXPECT_DOUBLE_EQ(roots[i], soln_roots[i]);
  }
}

TEST(SolveQuarticPolynomial, SolveReals) {
  // y = 3x^4 + 6x^3 - 123x^2 - 126x + 1080 = 0
  const long double a = 3.0;
  const long double b = 6.0;
  const long double c = -123.0;
  const long double d = -126.0;
  const long double e = 1080.0;

  long double roots[4];
  int num_roots = SolveQuarticReals(a, b, c, d, e, kTolerance, roots);
  EXPECT_EQ(num_roots, 4);

  // Check that each root is valid.
  std::sort(std::begin(roots), std::end(roots));
  const long double soln_roots[4] = { -6.0, -4.0, 3.0, 5.0 };
  for (int i = 0; i < num_roots; i++) {
    EXPECT_DOUBLE_EQ(roots[i], soln_roots[i]);
  }
}

}  // namespace theia
