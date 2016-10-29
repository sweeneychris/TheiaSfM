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

#ifndef THEIA_MATH_POLYNOMIAL_H_
#define THEIA_MATH_POLYNOMIAL_H_

#include <Eigen/Core>

namespace theia {

// All polynomials are assumed to be the form
//
//   sum_{i=0}^N polynomial(i) x^{N-i}.
//
// and are given by a vector of coefficients of size N + 1.

// Finds all roots of a polynomial. The real or imaginary output variable may be
// NULL if the respective output is not desired.
bool FindPolynomialRoots(const Eigen::VectorXd& polynomial,
                         Eigen::VectorXd* real,
                         Eigen::VectorXd* imaginary);

// Remove leading terms with zero coefficients.
Eigen::VectorXd RemoveLeadingZeros(const Eigen::VectorXd& polynomial_in);

// Evaluate the polynomial at x using the Horner scheme. This method is
// templated so that the polynomial may be evaluated at complex points.
template <typename T>
inline T EvaluatePolynomial(const Eigen::VectorXd& polynomial, const T& x) {
  T v = 0.0;
  for (int i = 0; i < polynomial.size(); ++i) {
    v = v * x + polynomial(i);
  }
  return v;
}

// Return the derivative of the given polynomial. It is assumed that
// the input polynomial is at least of degree zero.
Eigen::VectorXd DifferentiatePolynomial(const Eigen::VectorXd& polynomial);

// Multiplies the two polynoimals together.
Eigen::VectorXd MultiplyPolynomials(const Eigen::VectorXd& poly1,
                                    const Eigen::VectorXd& poly2);

// Adds two polynomials (of potentially different sizes) together.
Eigen::VectorXd AddPolynomials(const Eigen::VectorXd& poly1,
                               const Eigen::VectorXd& poly2);

// Performs polynomial division such that
// polynomial = divisor * quotient + remainder.
void DividePolynomial(const Eigen::VectorXd& polynomial,
                      const Eigen::VectorXd& divisor,
                      Eigen::VectorXd* quotient,
                      Eigen::VectorXd* remainder);

// Find the minimum value of the polynomial in the interval [x_min,
// x_max]. The minimum is obtained by computing all the roots of the
// derivative of the input polynomial. All real roots within the
// interval [x_min, x_max] are considered as well as the end points
// x_min and x_max. Since polynomials are differentiable functions,
// this ensures that the true minimum is found.
void MinimizePolynomial(const Eigen::VectorXd& polynomial, double x_min,
                        double x_max, double* optimal_x, double* optimal_value);

// Find roots of polynomials of the form a * x + b = 0. The real or imaginary
// output variable may be NULL if the respective output is not desired.
void FindLinearPolynomialRoots(const Eigen::VectorXd& polynomial,
                               Eigen::VectorXd* real,
                               Eigen::VectorXd* imaginary);

// Finds roots of polynomials of the form a * x^2 + b * x + c = 0. The real or
// imaginary output variable may be NULL if the respective output is not
// desired.
void FindQuadraticPolynomialRoots(const Eigen::VectorXd& polynomial,
                                  Eigen::VectorXd* real,
                                  Eigen::VectorXd* imaginary);

// Finds a single polynomial root iteratively with the starting position x0 and
// guaranteed precision of epsilon using Laguerre's method.
double FindRootIterativeLaguerre(const Eigen::VectorXd& polynomial,
                                 const double x0,
                                 const double epsilon,
                                 const int max_iterations);

// Find a root from the starting guess using Newton's method.
double FindRootIterativeNewton(const Eigen::VectorXd& polynomial,
                               const double x0,
                               const double epsilon,
                               const int max_iterations);
}  // namespace theia

#endif  // THEIA_MATH_POLYNOMIAL_H_
