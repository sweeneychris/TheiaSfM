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

#include "theia/math/polynomial.h"

#include <Eigen/Core>
#include <glog/logging.h>

#include <cmath>
#include <limits>

#include "theia/math/find_polynomial_roots_jenkins_traub.h"

namespace theia {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

bool FindPolynomialRoots(const VectorXd& polynomial,
                         VectorXd* real,
                         VectorXd* imaginary) {
  return FindPolynomialRootsJenkinsTraub(polynomial, real, imaginary);
}

// Remove leading terms with zero coefficients.
VectorXd RemoveLeadingZeros(const VectorXd& polynomial_in) {
  int i = 0;
  while (i < (polynomial_in.size() - 1) && polynomial_in(i) == 0) {
    ++i;
  }
  return polynomial_in.tail(polynomial_in.size() - i);
}

VectorXd DifferentiatePolynomial(const VectorXd& polynomial) {
  const int degree = polynomial.rows() - 1;
  CHECK_GE(degree, 0);

  // Degree zero polynomials are constants, and their derivative does
  // not result in a smaller degree polynomial, just a degree zero
  // polynomial with value zero.
  if (degree == 0) {
    return VectorXd::Zero(1);
  }

  VectorXd derivative(degree);
  for (int i = 0; i < degree; ++i) {
    derivative(i) = (degree - i) * polynomial(i);
  }

  return derivative;
}

VectorXd MultiplyPolynomials(const VectorXd& poly1, const VectorXd& poly2) {
  VectorXd multiplied_poly = VectorXd::Zero(poly1.size() + poly2.size() - 1);;
  for (int i = 0; i < poly1.size(); i++) {
    for (int j = 0; j < poly2.size(); j++) {
      multiplied_poly.reverse().operator()(i + j) +=
          poly1.reverse()(i) * poly2.reverse()(j);
    }
  }
  return multiplied_poly;
}

void DividePolynomial(const VectorXd& polynomial,
                      const VectorXd& divisor,
                      VectorXd* quotient,
                      VectorXd* remainder) {
  // If the divisor is higher degree than the polynomial then it cannot be
  // divided so we simply return the remainder.
  if (polynomial.size() < divisor.size()) {
    *quotient = VectorXd::Zero(1);
    *remainder = polynomial;
    return;
  }

  VectorXd numerator = RemoveLeadingZeros(polynomial);
  VectorXd denominator;
  *quotient = VectorXd::Zero(numerator.size() - divisor.size() + 1);
  while (numerator.size() >= divisor.size()) {
    denominator = VectorXd::Zero(numerator.size());
    denominator.head(divisor.size()) = divisor;

    const double quotient_scalar = numerator(0) / denominator(0);
    quotient->reverse().operator()(numerator.size() - divisor.size()) =
        quotient_scalar;
    denominator = denominator * quotient_scalar;
    numerator = numerator - denominator;
    // Sometimes there are floating point errors that result in a non-zero first
    // value.
    numerator(0) = 0;
    numerator = RemoveLeadingZeros(numerator);
  }
  *remainder = numerator;
}

VectorXd AddPolynomials(const VectorXd& poly1, const VectorXd& poly2) {
  if (poly1.size() > poly2.size()) {
    VectorXd sum = poly1;
    sum.tail(poly2.size()) += poly2;
    return sum;
  } else {
    VectorXd sum = poly2;
    sum.tail(poly1.size()) += poly1;
    return sum;
  }
}

void FindLinearPolynomialRoots(const VectorXd& polynomial,
                               VectorXd* real,
                               VectorXd* imaginary) {
  CHECK_EQ(polynomial.size(), 2);
  if (real != NULL) {
    real->resize(1);
    (*real)(0) = -polynomial(1) / polynomial(0);
  }

  if (imaginary != NULL) {
    imaginary->setZero(1);
  }
}

void FindQuadraticPolynomialRoots(const VectorXd& polynomial,
                                  VectorXd* real,
                                  VectorXd* imaginary) {
  CHECK_EQ(polynomial.size(), 3);
  const double a = polynomial(0);
  const double b = polynomial(1);
  const double c = polynomial(2);
  const double D = b * b - 4 * a * c;
  const double sqrt_D = sqrt(fabs(D));
  if (real != NULL) {
    real->setZero(2);
  }
  if (imaginary != NULL) {
    imaginary->setZero(2);
  }

  // Real roots.
  if (D >= 0) {
    if (real != NULL) {
      // Stable quadratic roots according to BKP Horn.
      // http://people.csail.mit.edu/bkph/articles/Quadratics.pdf
      if (b >= 0) {
        (*real)(0) = (-b - sqrt_D) / (2.0 * a);
        (*real)(1) = (2.0 * c) / (-b - sqrt_D);
      } else {
        (*real)(0) = (2.0 * c) / (-b + sqrt_D);
        (*real)(1) = (-b + sqrt_D) / (2.0 * a);
      }
    }
    return;
  }

  // Use the normal quadratic formula for the complex case.
  if (real != NULL) {
    (*real)(0) = -b / (2.0 * a);
    (*real)(1) = -b / (2.0 * a);
  }
  if (imaginary != NULL) {
    (*imaginary)(0) = sqrt_D / (2.0 * a);
    (*imaginary)(1) = -sqrt_D / (2.0 * a);
  }
}

void MinimizePolynomial(const VectorXd& polynomial,
                        const double x_min,
                        const double x_max,
                        double* optimal_x,
                        double* optimal_value) {
  // Find the minimum of the polynomial at the two ends.
  //
  // We start by inspecting the middle of the interval. Technically
  // this is not needed, but we do this to make this code as close to
  // the minFunc package as possible.
  *optimal_x = (x_min + x_max) / 2.0;
  *optimal_value = EvaluatePolynomial(polynomial, *optimal_x);

  const double x_min_value = EvaluatePolynomial(polynomial, x_min);
  if (x_min_value < *optimal_value) {
    *optimal_value = x_min_value;
    *optimal_x = x_min;
  }

  const double x_max_value = EvaluatePolynomial(polynomial, x_max);
  if (x_max_value < *optimal_value) {
    *optimal_value = x_max_value;
    *optimal_x = x_max;
  }

  // If the polynomial is linear or constant, we are done.
  if (polynomial.rows() <= 2) {
    return;
  }

  const VectorXd derivative = DifferentiatePolynomial(polynomial);
  VectorXd roots_real;
  if (!FindPolynomialRoots(derivative, &roots_real, NULL)) {
    LOG(WARNING) << "Unable to find the critical points of "
                 << "the interpolating polynomial.";
    return;
  }

  // This is a bit of an overkill, as some of the roots may actually
  // have a complex part, but its simpler to just check these values.
  for (int i = 0; i < roots_real.rows(); ++i) {
    const double root = roots_real(i);
    if ((root < x_min) || (root > x_max)) {
      continue;
    }

    const double value = EvaluatePolynomial(polynomial, root);
    if (value < *optimal_value) {
      *optimal_value = value;
      *optimal_x = root;
    }
  }
}

// An iterative solver to find the closest root based on an initial guess. We
// use Laguerre's method, which is a polynomial root finding method that
// converges to a root with very high certainty. For multiple roots, the
// convergence is linear, otherwise it is cubic.
double FindRootIterativeLaguerre(const VectorXd& polynomial,
                                 const double x0,
                                 const double epsilon,
                                 const int max_iter) {
  const double kSmallestValue = 1e-10;

  // Constant symbolic derivitives.
  const VectorXd f_prime = DifferentiatePolynomial(polynomial);
  const VectorXd f_prime_prime = DifferentiatePolynomial(f_prime);
  const double k = static_cast<double>(polynomial.size());

  double x = x0;

  for (int i = 0; i < max_iter; i++) {
    const double f_of_x = EvaluatePolynomial(polynomial, x);
    if (std::abs(f_of_x) < kSmallestValue) {
      break;
    }

    const double g = EvaluatePolynomial(f_prime, x) / f_of_x;
    const double h = g * g - EvaluatePolynomial(f_prime_prime, x) / f_of_x;
    const double denom_part = std::sqrt(std::abs((k - 1.0) * (k * h - g * g)));
    const double denom = (g < 0) ? g - denom_part : g + denom_part;
    const double delta =  k / denom;
    if (std::abs(delta) < epsilon) {
      break;
    }

    x -= delta;
  }

  return x;
}

double FindRootIterativeNewton(const Eigen::VectorXd& polynomial,
                               const double x0,
                               const double epsilon,
                               const int max_iterations) {
  double root = x0;
  const Eigen::VectorXd derivative = DifferentiatePolynomial(polynomial);
  double prev = std::numeric_limits<double>::max();
  for (int i = 0; i < max_iterations && std::abs(prev - root) > epsilon; i++) {
    prev = root;
    root -= EvaluatePolynomial(polynomial, root) /
            EvaluatePolynomial(derivative, root);
  }
  return root;
}

}  // namespace theia
