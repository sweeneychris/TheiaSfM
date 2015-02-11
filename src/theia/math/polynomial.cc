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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <glog/logging.h>

#include <cmath>
#include <complex>
#include <queue>
#include <vector>

#include "theia/math/sturm_chain.h"
#include "theia/util/random.h"

namespace theia {

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

// Holds an interval for sturm chains. The difference between the lower and
// upper sturm number is the number of roots contained in that interval.
struct SturmBound {
  double lower_bound;
  int lower_sturm_number;
  double upper_bound;
  int upper_sturm_number;
};

// Balancing function as described by B. N. Parlett and C. Reinsch,
// "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors".
// In: Numerische Mathematik, Volume 13, Number 4 (1969), 293-304,
// Springer Berlin / Heidelberg. DOI: 10.1007/BF02165404
void BalanceCompanionMatrix(MatrixXd* companion_matrix_ptr) {
  CHECK_NOTNULL(companion_matrix_ptr);
  MatrixXd& companion_matrix = *companion_matrix_ptr;
  MatrixXd companion_matrix_offdiagonal = companion_matrix;
  companion_matrix_offdiagonal.diagonal().setZero();

  const int degree = companion_matrix.rows();

  // gamma <= 1 controls how much a change in the scaling has to
  // lower the 1-norm of the companion matrix to be accepted.
  //
  // gamma = 1 seems to lead to cycles (numerical issues?), so
  // we set it slightly lower.
  const double gamma = 0.9;

  // Greedily scale row/column pairs until there is no change.
  bool scaling_has_changed;
  do {
    scaling_has_changed = false;

    for (int i = 0; i < degree; ++i) {
      const double row_norm = companion_matrix_offdiagonal.row(i).lpNorm<1>();
      const double col_norm = companion_matrix_offdiagonal.col(i).lpNorm<1>();

      // Decompose row_norm/col_norm into mantissa * 2^exponent,
      // where 0.5 <= mantissa < 1. Discard mantissa (return value
      // of frexp), as only the exponent is needed.
      int exponent = 0;
      std::frexp(row_norm / col_norm, &exponent);
      exponent /= 2;

      if (exponent != 0) {
        const double scaled_col_norm = std::ldexp(col_norm, exponent);
        const double scaled_row_norm = std::ldexp(row_norm, -exponent);
        if (scaled_col_norm + scaled_row_norm < gamma * (col_norm + row_norm)) {
          // Accept the new scaling. (Multiplication by powers of 2 should not
          // introduce rounding errors (ignoring non-normalized numbers and
          // over- or underflow))
          scaling_has_changed = true;
          companion_matrix_offdiagonal.row(i) *= std::ldexp(1.0, -exponent);
          companion_matrix_offdiagonal.col(i) *= std::ldexp(1.0, exponent);
        }
      }
    }
  } while (scaling_has_changed);

  companion_matrix_offdiagonal.diagonal() = companion_matrix.diagonal();
  companion_matrix = companion_matrix_offdiagonal;
}

void BuildCompanionMatrix(const VectorXd& polynomial,
                          MatrixXd* companion_matrix_ptr) {
  CHECK_NOTNULL(companion_matrix_ptr);
  MatrixXd& companion_matrix = *companion_matrix_ptr;

  const int degree = polynomial.size() - 1;

  companion_matrix.resize(degree, degree);
  companion_matrix.setZero();
  companion_matrix.diagonal(-1).setOnes();
  companion_matrix.col(degree - 1) = -polynomial.reverse().head(degree);
}

// Remove leading terms with zero coefficients.
VectorXd RemoveLeadingZeros(const VectorXd& polynomial_in) {
  int i = 0;
  while (i < (polynomial_in.size() - 1) && polynomial_in(i) == 0) {
    ++i;
  }
  return polynomial_in.tail(polynomial_in.size() - i);
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

}  // namespace

bool FindPolynomialRoots(const VectorXd& polynomial_in,
                         VectorXd* real, VectorXd* imaginary) {
  if (polynomial_in.size() == 0) {
    LOG(ERROR) << "Invalid polynomial of size 0 passed to FindPolynomialRoots";
    return false;
  }

  VectorXd polynomial = RemoveLeadingZeros(polynomial_in);
  const int degree = polynomial.size() - 1;

  // Is the polynomial constant?
  if (degree == 0) {
    LOG(WARNING) << "Trying to extract roots from a constant "
                 << "polynomial in FindPolynomialRoots";
    // We return true with no roots, not false, as if the polynomial is constant
    // it is correct that there are no roots. It is not the case that they were
    // there, but that we have failed to extract them.
    return true;
  }

  // Linear
  if (degree == 1) {
    FindLinearPolynomialRoots(polynomial, real, imaginary);
    return true;
  }

  // Quadratic
  if (degree == 2) {
    FindQuadraticPolynomialRoots(polynomial, real, imaginary);
    return true;
  }

  // The degree is now known to be at least 3. For cubic or higher
  // roots we use the method of companion matrices.

  // Divide by leading term
  const double leading_term = polynomial(0);
  polynomial /= leading_term;

  // Build and balance the companion matrix to the polynomial.
  MatrixXd companion_matrix(degree, degree);
  BuildCompanionMatrix(polynomial, &companion_matrix);
  BalanceCompanionMatrix(&companion_matrix);

  // Find its (complex) eigenvalues.
  Eigen::EigenSolver<MatrixXd> solver(companion_matrix, false);
  if (solver.info() != Eigen::Success) {
    LOG(ERROR) << "Failed to extract eigenvalues from companion matrix.";
    return false;
  }

  // Output roots
  if (real != NULL) {
    *real = solver.eigenvalues().real();
  } else {
    LOG(WARNING) << "NULL pointer passed as real argument to "
                 << "FindPolynomialRoots. Real parts of the roots will not "
                 << "be returned.";
  }
  if (imaginary != NULL) {
    *imaginary = solver.eigenvalues().imag();
  }
  return true;
}

bool FindRealPolynomialRoots(const VectorXd& coeffs,
                             VectorXd* real_roots) {
  VectorXd real_roots_temp;
  VectorXd complex_roots;
  if (FindPolynomialRoots(coeffs, &real_roots_temp, &complex_roots)) {
    int num_real_roots = 0;
    real_roots->resize(complex_roots.size());
    for (int i = 0; i < complex_roots.size(); i++) {
      if (complex_roots(i) == 0.0) {
        (*real_roots)(num_real_roots) = real_roots_temp(i);
        num_real_roots++;
      }
    }
    real_roots->conservativeResize(num_real_roots);
    return true;
  }
  return false;
}

bool FindRealPolynomialRootsSturm(const VectorXd& coeffs,
                                  VectorXd* real_roots) {
  VectorXd polynomial = RemoveLeadingZeros(coeffs);
  const int degree = polynomial.size() - 1;
  if (degree <= 4) {
    return FindRealPolynomialRoots(polynomial, real_roots);
  }

  InitRandomGenerator();

  // Create the sturm chain.
  polynomial = polynomial / polynomial(0);
  SturmChain sturm_chain(polynomial);

  // Find finite bounds.
  SturmBound initial_bound;
  sturm_chain.ComputeRootBounds(&initial_bound.lower_bound,
                                &initial_bound.upper_bound);
  initial_bound.lower_sturm_number =
      sturm_chain.NumSignChanges(initial_bound.lower_bound);
  initial_bound.upper_sturm_number =
      sturm_chain.NumSignChanges(initial_bound.upper_bound);

  // Explore all intervals using a simple bisection method. Refine the intervals
  // until there is exactly one root within the interval, then use a root
  // finding method to extract that roo.
  static const double kEpsilon = 1e-8;
  static const double kMaxIter = 10;
  std::queue<SturmBound> bounds;
  bounds.emplace(initial_bound);
  std::vector<double> computed_roots;
  while (!bounds.empty()) {
    const SturmBound& current_interval = bounds.front();
    const double midpoint =
        (current_interval.upper_bound + current_interval.lower_bound) / 2.0;
    const int num_roots = current_interval.lower_sturm_number -
                          current_interval.upper_sturm_number;

    if (num_roots == 1 || current_interval.upper_bound - midpoint < kEpsilon) {
      // Try to find the root with a starting point at the middle of the
      // interval.
      double root = FindRealRootIterative(coeffs, midpoint, kEpsilon, kMaxIter);

      // If the root ends up not being within the bounds then retry with a
      // random initialization.
      while (root < current_interval.lower_bound ||
             root > current_interval.upper_bound) {
        root = FindRealRootIterative(coeffs,
                                     RandDouble(current_interval.lower_bound,
                                                current_interval.upper_bound),
                                     kEpsilon, kMaxIter);
      }
      computed_roots.emplace_back(root);
    } else if (num_roots > 1) {
      const int midpoint_sturm_number = sturm_chain.NumSignChanges(midpoint);

      // Split the interval in half and add the new roots.
      SturmBound new_lower_interval = current_interval;
      new_lower_interval.upper_bound = midpoint;
      new_lower_interval.upper_sturm_number = midpoint_sturm_number;
      bounds.emplace(new_lower_interval);

      SturmBound new_upper_interval = current_interval;
      new_upper_interval.lower_bound = midpoint;
      new_upper_interval.lower_sturm_number = midpoint_sturm_number;
      bounds.emplace(new_upper_interval);
    }
    bounds.pop();
  }

  *real_roots =
      Eigen::Map<VectorXd>(computed_roots.data(), computed_roots.size());

  return true;
}


// An iterative solver to find the closest root based on an initial guess. We
// use Laguerre's method, which is a polynomial root finding method that
// converges to a root with very high certainty. For multiple roots, the
// convergence is linear, otherwise it is cubic.
double FindRealRootIterative(const VectorXd& polynomial,
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
      multiplied_poly.reverse()(i + j) +=
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
    quotient->reverse()(numerator.size() - divisor.size()) =
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

void MinimizePolynomial(const VectorXd& polynomial, const double x_min,
                        const double x_max, double* optimal_x,
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

}  // namespace theia
