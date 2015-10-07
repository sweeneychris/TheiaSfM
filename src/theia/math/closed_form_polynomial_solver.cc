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

#include "theia/math/closed_form_polynomial_solver.h"

#include <glog/logging.h>

#include <complex>
#include <cmath>

namespace theia {
namespace {

// Solve depressed cubic using Cardano's method.
int SolveDepressedCubic(const double p, const double q,
                        std::complex<double>* roots) {
  if (p == 0.0) {
    roots[0] = std::pow(-1.0 * q, 1.0 / 3.0);
    return 1;
  }

  std::complex<double> cubic_root_of_unity(-0.5, 0.5 * std::sqrt(3.0));
  std::complex<double> temp = q * q / 4.0 + p * p * p / 27.0;
  std::complex<double> sqrt_t = std::sqrt(temp);
  std::complex<double> u = std::pow(-0.5 * q + sqrt_t, 1.0 / 3.0);
  std::complex<double> v = std::pow(-0.5 * q - sqrt_t, 1.0 / 3.0);
  roots[0] = u + v;
  roots[1] =
      u * cubic_root_of_unity + v * cubic_root_of_unity * cubic_root_of_unity;
  roots[2] =
      u * cubic_root_of_unity * cubic_root_of_unity + v * cubic_root_of_unity;
  return 3;
}
}  // namespace

// Provides solutions to the equation a*x^2 + b*x + c = 0.
int SolveQuadraticReals(const double a, const double b, const double c,
                        double* roots) {
  std::complex<double> complex_roots[2];
  int num_complex_solutions = SolveQuadratic(a, b, c, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    roots[num_real_solutions++] = complex_roots[i].real();
  }
  return num_real_solutions;
}

int SolveQuadraticReals(const double a, const double b, const double c,
                        const double tolerance, double* roots) {
  std::complex<double> complex_roots[2];
  int num_complex_solutions = SolveQuadratic(a, b, c, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    if (std::abs(complex_roots[i].imag()) < tolerance) {
      roots[num_real_solutions++] = complex_roots[i].real();
    }
  }
  return num_real_solutions;
}

int SolveQuadratic(const double a, const double b, const double c,
                   std::complex<double>* roots) {
  // If the equation is actually linear.
  if (a == 0.0) {
    roots[0] = -1.0 * c / b;
    return 1;
  }

  const double D = b * b - 4 * a * c;
  const double sqrt_D = std::sqrt(std::abs(D));

  // Real roots.
  if (D >= 0) {
      // Stable quadratic roots according to BKP Horn.
    // http://people.csail.mit.edu/bkph/articles/Quadratics.pdf
    if (b >= 0) {
      roots[0] = (-b - sqrt_D) / (2.0 * a);
      roots[1] = (2.0 * c) / (-b - sqrt_D);
    } else {
      roots[0] = (2.0 * c) / (-b + sqrt_D);
      roots[1] = (-b + sqrt_D) / (2.0 * a);
    }
    return 2;
  }

  // Use the normal quadratic formula for the complex case.
  roots[0].real(-b / (2.0 * a));
  roots[1].real(-b / (2.0 * a));
  roots[0].imag(sqrt_D / (2.0 * a));
  roots[1].imag(-sqrt_D / (2.0 * a));
  return 2;
}

// Provides solutions to the equation a*x^3 + b*x^2 + c*x + d = 0 using Cardan's
// method.
int SolveCubicReals(const double a, const double b, const double c,
                    const double d, double* roots) {
  std::complex<double> complex_roots[3];
  int num_complex_solutions = SolveCubic(a, b, c, d, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    roots[num_real_solutions++] = complex_roots[i].real();
  }
  return num_real_solutions;
}

int SolveCubicReals(const double a, const double b, const double c,
                    const double d, const double tolerance, double* roots) {
  std::complex<double> complex_roots[3];
  int num_complex_solutions = SolveCubic(a, b, c, d, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    if (std::abs(complex_roots[i].imag()) < tolerance) {
      roots[num_real_solutions++] = complex_roots[i].real();
    }
  }
  return num_real_solutions;
}

int SolveCubic(const double a, const double b, const double c, const double d,
               std::complex<double>* roots) {
  if (a == 0.0) {
    return SolveQuadratic(b, c, d, roots);
  }

  // Solve by first reducing the problem to a depressed cubic.
  double p = (3.0 * a * c - b * b) / (3.0 * a * a);
  double q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) /
             (27.0 * a * a * a);
  int num_solutions = SolveDepressedCubic(p, q, roots);
  // Transform solution back to normal params.
  roots[0] -= b / (3.0 * a);
  roots[1] -= b / (3.0 * a);
  roots[2] -= b / (3.0 * a);
  return num_solutions;
}

// Provides solutions to the equation a*x^4 + b*x^3 + c*x^2 + d*x + e = 0 using
// Ferrari's method to reduce to problem to a depressed cubic.
int SolveQuarticReals(const long double a, const long double b,
                      const long double c, const long double d,
                      const long double e, long double* roots) {
  std::complex<long double> complex_roots[4];
  int num_complex_solutions = SolveQuartic(a, b, c, d, e, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    roots[num_real_solutions++] = complex_roots[i].real();
  }
  return num_real_solutions;
}

int SolveQuarticReals(const long double a, const long double b,
                      const long double c, const long double d,
                      const long double e, const long double tolerance,
                      long double* roots) {
  std::complex<long double> complex_roots[4];
  int num_complex_solutions = SolveQuartic(a, b, c, d, e, complex_roots);
  int num_real_solutions = 0;
  for (int i = 0; i < num_complex_solutions; i++) {
    if (std::abs(complex_roots[i].imag()) < tolerance) {
      roots[num_real_solutions++] = complex_roots[i].real();
    }
  }
  return num_real_solutions;
}

int SolveQuartic(const long double a, const long double b, const long double c,
                 const long double d, const long double e,
                 std::complex<long double>* roots) {
  const long double a_pw2 = a * a;
  const long double b_pw2 = b * b;
  const long double a_pw3 = a_pw2 * a;
  const long double b_pw3 = b_pw2 * b;
  const long double a_pw4 = a_pw3 * a;
  const long double b_pw4 = b_pw3 * b;

  const long double alpha = -3.0l * b_pw2 / (8.0l * a_pw2) + c / a;
  const long double beta =
      b_pw3 / (8.0l * a_pw3) - b * c / (2.0l * a_pw2) + d / a;
  const long double gamma =
      -3.0l * b_pw4 / (256.0l * a_pw4) + b_pw2 * c / (16.0l * a_pw3) -
      b * d / (4.0l * a_pw2) + e / a;

  const long double alpha_pw2 = alpha * alpha;
  const long double alpha_pw3 = alpha_pw2 * alpha;

  const std::complex<long double> P(-alpha_pw2 / 12.0l - gamma, 0);
  const std::complex<long double> Q(
      -alpha_pw3 / 108.0l + alpha * gamma / 3.0l - std::pow(beta, 2.0l) / 8.0l,
      0);
  const std::complex<long double> R =
      -Q / 2.0l +
      std::sqrt(std::pow(Q, 2.0) / 4.0l + std::pow(P, 3.0l) / 27.0l);

  const std::complex<long double> U = std::pow(R, (1.0l / 3.0l));
  std::complex<long double> y;

  const long double kEpsilon = 1e-8;
  if (std::abs(U.real()) < kEpsilon) {
    y = -5.0l * alpha / 6.0l - std::pow(Q, (1.0l / 3.0l));
  } else {
    y = -5.0l * alpha / 6.0l - P / (3.0l * U) + U;
  }

  const std::complex<long double> w = std::sqrt(alpha + 2.0l * y);

  roots[0] =
      -b / (4.0l * a) +
      0.5l * (w + std::sqrt(-(3.0l * alpha + 2.0l * y + 2.0l * beta / w)));
  roots[1] =
      -b / (4.0l * a) +
      0.5l * (w - std::sqrt(-(3.0l * alpha + 2.0l * y + 2.0l * beta / w)));
  roots[2] =
      -b / (4.0l * a) +
      0.5l * (-w + std::sqrt(-(3.0l * alpha + 2.0l * y - 2.0l * beta / w)));
  roots[3] =
      -b / (4.0l * a) +
      0.5l * (-w - std::sqrt(-(3.0l * alpha + 2.0l * y - 2.0l * beta / w)));

  return 4;
}
}  // namespace theia
